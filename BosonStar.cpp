#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#ifndef BOSONSTAR_CPP_
#define BOSONSTAR_CPP_

#include "BosonStar.h"
#include "DimensionMacros.h"
#include "mathutils.h"
#include <sstream>
#include<iomanip>

using namespace std;


//overloads for addition/subtraction/scalar multiplication of fieldstate vars, useful for simplifying RK4 evolution code
FieldState operator+(const FieldState& s1, const FieldState& s2)
{
    return (FieldState){s1.A + s2.A, s1.X + s2.X, s1.phi + s2.phi, s1.eta + s2.eta};
}

FieldState operator-(const FieldState& s1, const FieldState& s2)
{
    return (FieldState){s1.A - s2.A, s1.X - s2.X, s1.phi - s2.phi, s1.eta - s2.eta};
}

FieldState operator*(double c, const FieldState& s)
{
    return (FieldState){c * s.A, c * s.X, c * s.phi, c * s.eta};
}

//potential and its  derivative (wrt |phi|^2)
double BosonStar::V( const double A)
{
    if (!solitonic)
        return mu * mu * A * A + lambda * pow(A,4);

    else
        return mu * mu * A * A * pow((1. - 2. * pow(A / sigma, 2)), 2) + lambda * pow(A,4);

}

double BosonStar::dV( const double A)
{
    if (!solitonic)
        return mu * mu + 2 * A * A * lambda;

    else
        return mu * mu - 8. * mu * mu * pow(A / sigma, 2) + 12. * mu * mu * pow(A / sigma, 4) + 2 * A * A * lambda ;

}

//read parameters in from file BSParams.par. TODO: default values
void BosonStar::read_parameters(bool quiet)
{
    ifstream params{ "BSParams.par" };

    // Print an error and exit if file cannot open
    if (!params)
    {
        std::cerr << "Could not open BSParams.par\n";
        abort();
    }

    string current_line{};

    while (getline(params, current_line))
    {
        fill_parameter(current_line, "mu = ", mu, quiet);
        fill_parameter(current_line, "lambda = ", lambda, quiet);
        fill_parameter(current_line, "solitonic = ", solitonic, quiet);
        fill_parameter(current_line, "G = ", G, quiet);
        fill_parameter(current_line, "sigma = ", sigma, quiet);
        fill_parameter(current_line, "eigen = ", eigen, quiet);
        fill_parameter(current_line, "alpha_central = ", alpha_central, quiet);
        fill_parameter(current_line, "A_central = ", A_central, quiet);
        fill_parameter(current_line, "frequency_guess = ", frequency_guess, quiet);
        fill_parameter(current_line, "freq_epsilon = ", freq_epsilon, quiet);
        fill_parameter(current_line, "isotropic = ", isotropic, quiet);
        fill_parameter(current_line, "r_match_fac = ", r_match_fac, quiet);

        fill_parameter(current_line, "R = ", R, quiet);
        fill_parameter(current_line, "courant_factor = ", courant_factor, quiet);
        fill_parameter(current_line, "n_gridpoints = ", n_gridpoints, quiet);
        fill_parameter(current_line, "stop_time = ", stop_time, quiet);

    }
}

//right-hand side of the EKG system of equations
//IDEA: add term in alpha to eta (tentative... for now)
FieldState BosonStar::state_RHS(const double radius, const long double frequency, FieldState  s, bool asymptotic_region, bool given_A)
{
    //enforce minimum radius epsilon if needed
    double epsilon = 0.0000001;
    double r = ((radius == 0.) ? epsilon : radius);

    //helper terms that appear multiple times / should be zeroed at r = 0
    double T1 = 0.5 * (s.X * s.X - 1) / r;
    double T2 =  s.eta * s.eta + frequency * frequency * s.A * s.A / exp(2 * s.phi) ;
    double T3 = s.eta / r;

    double F1 = 2. * PI * G * r * s.X * s.X;

    //zero out terms that should be zeroed at origin as eta = 0, X = 1 there
    if (radius == 0.)
        {T1 = 0; T3 = 0;}

    double dPhi = T1 + F1 * (T2 - V(s.A));

    //in the asymptotic region, evolve phi and X normally but do not update A, eta (these will be hard-coded to asymptotic expressions)
    if (asymptotic_region)
    {
        return  (FieldState) {0., s.X * ( F1 * (T2 + V(s.A)) - T1 ), dPhi, 0.};
    }
    else if (given_A)
    {
        //return RHS of field state variables outside of asymptotic region, A excepted
        return  (FieldState) {0, s.X * ( F1 * (T2 + V(s.A)) - T1 ), dPhi,
        -2. * T3 - s.eta * dPhi + s.X * s.A * (dV(s.A) - frequency * frequency  / exp(2 * s.phi))};
    }

    else
    {
        //return RHS of field state variables outside of asymptotic region
        return  (FieldState) {s.X * s.eta, s.X * ( F1 * (T2 + V(s.A)) - T1 ), dPhi,
        -2. * T3 - s.eta * dPhi + s.X * s.A * (dV(s.A) - frequency * frequency  / exp(2 * s.phi))};
    }

   //test for playing with RK4 convergence-- seems that 1/r terms can spoil 4th-order convergence (but we still get tight 3rd-order)
   //return  (FieldState) {s.A / r + s.X, s.X + s.phi, s.phi + s.eta, s.eta + s.A};

}

void BosonStar::rk4_solve (const long double freq)
{

    state = {(FieldState){A_central, 1.0, log(alpha_central), 0.0 }};
    radius_array = {0.};

    blowup_point = n_gridpoints; //make this 1 larger than max possible value to start, in case solution does not break

    state.resize(n_gridpoints);
    radius_array.resize(n_gridpoints);

    //cout << " \nInitialized state "  << endl;

    double dr = R / (n_gridpoints - 1);
    //double dt = dr * courant_factor;

    //inter-level state values for RK4 evolution
    FieldState s1, s2, s3, s4;

    //fill in grid using RK4 evolution
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        s1 = state_RHS(r, freq, state[j], 0);
        s2 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s1, 0);
        s3 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s2, 0);
        s4 = state_RHS(r + dr, freq, state[j] + dr * s3, 0);

        //update state variables and radius array
        state[j + 1] = state[j] + (dr / 6.) * (s1 + 2 * s2 + 2 * s3 + s4);
        radius_array[j+1] = (j + 1) * dr;

        //cout << "A = " << state[j].A << ", X = " << state[j].X << ", phi = " << state[j].phi << ", eta = " << state[j].A << ", m = " <<  r  / 2. * (1 - (1 / (state[j].X * state[j].X)))<< endl;

        if (isnan(state[j].A) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].eta))
        {
            //cerr << "State values have become nan on step " << j << endl;
            blowup_point = j;
            return;
        }
    }

    //cout << "\n" << "Finished RK4 evolution" << endl;

}

double BosonStar::m(int j)
{
    if (j >= n_gridpoints || j < 0)
        {
            cerr << "ERROR: invalid index passed to m(r)" << endl;
            exit(1);
        }
    return radius_array[j]  / 2. * (1 - (1 / (state[j].X * state[j].X)));
}

//writes field values to BSdata.dat
void BosonStar::write_field()
{

    std::ofstream data_file{ "BSdata.dat" };
    //double dr = R / (n_gridpoints - 1);

    // If we couldn't open the output file stream for writing
    if (!data_file)
    {
        // Print an error and exit
        std::cerr << "BSdata.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 0; j < n_gridpoints; j++)
    {
        data_file << std::setprecision (10) << radius_array[j] << "   " << state[j].A << "    " << state[j].X << "    " << state[j].phi << "    " << state[j].eta << "    " << m(j) << endl;
    }


}

void BosonStar::double_resolution()
{
    n_gridpoints = 2 * n_gridpoints - 1;
}

//convergence test for the RK4 solver
void BosonStar::convergence_test(long double freq)
{
    //FieldState vectors for low, medium, high resolution, each off by factor of 2
    std::vector<FieldState> state_l, state_m, state_h;

    if (freq == 0.) freq = omega_pre_rescale; //if no frequency provided (for rk4 solver) use omega

    //fill in state vectors at low, mid, high resolution
    solve(1); // rk4_solve(freq);

    cout << setprecision(10) << static_cast<double>(omega) << endl;

    state_l = state;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1); //rk4_solve(freq);
    state_m = state;

    cout <<  setprecision(10) << static_cast<double>(omega) << endl;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1); //rk4_solve(freq);
    state_h = state;

    cout << setprecision(10) << static_cast<double>(omega) << endl;

    //write to file
    std::ofstream conv_file{ "conv.dat" };

    // If we couldn't open the output file stream for writing
    if (!conv_file)
    {
        // Print an error and exit
        std::cerr << "Error: conv.dat could not be opened for writing\n";
        exit(1);
    }

    //write change in A between med/low and high/med resolution to file.
    for (int j = 0; j < (n_gridpoints +3) / 4; j++)
    {
        conv_file << radius_array[4*j] << "   " << state_m[2*j].A -  state_l[j].A << "    " << 8.*(state_h[4*j].A -  state_m[2*j].A) << endl;
    }
}

//returns number of zero crossings in A for BS solution. Must call rk4_solve() first! Could optimize by doing during RK4 if needed?
int BosonStar::count_zero_crossings()
{
    int zero_crossings = 0;

    //check for zero crossings up to blowup point (neglecting one extra point as we sometimes get spurious zero crossings on explosion)
    for (int j = 1; j < blowup_point - 1; j++)
    {
        if ((state[j].A == 0.) || (state[j].A > 0. && state[j - 1].A < 0.) || (state[j].A < 0. && state[j - 1].A > 0.) )
            {zero_crossings++;}
    }

    //cout << "\n" << zero_crossings << " zero crossings found" << endl;

    return zero_crossings;
}

//interval bisection algorithm-- initial_guess must be greater than true value for search to work
long double BosonStar::find_frequency(bool quiet)
{
    //lower bound on frequency
    long double lower_guess = frequency_guess;

    //tolerance for uncertainty in frequency, as well as minimum frequency to try
    long double epsilon = freq_epsilon;

    //at first, compute BS models, halving frequency each time until we find one with fewer/equal zero crossings to the number desired
    rk4_solve(lower_guess);
    if (count_zero_crossings() <= eigen && !quiet)
    {
        cerr << "WARNING: initial guess of "<< frequency_guess << " too small for A_central = " << A_central << endl;
        //exit(1);
    }

    while (count_zero_crossings() > eigen && lower_guess > epsilon)
    {
        lower_guess /= 2;
        rk4_solve(lower_guess);
    }

    if (lower_guess <= epsilon)
    {
        cerr << "ERROR: could not find suitable lower frequency bound for A_central = " << A_central << endl;
        //exit(1);
    }

    //upper bound on frequency
    long double upper_guess = lower_guess * 2.;
    long double midpoint = 0.5 * (upper_guess + lower_guess);

    int num_solves = 0;

    //now use interval bisection to converge to correct frequency
    while ( (upper_guess - lower_guess) > epsilon && num_solves < 100)
    {
        //replace upper/lower bound with midpoint if it has greater/as many or fewer zero crossings than desired, so both ultimately converge on boundary value between eigen and eigen + 1 zero crossings as needed
        midpoint = 0.5 * (upper_guess + lower_guess);
        rk4_solve(midpoint);
        if (count_zero_crossings() > eigen)
            upper_guess = midpoint;
        else
            lower_guess = midpoint;

        num_solves++;

       // cout << static_cast<double>(upper_guess - lower_guess) << endl;
    }

    omega = lower_guess; //use lower_guess for frequency so we are guarenteed to have right # zero crossings in principle...
    rk4_solve(omega);

    if (count_zero_crossings() != eigen)
        cerr << "WARNING: zero crossing count may be incorrect for A_central = " << A_central << ", found " << count_zero_crossings() << endl;


    if (!quiet) cout << " \nFound solution in eigenstate " << eigen  << " with frequency " << static_cast<double>(omega) << endl;


    return lower_guess;


}

//returns the index of the last minimum of |A|, searching inwards from the blowup point.
int BosonStar::find_last_minimum()
{

    for (int j = blowup_point - 2; j > 1; j--)
    {
     if (abs(state[j].A) <= abs(state[j - 1].A) && abs(state[j].A) <= abs(state[j + 1].A) && state[j].A < 0.05 * A_central  )
        return j;
    }

    //case where our frequency guess is good enough that the solution does not blow up within the grid
    if (blowup_point == n_gridpoints && abs(state[n_gridpoints - 1].A) <= abs(state[n_gridpoints - 2].A ))
        return (n_gridpoints - 1);


    return -1; //signifies that next function should skip
}
//fill up the solution in the region after the blowup point, using asymptotic matching for A and eta and integrating phi and X. Returns 1 if successful
bool BosonStar:: fill_asymptotic(bool quiet)
{

    double min_index = find_last_minimum();

    if (min_index == -1 && !quiet)
        cerr << "ERROR: No local minimum in |A| found for A_central = " << A_central << endl;

    if (min_index == -1)
        return 0;

    //double r_match_fac = 0.75;
    int j_match = round( r_match_fac * min_index); //index at which matching takes place

    //cout<< "Found last minimum at index j= " << j_match << endl;

    double dr = R / (n_gridpoints - 1);
    double r_match = dr * j_match; //matching radius

    double phi_match = state[j_match].phi; //match values in current gauge
    double A_match = state[j_match].A;
    double eta_match = state[j_match].eta;
    double deta_match = (state[j_match].eta - state[j_match - 1].eta) / dr; //estimate for derivative of eta at r_match, used to crudely fit exponential falloff for eta to first order

    //ensures continuity of A
    double A_factor = A_match * exp(sqrt( 1 - pow(omega / exp(phi_match), 2)) * r_match) * r_match;

    //fit exponential of form B * exp(-k * r) to eta in exponential region. If eta is already zero, just keep it there (will likely not happen in practice)
    double k, B;
    if (eta_match != 0.)
    {
        k = abs(deta_match / eta_match);
        B = exp(k * r_match) * eta_match;
    }
    else
    {
        k = 0.; B = 0.;
    }

    FieldState s1, s2, s3, s4;
    for (int j = j_match; j < n_gridpoints - 1; j++)
    {
        long double r = j * dr;

        s1 = state_RHS(r, omega, state[j], 1);
        s2 = state_RHS(r + dr / 2.,  omega, state[j] + 0.5 * dr * s1, 1);
        s3 = state_RHS(r + dr / 2.,  omega, state[j] + 0.5 * dr * s2, 1);
        s4 = state_RHS(r + dr,  omega, state[j] + dr * s3, 1);

        //update state variables and radius array
        state[j + 1] = state[j] + (dr / 6.) * (s1 + 2 * s2 + 2 * s3 + s4);
        radius_array[j + 1] = (j + 1) * dr;

        //cout << A_factor * exp(-sqrt( 1 - pow(omega / exp(phi_match), 2)) * radius_array[j + 1]) / radius_array[j + 1]  << endl;

        //fix asymptotic values for A, eta
        state[j + 1].A = A_factor * exp(-sqrt( 1 - pow(omega / exp(phi_match), 2)) * radius_array[j + 1]) / radius_array[j + 1];
        state[j + 1].eta = - (1 /radius_array[j + 1] + sqrt( 1 - pow(omega / exp(phi_match), 2)) ) * state[j + 1].A;

        //cout << "A = " << state[j].A << ", X = " << state[j].X << ", phi = " << state[j].phi << ", eta = " << state[j].A << ", m = " <<  r  / 2. * (1 - (1 / (state[j].X * state[j].X)))<< endl;

        if (isnan(state[j].A) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].eta))
        {
            cerr << "WARNING: During asymptotic phase, state values became nan on step " << j <<  " with A_central = " << A_central << endl;
            return 0;
        }
    }

    //finally, enforce Schwarzchild condition Phi = - ln(X) at grid boundary to approximate Phi(infty) = 0
    double phi_shift = -log(state[n_gridpoints - 1].X) - state[n_gridpoints - 1].phi;

    //offset all phi values by phi_shift, equivalent to rescaling lapse by exp(phi_shift)
    omega_pre_rescale = omega;
    rescale_lapse (phi_shift);

    //fills total mass
    M = m(n_gridpoints - 1);

    //fills radius
    int j = 0;
    while (j < n_gridpoints - 1 && m(j) < 0.99 * M)
    {
        //just uses linear interpolation to find r_99 for now
        if ( m(j + 1) != m(j) )
            r_99 = dr * j + dr * (0.99 * M - m(j) ) / (m(j + 1) - m(j));

        j++;
    }

     noether_charge = get_noether_charge();
     binding_energy = M - mu * noether_charge;
     compactness = M / r_99;

     if (!quiet)
     {
        cout << "\nFinal frequency after enforcing lapse condition is " << static_cast<double>(omega) << endl;
        cout << "\nMass is M =  " << M << endl;
        cout << "\nBinding energy is E = " << binding_energy << endl;
     }




     return 1;
}

//rescales lapse by e^(phi_shift) as well as the frequency
void BosonStar::rescale_lapse (double phi_shift)
{
     for (FieldState& s: state)
        s.phi += phi_shift;

     omega *= exp(phi_shift);
}

//simply calls find_frequency and fill_asymptotic in one; return 1 iff successful
bool BosonStar::solve(bool quiet)
{
    find_frequency(quiet);
    return fill_asymptotic(quiet);
}

//returns the noether charge associated with the model. Must have computed model + frequency first (polar)
double BosonStar::get_noether_charge()
{

    double Q = 0.; //charge to return
    double dr = R / (n_gridpoints - 1);

    //integrate charge over radius
    for (int j = 0; j < n_gridpoints; j++)
    {
        double J_0 = 2 * omega * state[j].A * state[j].A; //0th component of the conserved current covector
        double r = j * dr;
        Q += r * r * dr * state[j].X * J_0 * exp(-1. * state[j].phi);
    }

    Q *= 2 * M_PI;
    return Q;
}

//returns right-hand side for the ODE that f solves. r is areal radius
double BosonStar::f_RHS(double r, double f)
{
    //indices of gridpoints bounding the desired value of r
    int j_low = floor((r / R) * (n_gridpoints - 1));
    int j_high = ceil((r / R) * (n_gridpoints - 1));

    double dr = R / (n_gridpoints - 1);
    double gap_frac = (r - dr * j_low) / dr; //portion of reached between gridpoints

    double X = (1 - gap_frac) * state[j_low].X + gap_frac * state[j_high].X; //linearly interpolates X; may need better interpolation method!

    int j0 = bound(j_low - 1, 0, n_gridpoints - 4);

    X = cubic_interp(r, state[j0].X, state[j0 + 1].X, state[j0 + 2].X, state[j0 + 3].X, j0, dr); //may need smth better at 0; exploit symmetry!

    double r2 = (r == 0.) ? 10e-10 : r;
    double rhs = f * (X - 1.) / r2;

    return rhs;

}

//returns areal radius corresponding to a given isotropic index
double BosonStar::r_areal(int j_iso)
{

    if (j_iso == 0)
        return 0.;

    double dr = R / (n_gridpoints - 1);
    double r_iso = dr * j_iso;

    int j_areal = 0; //index of the upper bound of r(R)

    while (j_areal < n_gridpoints && r_iso > r_iso_array[j_areal])
        j_areal++;

    //use asymptotic expression at large r (where we'd be out of the areal array range)
    if(j_areal == n_gridpoints)
    {
        //cerr << "WARNING: extrapolating to find areal radius corresponding to isotropic index " << j_iso << endl;
        //exit(1);

        return r_iso + M + M * M / (4. * r_iso);
    }

    //upper and lower bounds of r_areal
    //double r_low = (j_areal - 1) * dr;
    //double r_high = j_areal * dr;

    int j0 = bound(j_areal - 2, -2, n_gridpoints - 4);

    vector<int> j{0, 1, 2, 3};
    vector<double> r(4);
    vector<double> r_areal(4);
    for (int& index : j)
        {r[index] = r_iso_array[j0 + index];
        r_areal[index] = (j0 + index) * dr;
        }

    //use symmetry thru z = 0 at inner boundary
    if (j0 == -1)
        r = {r_iso_array[1], r_iso_array[0], r_iso_array[1], r_iso_array[2]};

    if (j0 == -2)
        r = {r_iso_array[2], r_iso_array[1], r_iso_array[0], r_iso_array[1]};


    return lagrange_interp(r_iso, r, r_areal);



    //interpolate via cubic to find r_areal, using nearest 4 known positions of r_iso (r_areal)
    /*return (r_iso - r[1]) * (r_iso - r[2])  * (r_iso - r[3]) * (j0 + j[0]) * dr / ( (r[0] - r[1]) *  (r[0] - r[2]) * (r[0] - r[3]))
         + (r_iso - r[0]) * (r_iso - r[2])  * (r_iso - r[3]) * (j0 + j[1]) * dr / ( (r[1] - r[0]) *  (r[1] - r[2]) * (r[1] - r[3]))
         + (r_iso - r[0]) * (r_iso - r[1])  * (r_iso - r[3]) * (j0 + j[2]) * dr / ( (r[2] - r[0]) *  (r[2] - r[1]) * (r[2] - r[3]))
         + (r_iso - r[0]) * (r_iso - r[1])  * (r_iso - r[2]) * (j0 + j[3]) * dr / ( (r[3] - r[0]) *  (r[3] - r[1]) * (r[3] - r[2]));*/

}

//fill out arrays for radius, field amplitude, and conformal factor in isotropic coordinates. Must have solved BS model first!
void BosonStar::fill_isotropic_arrays()
{
    double dr = R / (n_gridpoints - 1);

    double f1, f2, f3, f4; //, g1, g2, g3, g4; //intermediate rk4 values

    //local array to store values of f, the ratio between isotropic and areal radii; assumed = 1 at the origin for now and later to be rescaled
    vector<double> f_array(n_gridpoints);
    vector<double> r_areal_array(n_gridpoints);
    f_array[0] = 1.;

    //g is just f on a grid uniform in isotropic rather than areal radius. Used to avoid need to apply grid inversion to find r_areal(r_iso), which adds to error.
    //vector<double> g_array(n_gridpoints);
    //g_array[0] = 1.;

    r_iso_array.resize(n_gridpoints);
    psi_iso_array.resize(n_gridpoints);
    phi_iso_array.resize(n_gridpoints);
    A_iso_array.resize(n_gridpoints);

    r_iso_array[0] = 0;

    //evolve f using RK4 evolution and fill out R_iso array
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        f1 = f_RHS(r, f_array[j]);
        f2 = f_RHS(r + dr / 2., f_array[j] + 0.5 * dr * f1);
        f3 = f_RHS(r + dr / 2.,  f_array[j] + 0.5 * dr * f2);
        f4 = f_RHS(r + dr, f_array[j] + dr * f3);

        //update f and isotropic radius array
        f_array[j + 1] = f_array[j] + (dr / 6.) * (f1 + 2 * f2 + 2 * f3 + f4);
        r_iso_array[j + 1] = f_array[j + 1] * (r + dr);

    }

    //rescale factors to match R_iso to asymptotic Schwarzchild at boundary
    double iso_scale_factor = (R - M - M *  M / (4. * ( R -  M))) / r_iso_array[n_gridpoints - 1];

    //rescale entire isotropic radius array and f array
    for (double& r_iso: r_iso_array)
        r_iso *= iso_scale_factor;

    for (double& f: f_array)
        f *= iso_scale_factor;

    //fill in A, phi, psi arrays
    for (int j = 0; j < n_gridpoints; j++)
    {
        r_areal_array[j] = r_areal(j);

        //upper/lower index bounds on areal radius
        int j_areal_low = floor(r_areal_array[j] / dr);
        int j_areal_high = ceil (r_areal_array[j] / dr);

        //double gap_frac = (r_areal_array[j] - dr * j_areal_low) / dr; //portion between gap between gridpoints

        double f;

        //linearly interpolate from the radial array data where we won't go OOB in radial array
        if (j_areal_high <= n_gridpoints - 1)
        {
            //phi_iso_array[j]  = (1. - gap_frac) * state[j_areal_low].phi + gap_frac * state[j_areal_high].phi;
            //A_iso_array[j]  = (1. - gap_frac) * state[j_areal_low].A + gap_frac * state[j_areal_high].A;
            //f = (1. - gap_frac) * f_array[j_areal_low] + gap_frac * f_array[j_areal_high];

            //cubic spline interpolation for all but boundaries of domain
            int j0 = bound(j_areal_low - 1, 0, n_gridpoints - 4);

            /*double r = r_areal_array[j];

            vector<double> phi_vals(4);
            vector<double> A_vals(4);
            vector<double> f_vals(4);*/


            phi_iso_array[j] = cubic_interp(r_areal_array[j], state[j0].phi, state[j0 + 1].phi, state[j0 + 2].phi, state[j0 + 3].phi, j0, dr );
            A_iso_array[j] = cubic_interp(r_areal_array[j], state[j0].A, state[j0 + 1].A, state[j0 + 2].A, state[j0 + 3].A, j0, dr );
            f = cubic_interp(r_areal_array[j], f_array[j0], f_array[j0 + 1], f_array[j0 + 2], f_array[j0 + 3], j0, dr );
        }


        else //extrapolation case; will kick in near boundary
        {
            phi_iso_array[j] = state[n_gridpoints - 1].phi  + (r_areal_array[j] - R) * (state[n_gridpoints - 1].phi - state[n_gridpoints - 2].phi ) / dr;
            A_iso_array[j] = state[n_gridpoints - 1].A  + (r_areal_array[j] - R) * (state[n_gridpoints - 1].A - state[n_gridpoints - 2].A ) / dr;
            f = f_array[n_gridpoints - 1]  + (r_areal_array[j] - R) * (f_array[n_gridpoints - 1] - f_array[n_gridpoints - 2] ) / dr;
            //f = cubic_interp(r_areal_array[j], f_array[n_gridpoints - 4], f_array[n_gridpoints - 3], f_array[n_gridpoints - 2], f_array[n_gridpoints - 1], n_gridpoints - 4, dr );
            //f = pow((1 + M / (2 * dr * j)), -2.);
        }

        //linearly transition between interpolated and asymptotic f in this region. u_frac must be small enough that we have non-extrapolated r_areal values there
        double l_frac = 0.8;
        double u_frac = 0.9;

        if (j > n_gridpoints * l_frac)
        {
            double portion_crossed = min((j - n_gridpoints * l_frac) / (n_gridpoints * (u_frac - l_frac)), 1.);

            f = f * (1 - portion_crossed) + (pow((1 + M / (2. * dr * j)), -2.)) * portion_crossed;
        }
    //psi is simply 1 / sqrt(f)
    psi_iso_array[j] = 1. / sqrt(/*g_array[j]*/f);
    }
}

//writes isotropic values to isotropic.dat
void BosonStar::write_isotropic()
{

    std::ofstream data_file{ "isotropic.dat" };
    double dr = R / (n_gridpoints - 1);

    // If we couldn't open the output file stream for writing
    if (!data_file)
    {
        // Print an error and exit
        std::cerr << "isotropic.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 1; j < n_gridpoints - 1; j++)
    {
        data_file << std::setprecision (10) << radius_array[j] << "   " << r_iso_array[j] << "    " << r_areal(j) << "    " << psi_iso_array[j]  <<  "    " << A_iso_array[j]  <<  "    " << (-2. * psi_iso_array[j] + psi_iso_array[j - 1] + psi_iso_array[j + 1] ) / (dr * dr) << endl;
    }
}

//attempts to construct solution corresponding to a given frequency by finding A.
//currently assumes omega decreasing with A
//tried to use for cycle_models; currently defunct
bool BosonStar::solve_finding_A(long double freq, double A_guess, double A_range, bool quiet)
{
    long double A_upper = A_guess + A_range;
    long double A_lower = A_guess - A_range;

    /*int max_blowup_point = 0;

    long double A_best = A_guess;
    for (int n = 0; n < 1000; ++n)
    {
        int sign = 2 * (n % 2) - 1;
        A_central = A_guess + sign * A_range * n;
        rk4_solve(freq);

        if (blowup_point > max_blowup_point)
        {
            max_blowup_point = blowup_point;
            A_best = A_central;
            cout << "blowup point = " << blowup_point << " at A = " << A_best<< endl;;
        }

    }*/


    //interval bisection to find A, assuming omega is decreasing with A. may need to first figure out sign of d_A(omega)/d_omega
    while (A_upper - A_lower > freq_epsilon)
    {
        long double midpoint = 0.5 * (A_upper + A_lower);


        rk4_solve(midpoint);
        if (count_zero_crossings() > eigen)
            A_upper = midpoint;
        else
            A_lower = midpoint;

        //cout << count_zero_crossings() << endl;
    }

    A_central = A_upper;
    rk4_solve(freq);
    fill_asymptotic(quiet);

    return 1;

}


void BosonStar::fill_given_A(const long double freq)
{


    state.resize(n_gridpoints);
    radius_array.resize(n_gridpoints);

    //copy old A-values to use later in RK4 loop
    vector<double> A_vals(n_gridpoints);
    vector<double> eta_vals(n_gridpoints);

    for (int j = 0; j < n_gridpoints; j++)
        {A_vals[j] = state[j].A; eta_vals[j] = state[j].eta;}

    blowup_point = n_gridpoints; //make this 1 larger than max possible value to start, in case solution does not break

    double dr = R / (n_gridpoints - 1);

    //inter-level state values for RK4 evolution
    FieldState s1, s2, s3, s4;

    //fill in grid using RK4 evolution
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        s1 = state_RHS(r, freq, state[j], 1, 1);

        state[j].A = 0.5 * (A_vals[j + 1] + A_vals[j]); //use linearly interpolated A-values for intermediate stages
        state[j].eta = 0.5 * (eta_vals[j + 1] + eta_vals[j]);

        s2 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s1, 1, 1);
        s3 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s2, 1, 1);

        state[j].A = A_vals[j + 1];
        state[j].eta = eta_vals[j + 1];

        s4 = state_RHS(r + dr, freq, state[j] + dr * s3, 1, 1);

        //update state variables and radius array
        state[j + 1] = state[j] + (dr / 6.) * (s1 + 2 * s2 + 2 * s3 + s4);

        state[j + 1].A = A_vals[j];
        state[j + 1].eta = eta_vals[j];

        radius_array[j+1] = (j + 1) * dr;


        if (isnan(state[j].A) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].eta))
        {
            //cerr << "State values have become nan on step " << j << endl;
            blowup_point = j;
            cout << "WARNING: obtained nan's during fill_given_A routine." << endl;
            return;
        }
    }

    state[0].A = A_vals[0];
    state[0].eta = eta_vals[0];

    //cout << "\n" << "Finished RK4 evolution" << endl;
}


//reads in data files from Uli's thin shell code
void BosonStar::read_thinshell()
{
    ifstream A_file("A.dat");
    ifstream X_file("rX.dat");
    ifstream m_file("m.dat");
    ifstream eta_file("thet.dat");
    ifstream phi_file("Phi.dat");
    ifstream info_file("output.dat");

    solitonic = 1;

    if (!A_file.is_open() ||!X_file.is_open() ||!phi_file.is_open() ||!info_file.is_open() || !m_file.is_open() || !eta_file.is_open() )
    {
        cerr << "Error reading thinshell files!" << endl;
        exit(1);
    }

    int line_count = 0;
    string line1, lineX, lineA, linePhi, lineInfo, lineThet, linem;
    while (std::getline(A_file, line1)) {
        ++line_count;
    }
    cout << line_count << endl;
    A_file.clear();
    A_file.seekg(0, ios::beg);

    //n_gridpoints = line_count;
    state.resize(n_gridpoints);
    radius_array.resize(line_count);

    int j = 0;
   /* while (std::getline(X_file, lineX))
    {
        double junk;
        std::istringstream iss(lineX);
        iss >> radius_array[j] >> junk;
           j++; //{  if(j > 0) cout << radius_array[j] - radius_array[j - 1] << endl; j++;}

    }*/

    //R = radius_array[line_count - 1];


    radius_array.resize(n_gridpoints);

    //holds the A, phi, and r values from the nonuniform r-y hybrid grid. Needed for interpolation
    //TODO: size properly, currently overestimated
    vector<double> A_vals(line_count);
    vector<double> phi_vals(line_count);
    vector<double> thet_vals(line_count);
    vector<double> m_vals(line_count);
    vector<double> r_vals(line_count);


    j = 0;

    while (std::getline(A_file, lineA))
    {
        std::istringstream iss(lineA);
        if(iss >> r_vals[j] >> A_vals[j])
            j++;
    }


   j = 0;

    while (std::getline(eta_file, lineThet))
    {
        std::istringstream iss(lineThet);
        if(iss >> r_vals[j] >> thet_vals[j])
            j++;
    }

    j = 0;

    while (std::getline(m_file, linem))
    {
        std::istringstream iss(linem);
        if(iss >> r_vals[j] >> m_vals[j])
            j++;
    }

    j = 0;

    double phi_offset;

    //read in info values including omega, M, r_99, C
    while (std::getline(info_file, lineInfo))
    {
        std::istringstream iss(lineInfo);
        if (!lineInfo.empty() && lineInfo[0] != '#') //ignore leading commented lines
        {
            double junk;
            if (!(iss >> A_central >> omega_pre_rescale >> junk >> omega >> phi_offset >> M >> junk >> junk >> compactness >> r_99 >> junk >> junk >> junk ))
                cout << "WARNING: reading thinshell output.dat may have failed " << endl;
        }

    }

    //read in phi-values, accounting for offset and (crudely) computing eta while we're at it
    while (std::getline(phi_file, linePhi))
    {
        std::istringstream iss(linePhi);
        if(iss >> r_vals[j] >> phi_vals[j] )
        {
            phi_vals[j] -= phi_offset; //enforce phi(infty) = 0
            //state[j].eta = (j == 0) ? 0 : (A_vals[j] - A_vals[j - 1] / ( r_vals[j] - r_vals[j - 1] ) / state[j].X ); //wrong currently bc of X
            j++;

            //cout << r_vals[j] << endl;
        }
    }



    double dr = R / (n_gridpoints - 1);
    //interpolate to fill uniform grid with A, phi, eta
    for (int k = 0; k < n_gridpoints; k++)
    {

        int l = 0;

        while (r_vals[l] < k * dr) //finds nearest r-value to uniform grid point
             l++;

    //interpolation order; must be odd!
    int interp_order = 3;

    //if (k * dr > 1.5 * r_99 ) interp_order = 1; //switch to linear interpolation well outside BS; turns out to be bad :(

    //extracts 4 surrounding grid points and their associated A, phi values
    int l0 = bound(l - interp_order + 1, 0, n_gridpoints - interp_order - 1);
    vector<int> L(interp_order + 1);
    for (int i = 0; i < interp_order + 1; i++ )
        L[i] = i;

    vector<double> r(interp_order + 1);
    vector<double> A(interp_order + 1);
    vector<double> m(interp_order + 1);
    vector<double> phi(interp_order + 1);
    vector<double> thet(interp_order + 1);

    for (int& index : L)
    {
        r[index] = r_vals[l0 + index];
        phi[index] = phi_vals[l0 + index];
        m[index] = m_vals[l0 + index];
        A[index] = A_vals[l0 + index];
        thet[index] = thet_vals[l0 + index];
    }

    double mass_val = lagrange_interp(k * dr, r, m);
    state[k].phi = lagrange_interp(k * dr, r, phi);

    //if (k * dr > 1.2 * r_99) //hard cutoff on M, phi prevents spurious oscillations when changing resolution at large r.
   // {
      //  mass_val = M;
        //state[k].phi = log(sqrt(1 - 2 * M / (k * dr)));
    //}

    state[k].A = lagrange_interp(k * dr, r, A);

    state[k].eta = lagrange_interp(k * dr, r, thet) * exp(-state[k].phi - phi_offset ) ;

    state[k].X = 1 /sqrt(1 - 2 * mass_val / (dr * k));
    radius_array[k] = dr * k;

   // double r_match = 28.9;
   // double interval = 0.5;

   // if (dr * k < r_match + 0.5 * interval && dr * k > r_match - 0.5 * interval)
        //{

       // }

    }
    state[0].X = 1; //interpolation for m seems to fail at r = 0

    //fill_given_A(omega);

    //for (int j = 0; j < n_gridpoints; j++)
       // state[j].eta *= exp(state[j].phi + phi_offset );

}

//TODO: resolving may need to only take place within [0,r_99]; add soft AMR support
void BosonStar::cycle_models(int n_stars, double A_0, double delta_A)
{
    read_parameters(0);
    ofstream data_file{"BosonStars.dat"};
    omega_pre_rescale = frequency_guess;

    //initial # of gridpoints (we'll increase for larger models)
    int N_gridpoints_init = n_gridpoints;

    //double refinement whenever we pass a refine_threshold
    int threshold_counter = 0;

    // mini BS: {0.3, 0.375, 0.425, 0.475, 0.525, 0.575, 0.62, 0.66, 0.7, 0.73} is good with n_gridpoints = 2000 from patams file
    //sigma = 0.2: {0.25,0.325,0.375, 0.425, 0.475, 0.525, 0.575};
    //sigma = 0.1: {0.05, 0.15, 0.2, 0.25,0.325,0.375, 0.425, 0.475, 0.525, 0.575};
    vector<double> refine_thresholds{0.05, 0.2, 0.25,0.325,0.375, 0.425, 0.475, 0.525, 0.575};
    //bool passed_last_threshold; //set to 1 after last threshold reached

    //frequency from previous guess for use in update
    long double omega_prev = frequency_guess;

    //holds all unrescaled omega / A / phi[0] values yet computed-
    vector<long double> omega_values(n_stars);
    vector<long double> A_values(n_stars);
    vector<long double> phi0_values(n_stars);

    long double phi0_prev = 0.;

    bool var_dA_method = 0;

    double guess_buffer = 0.05;

    double trouble_start = 10.0991;
    double trouble_end = 10.1004;

    int guess_layer = 0;
    const int max_guess_layer = 50;

    double f_start = frequency_guess;
    double f_search_interval = f_start / max_guess_layer;
    double dA_factor = 1;

    double alpha_start = -1.0;
    double alpha_end = -1.0;

    bool alpha_flag = 0;


    A_central = A_0;

    //double n_gridpoints_no_error = n_gridpoints; //stores number of gridpoints with exponential correction keeping the fractional part to avoid roundoff error accumulating

    for (int j = 0; j < n_stars; j++)
    {

        if (!var_dA_method)
            A_central +=  delta_A * dA_factor;

        if (j!= 0)
            {

                if(j > 2 && !alpha_flag)
                {
                    //frequency_guess = 1.2;//omega_values[j - 1] +  (omega_values[j - 1] - omega_values[j - 2]) + guess_buffer;

                    //this breaks with __float128 : (
                    alpha_central =  exp(2 * phi0_values[j - 1] - phi0_values[j - 2]);
                }

                alpha_central = exp(2 * state[0].phi - phi0_prev);
                omega_prev = omega;
                phi0_prev = state[0].phi;

            }

        if (j > 3 && var_dA_method)
        {


            double theta = atan(10 * (omega_values[j - 1] - omega_values[j - 2]) / (A_values[j - 1] - A_values[j - 2]));
            double theta_prev = atan(10 * (omega_values[j - 2] - omega_values[j - 3]) / (A_values[j - 2] - A_values[j - 3]));

            //account for branch cut in atan
            if (theta > 0 && theta_prev < 0)
                theta_prev += M_PI;

            double theta_new = theta + delta_A * (theta - theta_prev);

            cout << theta_new << endl;

            //cout << theta << endl;

            //linear extrapolation works fine for angles not too steep
            if (/*abs(theta) < 31. * M_PI / 32.*/ j > 20 )
                {
                    A_central = A_values[j - 1] + delta_A * cos(theta_new);
                    frequency_guess = omega_values[j - 1] +  sin(theta_new) * delta_A * 10. + guess_buffer ;
                }

            else
            {
                //cout << "Switched interp method" << endl;

                double A_central_extrap = A_values[j - 3] - 3 * A_values[j - 2] + 3 * A_values[j - 1];
                double omega_central_extrap = omega_values[j - 3] - 3 * omega_values[j - 2] + 3 * omega_values[j - 1];
                double extrap_mag = sqrt( pow(A_central_extrap - A_values[j - 1], 2) + pow(omega_central_extrap - omega_values[j - 1], 2) );

                A_central = A_values[j - 1]  + delta_A * (A_central_extrap - A_values[j - 1]) / extrap_mag;
                frequency_guess = omega_values[j - 1]  + delta_A * (omega_central_extrap - omega_values[j - 1]) / extrap_mag + guess_buffer;
            }
        }

       //if (j == 10)
            //var_dA_method = 1;


        if (A_central > trouble_start && alpha_start < 0)
            alpha_start = exp(phi0_values[j - 1]);

        alpha_flag = 0;
        //vertical searching of manually entered trouble range using variety of frequency guesses
       if (A_central > trouble_end && guess_layer <= max_guess_layer)
       {
            A_central = trouble_start;
            alpha_central = alpha_start;
            alpha_end = exp(phi0_values[j - 1]);
            frequency_guess -= f_search_interval;
            cout << frequency_guess << endl;
            dA_factor = 0.1;
            guess_layer++;
            alpha_flag = 1;
       }

       if (A_central > trouble_end && guess_layer > max_guess_layer)
            {frequency_guess = f_start; dA_factor = 1.; alpha_central = alpha_end; alpha_flag = 1;}




        //updates resolution according to # of threshold A-values passed
        while ( threshold_counter < refine_thresholds.size() - 1 && A_central > refine_thresholds[threshold_counter])
        {
            //n_gridpoints *= 2;

            threshold_counter++;

            cout << "n_gridpoints = " << n_gridpoints << " at threshold " << threshold_counter << endl;
        }

        if (threshold_counter > 0)
        {

            double A_upper = refine_thresholds[threshold_counter];
            double A_lower =  refine_thresholds[threshold_counter - 1];
            double k = log(2) / (A_upper - A_lower); //factor in exponential to ensure number of gridpoints doubles exponentially/smoothly between thresholds


            n_gridpoints = ceil(N_gridpoints_init * pow(2., threshold_counter - 1) * exp (k * (A_central - A_lower)));

            //n_gridpoints_no_error *= exp(k * delta_A);
            //n_gridpoints = ceil(n_gridpoints_no_error);

        }

        //cout << state[0].phi << endl;

         if (!solve(1) )
        {
            //cout << "Failed to solve with A = " << A_central << endl;
            omega = 0;
            //exit(1);
        }
        else data_file << std::setprecision (10) <<  A_central << "     " << M << "     " << r_99 << "     " << noether_charge <<  "     " << binding_energy
        << "     "  << static_cast<double>(omega) << "     " << static_cast<double>(omega_pre_rescale) << "     " << state[0].phi << endl;

        A_values[j] = A_central;
        omega_values[j]  = omega;
        phi0_values[j] = state[0].phi;

        //cout << frequency_guess << endl;


    }


}


#endif /* BOSONSTAR_CPP_ */
