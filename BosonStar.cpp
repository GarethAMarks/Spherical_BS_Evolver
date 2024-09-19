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
        return mu * mu * A * A;

    else
        return mu * mu * A * A * pow((1. - 2. * pow(A / sigma, 2)), 2);

}

double BosonStar::dV( const double A)
{
    if (!solitonic)
        return mu * mu;

    else
        return mu * mu - 8. * mu * mu * pow(A / sigma, 2) + 12. * mu * mu * pow(A / sigma, 4);

}


//helper function for read_parameters; searches current_line for line_start and fills in parameter
//TODO: add default parameters for not found case
/*template <typename T>
void fill_parameter (string& current_line, string line_start, T& parameter, bool quiet)
{
    if (current_line.find(line_start) != string::npos)
    {
        // Create a substring starting from the position right after line_start
        size_t pos = current_line.find(line_start);
        string rest_of_line = current_line.substr(pos + line_start.length());

        // Create a stringstream from the rest_of_line
        stringstream ss(rest_of_line);

        // Extract the double value from the stringstream
        if (ss >> parameter) {
            if (!quiet) cout << "Read in " << line_start << parameter << endl;
        } else if (!(ss >> parameter)) {
            cout << "WARNING: Failed to extract value for parameter " << line_start << endl;
        }
    }
}*/

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
        fill_parameter(current_line, "solitonic = ", solitonic, quiet);
        fill_parameter(current_line, "G = ", G, quiet);
        fill_parameter(current_line, "sigma = ", sigma, quiet);
        fill_parameter(current_line, "eigen = ", eigen, quiet);
        fill_parameter(current_line, "alpha_central = ", alpha_central, quiet);
        fill_parameter(current_line, "A_central = ", A_central, quiet);
        fill_parameter(current_line, "frequency_guess = ", frequency_guess, quiet);
        fill_parameter(current_line, "freq_epsilon = ", freq_epsilon, quiet);
        fill_parameter(current_line, "isotropic = ", isotropic, quiet);

        fill_parameter(current_line, "R = ", R, quiet);
        fill_parameter(current_line, "n_gridpoints = ", n_gridpoints, quiet);
        fill_parameter(current_line, "courant_factor = ", courant_factor, quiet);
        fill_parameter(current_line, "stop_time = ", stop_time, quiet);
    }
}

//right-hand side of the EKG system of equations
//IDEA: add term in alpha to eta (tentative... for now)
FieldState BosonStar::state_RHS(const double radius, const double frequency, FieldState  s, bool asymptotic_region)
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
    else
    {
        //return RHS of field state variables outside of asymptotic region
        return  (FieldState) {s.X * s.eta, s.X * ( F1 * (T2 + V(s.A)) - T1 ), dPhi,
        -2. * T3 - s.eta * dPhi + s.X * s.A * (dV(s.A) - frequency * frequency  / exp(2 * s.phi))};
    }

   //test for playing with RK4 convergence-- seems that 1/r terms can spoil 4th-order convergence (but we still get tight 3rd-order)
   //return  (FieldState) {s.A / r + s.X, s.X + s.phi, s.phi + s.eta, s.eta + s.A};

}

void BosonStar::rk4_solve (const double freq)
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
void BosonStar::convergence_test(double freq)
{
    //FieldState vectors for low, medium, high resolution, each off by factor of 2
    std::vector<FieldState> state_l, state_m, state_h;

    if (freq == 0.) freq = omega_pre_rescale; //if no frequency provided (for rk4 solver) use omega

    //fill in state vectors at low, mid, high resolution
    solve(1); // rk4_solve(freq);

    cout << setprecision(10) << omega << endl;

    state_l = state;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1); //rk4_solve(freq);
    state_m = state;

    cout <<  setprecision(10) << omega << endl;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1); //rk4_solve(freq);
    state_h = state;

    cout << setprecision(10) << omega << endl;

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
double BosonStar::find_frequency(bool quiet)
{
    //lower bound on frequency
    double lower_guess = frequency_guess;

    //tolerance for uncertainty in frequency, as well as minimum frequency to try
    double epsilon = freq_epsilon;

    //at first, compute BS models, halving frequency each time until we find one with fewer/equal zero crossings to the number desired
    rk4_solve(lower_guess);
    if (count_zero_crossings() <= eigen)
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
    double upper_guess = lower_guess * 2.;
    double midpoint = 0.5 * (upper_guess + lower_guess);

    //now use interval bisection to converge to correct frequency
    while ( (upper_guess - lower_guess) > epsilon)
    {
        //replace upper/lower bound with midpoint if it has greater/as many or fewer zero crossings than desired, so both ultimately converge on boundary value between eigen and eigen + 1 zero crossings as needed
        midpoint = 0.5 * (upper_guess + lower_guess);
        rk4_solve(midpoint);
        if (count_zero_crossings() > eigen)
            upper_guess = midpoint;
        else
            lower_guess = midpoint;
    }

    omega = lower_guess; //use lower_guess for frequency so we are guarenteed to have right # zero crossings in principle...
    rk4_solve(omega);

    if (count_zero_crossings() != eigen)
        cerr << "WARNING: zero crossing count may be incorrect for A_central = " << A_central << ", found " << count_zero_crossings() << endl;


    if (!quiet) cout << " \nFound solution in eigenstate " << eigen  << " with frequency " << omega << endl;


    return lower_guess;


}

//returns the index of the last minimum of |A|, searching inwards from the blowup point.
int BosonStar::find_last_minimum()
{

    for (int j = blowup_point - 2; j > 1; j--)
    {
     if (abs(state[j].A) <= abs(state[j - 1].A) && abs(state[j].A) <= abs(state[j + 1].A)  )
        return j;
    }

    //case where our frequency guess is good enough that the solution does not blow up within the grid
    if (blowup_point == n_gridpoints && abs(state[n_gridpoints - 1].A) <= abs(state[n_gridpoints - 2].A ))
        return (n_gridpoints - 1);

    cerr << "ERROR: No local minimum in |A| found for A_central = " << A_central << endl;
    return -1; //signifies that next function should skip
}
//fill up the solution in the region after the blowup point, using asymptotic matching for A and eta and integrating phi and X. Returns 1 if successful
bool BosonStar:: fill_asymptotic(bool quiet)
{

    double min_index = find_last_minimum();

    if (min_index == -1)
        return 0;

    double r_match_fac = 0.9;
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
        double r = j * dr;

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

     if (!quiet)
     {
        cout << "\nFinal frequency after enforcing lapse condition is " << omega << endl;
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
//Currently using a separate array to give r_areal(r_iso) via RK4, so this is deprecated!
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
        double l_frac = 0.7;
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
    //sigma = 0.1: {0.2, 0.275}
    vector<double> refine_thresholds{0.25,0.325,0.375, 0.425, 0.475, 0.525, 0.575};
    //bool passed_last_threshold; //set to 1 after last threshold reached

    //frequency from previous guess for use in update
    double omega_prev = frequency_guess;

    //holds all unrescaled omega / A / phi[0] values yet computed
    vector<double> omega_values(n_stars);
    vector<double> A_values(n_stars);
    vector<double> phi0_values(n_stars);

    double phi0_prev = 0.;

    bool var_dA_method = 0;

    double guess_buffer = 0.00005;


    //double n_gridpoints_no_error = n_gridpoints; //stores number of gridpoints with exponential correction keeping the fractional part to avoid roundoff error accumulating

    for (int j = 0; j < n_stars; j++)
    {

        if (!var_dA_method)
            A_central = A_0 + j * delta_A;


        //writes info about list of solved boson stars
        if (j!= 0)
            {
                //cout << omega_pre_rescale_prev << " " << omega_pre_rescale << endl;
                //frequency_guess = omega_pre_rescale + max(1., 2. * (omega_pre_rescale - omega_pre_rescale_prev));
                 //frequency_guess = omega + abs(omega_pre_rescale - omega) * 2;
                //cout << frequency_guess << endl;

                if(j > 2)
                {
                    frequency_guess = omega_values[j - 1] +  (omega_values[j - 1] - omega_values[j - 2]) + guess_buffer;
                    alpha_central = exp(2 * phi0_values[j - 1] - phi0_values[j - 2]);
                }

                //alpha_central = exp(2 * state[0].phi - phi0_prev);
                omega_prev = omega;
                phi0_prev = state[0].phi;

            }

        if (j > 2 && var_dA_method)
        {


            double theta = atan((omega_values[j - 1] - omega_values[j - 2]) / (A_values[j - 1] - A_values[j - 2]));


            //cout << theta << endl;

            //linear extrapolation works fine for angles not too steep
            if (/*abs(theta) < 31. * M_PI / 32.*/ j > 20 )
                {
                    A_central = A_values[j - 1] + delta_A * cos(theta);
                    frequency_guess = omega_values[j - 1] +  sin(theta) * delta_A + guess_buffer;
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

       // if (j == 10)
            //var_dA_method = 1;

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

         if (!solve(1) )
            cout << "Failed to solve with A = " << A_central << endl;

        A_values[j] = A_central;
        omega_values[j]  = omega;
        phi0_values[j] = state[0].phi;

        //cout << frequency_guess << endl;
        data_file << std::setprecision (10) <<  A_central << "     " << M << "     " << r_99 << "     " << noether_charge <<  "     " << binding_energy << "     "  << omega << "     " << omega_pre_rescale << "     " << frequency_guess << endl;

    }
}


#endif /* BOSONSTAR_CPP_ */
