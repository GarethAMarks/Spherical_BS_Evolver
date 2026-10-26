#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#ifndef DIRACSTAR_CPP_
#define DIRACSTAR_CPP_

#include "DiracStar.h"
#include "mathutils.h"
#include "mpi_helpers.h"
#include <sstream>
#include<iomanip>

using namespace std;

//overloads for addition/subtraction/scalar multiplication of DiracFieldState vars, useful for simplifying RK4 code
DiracFieldState operator+(const DiracFieldState& s1, const DiracFieldState& s2)
{
    return DiracFieldState{s1.F + s2.F, s1.G + s2.G, s1.X + s2.X, s1.phi + s2.phi};
}

DiracFieldState operator-(const DiracFieldState& s1, const DiracFieldState& s2)
{
    return DiracFieldState{s1.F - s2.F, s1.G - s2.G, s1.X - s2.X, s1.phi - s2.phi};
}

DiracFieldState operator*(double c, const DiracFieldState& s)
{
    return DiracFieldState{c * s.F, c * s.G, c * s.X, c * s.phi};
}

//read parameters in from file (default BSParams.par); uses the same keys as BosonStar/ProcaStar
void DiracStar::read_parameters(bool quiet, string filename)
{
    ifstream params{ filename };

    // Print an error and exit if file cannot open
    if (!params)
    {
        std::cerr << "Could not open " << filename << "\n";
        abort();
    }

    string current_line{};

    while (getline(params, current_line))
    {
        fill_parameter(current_line, "D = ", D, quiet);
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
        fill_parameter(current_line, "n_gridpoints = ", n_gridpoints, quiet);
        fill_parameter(current_line, "courant_factor = ", courant_factor, quiet);
        fill_parameter(current_line, "stop_time = ", stop_time, quiet);

        fill_parameter(current_line, "cycle_only = ", cycle_only, quiet);
        fill_parameter(current_line, "A0 = ", A0, quiet);
        fill_parameter(current_line, "dA = ", dA, quiet);
        fill_parameter(current_line, "n_stars = ", n_stars, quiet);
    }

    check_supported();
}

//the Dirac star solver is implemented for the free field in D = 4 only
void DiracStar::check_supported()
{
    if (D != 4.)
    {
        cerr << "ERROR: Dirac stars are only implemented for D = 4 (the number of spinors and their angular structure depend on D)" << endl;
        exit(1);
    }
    if (lambda != 0. || solitonic)
    {
        cerr << "ERROR: Dirac stars are only implemented for the free field: set lambda = 0 and solitonic = 0" << endl;
        exit(1);
    }
}

//number density n, energy density rho and radial pressure p_r of the two spinors at radius r (orthonormal frame);
//at r = 0 the regular limits G/r -> g0 sqrt(alpha0), F/r -> 0, F/r^2 -> (mu - omega/alpha0) g0 sqrt(alpha0) / 3 are used
void DiracStar::matter_terms(double r, long double frequency, const DiracFieldState& s, double& rho, double& p_r, double& n)
{
    double alpha = exp(s.phi);
    double w_a = static_cast<double>(frequency) / alpha;

    double Gr, Fr, FG_r3; //G/r, F/r and F G / r^3
    if (r == 0.)
    {
        Gr = A_central * sqrt(alpha);
        Fr = 0.;
        FG_r3 = (mu - w_a) * Gr * Gr / 3.;
    }
    else
    {
        Gr = s.G / r;
        Fr = s.F / r;
        FG_r3 = Fr * Gr / r;
    }

    double norm = 1. / (2. * M_PI * alpha);
    n = norm * (Fr * Fr + Gr * Gr);
    rho = w_a * n;
    p_r = norm * (w_a * (Fr * Fr + Gr * Gr) - mu * (Gr * Gr - Fr * Fr) + 2. * FG_r3);
}

//right-hand side of the Einstein-Dirac system
DiracFieldState DiracStar::state_RHS(const double radius, const long double frequency, DiracFieldState s, bool asymptotic_region)
{
    const double& r = radius;
    double rho, p_r, n;
    matter_terms(r, frequency, s, rho, p_r, n);

    double w_a = static_cast<double>(frequency) * exp(-s.phi);
    double T1 = (r == 0.) ? 0. : 0.5 * (s.X * s.X - 1.) / r;

    double dX = s.X * (4. * M_PI * G * r * s.X * s.X * rho - T1);
    double dPhi = T1 + 4. * M_PI * G * r * s.X * s.X * p_r;

    //in the asymptotic region, evolve phi and X normally but do not update F, G (these will be hard-coded to asymptotic expressions)
    if (asymptotic_region)
        return DiracFieldState{0., 0., dX, dPhi};

    double Gr = (r == 0.) ? A_central * exp(0.5 * s.phi) : s.G / r;
    double Fr = (r == 0.) ? 0. : s.F / r;

    double dF = s.X * (-Fr + (mu - w_a) * s.G);
    double dG = s.X * (Gr + (mu + w_a) * s.F);

    return DiracFieldState{dF, dG, dX, dPhi};
}

void DiracStar::rk4_solve (const long double freq)
{
    state = {DiracFieldState{0., 0., 1.0, log(alpha_central)}};
    radius_array = {0.};

    blowup_point = n_gridpoints; //make this 1 larger than max possible value to start, in case solution does not break

    state.resize(n_gridpoints);
    radius_array.resize(n_gridpoints);

    double dr = R / (n_gridpoints - 1);

    //inter-level state values for RK4 evolution
    DiracFieldState s1, s2, s3, s4;

    //fill in grid using RK4 scheme
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        s1 = state_RHS(r, freq, state[j], 0);
        s2 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s1, 0);
        s3 = state_RHS(r + dr / 2., freq, state[j] + 0.5 * dr * s2, 0);
        s4 = state_RHS(r + dr, freq, state[j] + dr * s3, 0);

        //update state variables and radius array
        state[j + 1] = state[j] + (dr / 6.) * (s1 + 2 * s2 + 2 * s3 + s4);
        radius_array[j + 1] = (j + 1) * dr;

        if (isnan(state[j].F) || isnan(state[j].G) || isnan(state[j].X) || isnan(state[j].phi)
            || isinf(state[j].G) || isinf(state[j].X))
        {
            blowup_point = j;
            return;
        }
    }
}

//aspect mass m(r) = r (1 - 1/X^2) / 2 (D = 4)
double DiracStar::m(int j)
{
    if (j >= n_gridpoints || j < 0)
        {
            cerr << "ERROR: invalid index passed to m(r)" << endl;
            exit(1);
        }

    return 0.5 * radius_array[j] * (1. - (1. / (state[j].X * state[j].X)));
}

//writes field values to DSdata.dat; f and g are the physical radial spinor functions F/(sqrt(alpha) r), G/(sqrt(alpha) r)
void DiracStar::write_field(string filename)
{
    std::ofstream data_file{filename};

    // If we couldn't open the output file stream for writing
    if (!data_file)
    {
        std::cerr << filename << " could not be opened for writing!\n";
        exit(1);
    }

    data_file << "# r F G X phi m f g\n";

    for (int j = 0; j < n_gridpoints - 1; ++j)
    {
        double r = radius_array[j];
        double sa = exp(0.5 * state[j].phi);
        double f = (j == 0) ? 0. : state[j].F / (sa * r);
        double g = (j == 0) ? A_central : state[j].G / (sa * r);

        data_file << std::setprecision (10) << r << "   " << state[j].F << "    " << state[j].G << "    " << state[j].X << "    "
                  << state[j].phi << "    " << m(j) << "    " << f << "    " << g << std::endl;
    }
}

//convergence test for the RK4 solver; writes differences in the physical upper spinor function g = G/(sqrt(alpha) r)
//between resolutions to conv_dirac.dat
void DiracStar::convergence_test(long double freq)
{
    std::vector<DiracFieldState> state_l, state_m, state_h;
    std::vector<double> r_l;

    if (freq == 0.) freq = omega_pre_rescale;

    solve(1);
    cout << setprecision(10) << static_cast<double>(omega) << endl;
    state_l = state; r_l = radius_array;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1);
    state_m = state;
    cout <<  setprecision(10) << static_cast<double>(omega) << endl;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1);
    state_h = state;
    cout << setprecision(10) << static_cast<double>(omega) << endl;

    std::ofstream conv_file{ "conv_dirac.dat" };

    if (!conv_file)
    {
        std::cerr << "Error: conv_dirac.dat could not be opened for writing\n";
        exit(1);
    }

    auto gp = [&](const DiracFieldState& s, double r) {return (r == 0.) ? A_central : s.G * exp(-0.5 * s.phi) / r;};

    for (int j = 0; j < (n_gridpoints + 3) / 4; j++)
    {
        double r = radius_array[4 * j];
        conv_file << r << "   " << gp(state_m[2*j], r) - gp(state_l[j], r) << "    " << 16.*(gp(state_h[4*j], r) - gp(state_m[2*j], r)) << endl;
    }
}

//returns number of zero crossings in G. Must call rk4_solve() first!
int DiracStar::count_zero_crossings()
{
    int zero_crossings = 0;

    //check for zero crossings up to blowup point (neglecting one extra point as we sometimes get spurious zero crossings on explosion)
    for (int j = 2; j < blowup_point - 1; j++)
    {
        if ((state[j].G == 0.) || (state[j].G > 0. && state[j - 1].G < 0.) || (state[j].G < 0. && state[j - 1].G > 0.) )
            {zero_crossings++;}
    }

    return zero_crossings;
}

//interval bisection algorithm-- initial_guess must be greater than true value for search to work
long double DiracStar::find_frequency(bool quiet)
{
    //lower bound on frequency
    long double lower_guess = frequency_guess;

    //tolerance for uncertainty in frequency, as well as minimum frequency to try
    long double epsilon = freq_epsilon;

    //at first, compute models, halving frequency each time until we find one with fewer/equal zero crossings to the number desired
    rk4_solve(lower_guess);
    if (count_zero_crossings() <= eigen && !quiet)
    {
        cerr << "WARNING: initial guess of "<< frequency_guess << " too small for A_central = " << A_central << endl;
    }

    while (count_zero_crossings() > eigen && lower_guess > epsilon)
    {
        lower_guess /= 2;
        rk4_solve(lower_guess);
    }

    if (lower_guess <= epsilon)
    {
        cerr << "ERROR: could not find suitable lower frequency bound for A_central = " << A_central << endl;
    }

    //upper bound on frequency
    long double upper_guess = lower_guess * 2.;
    long double midpoint = 0.5 * (upper_guess + lower_guess);

    int num_solves = 0;

    //now use interval bisection to converge to correct frequency
    while ( (upper_guess - lower_guess) > epsilon && num_solves < 100)
    {
        midpoint = 0.5 * (upper_guess + lower_guess);
        rk4_solve(midpoint);
        if (count_zero_crossings() > eigen)
            upper_guess = midpoint;
        else
            lower_guess = midpoint;

        num_solves++;
    }

    omega = lower_guess; //use lower_guess for frequency so we are guaranteed to have right # zero crossings in principle...
    rk4_solve(omega);

    if (count_zero_crossings() != eigen)
        cerr << "WARNING: zero crossing count may be incorrect for A_central = " << A_central << ", found " << count_zero_crossings() << endl;

    if (!quiet) cout << " \nFound solution in eigenstate " << eigen  << " with frequency " << static_cast<double>(omega) << endl;

    return lower_guess;
}

//returns the index of the last minimum of |G|, searching inwards from the blowup point.
int DiracStar::find_last_minimum()
{
    //G vanishes at the origin, so compare against its maximum magnitude rather than its central value
    double G_max = 0.;
    for (int j = 0; j < blowup_point - 1 && j < n_gridpoints; j++)
        G_max = max(G_max, abs(state[j].G));

    for (int j = blowup_point - 2; j > 1; j--)
    {
     if (abs(state[j].G) <= abs(state[j - 1].G) && abs(state[j].G) <= abs(state[j + 1].G) && abs(state[j].G) < 0.05 * G_max  )
        return j;
    }

    //case where our frequency guess is good enough that the solution does not blow up within the grid
    if (blowup_point == n_gridpoints && abs(state[n_gridpoints - 1].G) <= abs(state[n_gridpoints - 2].G ))
        return (n_gridpoints - 1);

    return -1; //signifies that next function should skip
}

//fill up the solution in the region after the blowup point, continuing F and G as decaying exponentials and integrating
//phi and X. Returns 1 if successful
bool DiracStar:: fill_asymptotic(bool quiet)
{
    int min_index = find_last_minimum();

    if (min_index == -1 && !quiet)
        cerr << "ERROR: No local minimum in |G| found for A_central = " << A_central << endl;

    if (min_index == -1)
        return 0;

    int j_match = round( r_match_fac * min_index); //index at which matching takes place

    double dr = R / (n_gridpoints - 1);
    double r_match = dr * j_match; //matching radius

    //far from the star F and G decay as exp(-k X r) with k^2 = mu^2 - omega^2 / alpha^2 (to leading order; the ratio F/G is
    //taken from the matching point so that both are continuous)
    double alpha_match = exp(state[j_match].phi);
    double k = state[j_match].X * sqrt( mu * mu - static_cast<double>(pow(omega / alpha_match, 2)));
    double F_match = state[j_match].F;
    double G_match = state[j_match].G;

    if (isnan(k))
    {
        if (!quiet) cerr << "ERROR: omega / alpha exceeds mu at the matching point for A_central = " << A_central << endl;
        return 0;
    }

    DiracFieldState s1, s2, s3, s4;
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

        //fix asymptotic values for F, G
        double decay = exp(-k * (radius_array[j + 1] - r_match));
        state[j + 1].F = F_match * decay;
        state[j + 1].G = G_match * decay;

        if (isnan(state[j].F) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].G))
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
    double cutoff_dev = 0.01;
    while (j < n_gridpoints - 1 && m(j) < (1 - cutoff_dev) * M)
    {
        //just uses linear interpolation to find r_99 for now
        if ( m(j + 1) != m(j) )
            r_99 = dr * j + dr * ((1 - cutoff_dev) * M - m(j) ) / (m(j + 1) - m(j));

        j++;
    }

     noether_charge = get_noether_charge();
     binding_energy = M - mu * noether_charge;
     compactness = M / r_99;

     if (!quiet)
     {
        cout << "\nFinal frequency after enforcing lapse condition is " << static_cast<double>(omega) << endl;
        cout << "\nMass is M =  " << M << endl;
        cout << "\nNoether charge (particle number of both spinors) is N = " << noether_charge << endl;
        cout << "\nBinding energy is E = " << binding_energy << endl;
     }

     return 1;
}

//rescales lapse by k = e^(phi_shift) as well as the frequency. F and G carry a factor alpha^{1/2}, so they scale by sqrt(k)
//(which leaves the physical spinor, n, rho and p_r unchanged)
void DiracStar::rescale_lapse (double phi_shift)
{
    double fac = exp(0.5 * phi_shift);

    for (DiracFieldState& s: state)
    {
        s.phi += phi_shift;
        s.F *= fac;
        s.G *= fac;
    }

    omega *= exp(phi_shift);
}

//simply calls find_frequency and fill_asymptotic in one; return 1 iff successful
bool DiracStar::solve(bool quiet)
{
    find_frequency(quiet);
    return fill_asymptotic(quiet);
}

//returns the total particle number of both spinors, Q = int n sqrt(gamma) d^3x = int 2 (F^2 + G^2) X / alpha dr.
//Must have computed model + frequency first (polar)
double DiracStar::get_noether_charge()
{
    double Q = 0.; //charge to return
    double dr = R / (n_gridpoints - 1);

    //trapezoid rule (the integrand vanishes at r = 0)
    for (int j = 1; j < n_gridpoints; j++)
    {
        double w = (j == n_gridpoints - 1) ? 0.5 : 1.;
        Q += w * dr * 2. * (state[j].F * state[j].F + state[j].G * state[j].G) * state[j].X * exp(-state[j].phi);
    }

    return Q;
}

//get compactness corresponding to Collodel et. al's C4 definition, defined using the radius within which X is within 0.999 of its Schwarzchild value. Must have computed M!
double DiracStar::get_c4()
{
    double r4 = R;
    double cutoff_dev = 0.001; //relative deviation from Schwarzschild at which we define the boundary

    for (int j = n_gridpoints - 3; j > 0; j--)
    {
        double dr = R / (n_gridpoints - 1);
        double r = j * dr;

        double BH_X_sq = 1. /(1. - 2. * M / r);
        double DS_X_sq = state[j].X * state[j].X;

        if (abs((DS_X_sq - BH_X_sq) / BH_X_sq) > cutoff_dev )
        {
            //linearly interpolate to estimate exactly where the relative deviation = cutoff_dev
            double f1 = abs((DS_X_sq - BH_X_sq) / BH_X_sq);

            double BH_X2_sq = 1. /(1. - 2. * M / (r + dr) );
            double f2 = abs((state[j + 1].X * state[j + 1].X  - BH_X2_sq) / BH_X2_sq);

            r4 = r + dr * (cutoff_dev - f1) / (f2 - f1);
            break;
        }

        if (j == 1)
            cout <<"WARNING: get_c4() failed to establish deviation from Schwarzschild BH!" << endl;
    }

    return (M / r4);
}

//return the maximum value of the compactness as defined by c = m(r) / r with m the aspect mass.
double DiracStar::get_cmax()
{
    double cmax = 0.;

    for (int j = 1; j < n_gridpoints - 1; j++)
    {
        double dr = R / (n_gridpoints - 1);
        double r = j * dr;
        double c = m(j) / r;

        if (c > cmax )
            cmax = c;
    }

    return cmax;
}

//returns right-hand side for the ODE that f (ratio of isotropic to areal radius) solves. r is areal radius
double DiracStar::f_RHS(double r, double f)
{
    double dr = R / (n_gridpoints - 1);
    int j_low = floor((r / R) * (n_gridpoints - 1));
    int j0 = bound(j_low - 1, 0, n_gridpoints - 4);

    double X = cubic_interp(r, state[j0].X, state[j0 + 1].X, state[j0 + 2].X, state[j0 + 3].X, j0, dr);

    double r2 = (r == 0.) ? 10e-10 : r;
    return f * (X - 1.) / r2;
}

//returns areal radius corresponding to a given isotropic index
double DiracStar::r_areal(int j_iso)
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
        return r_iso + M + M * M / (4. * r_iso);

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
}

//fill out arrays for radius, spinor functions and conformal factor in isotropic coordinates. Must have solved model first!
//The physical radial spinor functions f = F/(sqrt(alpha) r), g = G/(sqrt(alpha) r) are orthonormal-frame components and so
//are unchanged by the radial relabelling r -> R (the radial tetrad leg is the same unit vector).
void DiracStar::fill_isotropic_arrays()
{
    double dr = R / (n_gridpoints - 1);

    double f1, f2, f3, f4; //intermediate rk4 values

    //physical spinor functions on the areal grid
    vector<double> f_phys(n_gridpoints), g_phys(n_gridpoints);
    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = j * dr;
        double sa = exp(0.5 * state[j].phi);
        f_phys[j] = (j == 0) ? 0. : state[j].F / (sa * r);
        g_phys[j] = (j == 0) ? A_central : state[j].G / (sa * r);
    }

    //local array to store values of f, the ratio between isotropic and areal radii; assumed = 1 at the origin for now and later to be rescaled
    vector<double> f_array(n_gridpoints);
    vector<double> r_areal_array(n_gridpoints);
    f_array[0] = 1.;

    r_iso_array.resize(n_gridpoints);
    psi_iso_array.resize(n_gridpoints);
    phi_iso_array.resize(n_gridpoints);
    f_iso_array.resize(n_gridpoints);
    g_iso_array.resize(n_gridpoints);

    r_iso_array[0] = 0;

    //evolve f using RK4 evolution and fill out R_iso array
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        f1 = f_RHS(r, f_array[j]);
        f2 = f_RHS(r + dr / 2., f_array[j] + 0.5 * dr * f1);
        f3 = f_RHS(r + dr / 2.,  f_array[j] + 0.5 * dr * f2);
        f4 = f_RHS(r + dr, f_array[j] + dr * f3);

        f_array[j + 1] = f_array[j] + (dr / 6.) * (f1 + 2 * f2 + 2 * f3 + f4);
        r_iso_array[j + 1] = f_array[j + 1] * (r + dr);
    }

    //rescale factors to match R_iso to asymptotic Schwarzchild at boundary
    double iso_scale_factor = (R - M - M *  M / (4. * ( R -  M))) / r_iso_array[n_gridpoints - 1];

    for (double& r_iso: r_iso_array)
        r_iso *= iso_scale_factor;

    for (double& f: f_array)
        f *= iso_scale_factor;

    //fill in spinor, phi, psi arrays
    for (int j = 0; j < n_gridpoints; j++)
    {
        r_areal_array[j] = r_areal(j);

        int j_areal_low = floor(r_areal_array[j] / dr);
        int j_areal_high = ceil (r_areal_array[j] / dr);

        double f;

        //cubic interpolation where we won't go OOB in radial array
        if (j_areal_high <= n_gridpoints - 1)
        {
            int j0 = bound(j_areal_low - 1, 0, n_gridpoints - 4);

            phi_iso_array[j] = cubic_interp(r_areal_array[j], state[j0].phi, state[j0 + 1].phi, state[j0 + 2].phi, state[j0 + 3].phi, j0, dr );
            f_iso_array[j] = cubic_interp(r_areal_array[j], f_phys[j0], f_phys[j0 + 1], f_phys[j0 + 2], f_phys[j0 + 3], j0, dr );
            g_iso_array[j] = cubic_interp(r_areal_array[j], g_phys[j0], g_phys[j0 + 1], g_phys[j0 + 2], g_phys[j0 + 3], j0, dr );
            f = cubic_interp(r_areal_array[j], f_array[j0], f_array[j0 + 1], f_array[j0 + 2], f_array[j0 + 3], j0, dr );
        }

        else //extrapolation case; will kick in near boundary
        {
            int N = n_gridpoints;
            double dR = r_areal_array[j] - R;
            phi_iso_array[j] = state[N - 1].phi  + dR * (state[N - 1].phi - state[N - 2].phi ) / dr;
            f_iso_array[j] = f_phys[N - 1]  + dR * (f_phys[N - 1] - f_phys[N - 2] ) / dr;
            g_iso_array[j] = g_phys[N - 1]  + dR * (g_phys[N - 1] - g_phys[N - 2] ) / dr;
            f = f_array[N - 1]  + dR * (f_array[N - 1] - f_array[N - 2] ) / dr;
        }

        //linearly transition between interpolated and asymptotic f in this region, as in BosonStar
        double l_frac = 0.8;
        double u_frac = 0.9;

        if (j > n_gridpoints * l_frac)
        {
            double portion_crossed = min((j - n_gridpoints * l_frac) / (n_gridpoints * (u_frac - l_frac)), 1.);

            f = f * (1 - portion_crossed) + (pow((1 + M / (2. * dr * j)), -2.)) * portion_crossed;
        }

        //psi is simply 1 / sqrt(f)
        psi_iso_array[j] = 1. / sqrt(f);
    }
}

//writes isotropic values to isotropic_dirac.dat
void DiracStar::write_isotropic()
{
    std::ofstream data_file{ "isotropic_dirac.dat" };

    if (!data_file)
    {
        std::cerr << "isotropic_dirac.dat could not be opened for writing!\n";
        exit(1);
    }

    data_file << "# r r_iso r_areal psi_iso f_iso g_iso phi_iso\n";

    for (int j = 1; j < n_gridpoints - 1; j++)
    {
        data_file << std::setprecision (10) << radius_array[j] << "   " << r_iso_array[j] << "    " << r_areal(j) << "    " << psi_iso_array[j]  <<  "    "
        << f_iso_array[j] << "    " << g_iso_array[j] << "    " << phi_iso_array[j] << std::endl;
    }
}

//constructs a family of Dirac star models with A_central = A_0 + delta_A, A_0 + 2 delta_A, ..., writing to DiracStars.dat.
//Same columns as BosonStars.dat; noether_charge is the particle number of both spinors.
void DiracStar::cycle_models(int n_stars, double A_0, double delta_A)
{
    read_parameters(0);
    ofstream data_file{"DiracStars.dat"};
    data_file << "# A_central M r_99 noether_charge binding_energy omega omega_unrescaled phi0 C4 Cmax\n";
    omega_pre_rescale = frequency_guess;

    //phi(0) from previous models, used to extrapolate a guess for the central lapse
    vector<double> phi0_values;

    A_central = A_0;

    for (int j = 0; j < n_stars; j++)
    {
        A_central += delta_A;

        //initial central lapse guess: linear extrapolation in phi(0) from previous successful models
        if (phi0_values.size() == 1)
            alpha_central = exp(phi0_values.back());
        else if (phi0_values.size() >= 2)
            alpha_central = exp(2. * phi0_values[phi0_values.size() - 1] - phi0_values[phi0_values.size() - 2]);

        if (!solve(1))
        {
            cout << "Failed to solve with A = " << A_central << endl;
            omega = 0;
            continue;
        }

        data_file << std::setprecision (10) <<  A_central << "     " << M << "     " << r_99 << "     " << noether_charge <<  "     " << binding_energy
        << "     "  << static_cast<double>(omega) << "     " << static_cast<double>(omega_pre_rescale) << "     " << state[0].phi << "     " << get_c4()
        << "     " << get_cmax() << endl;

        phi0_values.push_back(state[0].phi);
    }
}

void DiracStar::broadcast()
{
#ifdef USE_MPI
    // Broadcast scalars that solve() sets.
    double buf[7] = {R, (double)omega, (double)omega_pre_rescale, M, r_99, noether_charge, A_central};
    bss_mpi::bcast_doubles(buf, 7);
    R               = buf[0];
    omega           = (long double)buf[1];
    omega_pre_rescale = (long double)buf[2];
    M               = buf[3];
    r_99            = buf[4];
    noether_charge  = buf[5];
    A_central       = buf[6];

    int ibuf[2] = {n_gridpoints, blowup_point};
    bss_mpi::bcast_ints(ibuf, 2);
    n_gridpoints = ibuf[0];
    blowup_point = ibuf[1];

    auto bcast_vec = [](std::vector<double>& v)
    {
        int sz = (int)v.size();
        bss_mpi::bcast_ints(&sz, 1);
        v.resize(sz);
        if (sz > 0)
            bss_mpi::bcast_doubles(v.data(), sz);
    };

    // state: vector<DiracFieldState>, each is 4 contiguous doubles (F, G, X, phi).
    {
        int sz = (int)state.size();
        bss_mpi::bcast_ints(&sz, 1);
        state.resize(sz);
        if (sz > 0)
            bss_mpi::bcast_doubles(reinterpret_cast<double*>(state.data()), sz * 4);
    }

    bcast_vec(psi_iso_array);
    bcast_vec(f_iso_array);
    bcast_vec(g_iso_array);
    bcast_vec(phi_iso_array);
    bcast_vec(r_iso_array);
    bcast_vec(radius_array);
#endif
}

#endif /* DIRACSTAR_CPP_ */
