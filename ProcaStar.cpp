#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#ifndef PROCASTAR_CPP_
#define PROCASTAR_CPP_

#include "ProcaStar.h"
#include "mathutils.h"
#include "mpi_helpers.h"
#include <sstream>
#include<iomanip>

using namespace std;


//overloads for addition/subtraction/scalar multiplication of ProcaFieldState vars, useful for simplifying RK4 code
ProcaFieldState operator+(const ProcaFieldState& s1, const ProcaFieldState& s2)
{
    return ProcaFieldState{s1.f + s2.f, s1.X + s2.X, s1.phi + s2.phi, s1.E + s2.E};
}

ProcaFieldState operator-(const ProcaFieldState& s1, const ProcaFieldState& s2)
{
    return ProcaFieldState{s1.f - s2.f, s1.X - s2.X, s1.phi - s2.phi, s1.E - s2.E};
}

ProcaFieldState operator*(double c, const ProcaFieldState& s)
{
    return ProcaFieldState{c * s.f, c * s.X, c * s.phi, c * s.E};
}

//potential and its derivatives wrt Psi = A_a Abar^a; same family as BosonStar with |phi|^2 -> Psi
double ProcaStar::V( const double Psi)
{
    if (!solitonic)
        return mu * mu * Psi + lambda * Psi * Psi;

    else
        return mu * mu * Psi * pow((1. - 2. * Psi / (sigma * sigma)), 2) + lambda * Psi * Psi;
}

double ProcaStar::dV( const double Psi)
{
    if (!solitonic)
        return mu * mu + 2. * Psi * lambda;

    else
        return mu * mu - 8. * mu * mu * Psi / (sigma * sigma) + 12. * mu * mu * Psi * Psi / pow(sigma, 4) + 2. * Psi * lambda;
}

double ProcaStar::ddV( const double Psi)
{
    if (!solitonic)
        return 2. * lambda;

    else
        return - 8. * mu * mu / (sigma * sigma) + 24. * mu * mu * Psi / pow(sigma, 4) + 2. * lambda;
}

//read parameters in from file (default BSParams.par); uses the same keys as BosonStar
void ProcaStar::read_parameters(bool quiet, string filename)
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

        fill_parameter(current_line, "cycle_only = ", cycle_only, quiet);
        fill_parameter(current_line, "A0 = ", A0, quiet);
        fill_parameter(current_line, "dA = ", dA, quiet);
        fill_parameter(current_line, "n_stars = ", n_stars, quiet);
    }
}

//solves the algebraic (radial) field equation V'(Psi) g = -omega E / alpha^2 for g, with Psi = -f^2/alpha^2 + g^2/X^2.
//Also returns W = V'(Psi) and Psi by reference. Returns nan if no solution is found (only possible with self-interactions).
double ProcaStar::solve_g(const long double frequency, const ProcaFieldState& s, double& W, double& Psi)
{
    double alpha2 = exp(2. * s.phi);
    double X2 = s.X * s.X;
    double src = static_cast<double>(-frequency * s.E / alpha2); // W * g = src
    double Psi_t = -s.f * s.f / alpha2; //temporal contribution to Psi

    //free field: explicit
    if (!self_interacting())
    {
        W = mu * mu;
        double g = src / W;
        Psi = Psi_t + g * g / X2;
        return g;
    }

    if (src == 0.)
    {
        Psi = Psi_t; W = dV(Psi);
        return 0.;
    }

    //Newton iteration on h(g) = V'(Psi(g)) g - src, starting from the solution neglecting g^2 in Psi
    double W_t = dV(Psi_t);
    double g = (W_t != 0.) ? src / W_t : src / (mu * mu);

    for (int it = 0; it < max_newton_iterations; it++)
    {
        Psi = Psi_t + g * g / X2;
        W = dV(Psi);
        double h = W * g - src;
        double dh = W + 2. * g * g * ddV(Psi) / X2;

        if (dh == 0. || isnan(dh) || isnan(h))
            {g_solve_failed = 1; return NAN;}

        double step = h / dh;
        g -= step;

        if (abs(step) <= newton_tol * abs(g))
        {
            Psi = Psi_t + g * g / X2;
            W = dV(Psi);
            return g;
        }
    }

    g_solve_failed = 1;
    return NAN; //no convergence: treat as breakdown of the solution
}

//right-hand side of the Einstein-Proca system of equations
ProcaFieldState ProcaStar::state_RHS(const double radius, const long double frequency, ProcaFieldState s, bool asymptotic_region)
{
    //enforce minimum radius epsilon if needed
    double epsilon = 0.0000001;
    double r = ((radius == 0.) ? epsilon : radius);

    double W, Psi;
    double g = solve_g(frequency, s, W, Psi);
    double alpha2 = exp(2. * s.phi);

    //helper terms that appear multiple times / should be zeroed at r = 0
    double T1 = 0.5 * (D - 3.) * (s.X * s.X - 1) / r;
    double T3 = s.E / r;

    double F1 = 4. * M_PI * G * r * s.X * s.X / (D - 2.);

    //zero out terms that should be zeroed at origin as E = 0, X = 1 there
    if (radius == 0.)
        {T1 = 0; T3 = 0;}

    //twice the energy density and radial pressure
    double E_term = s.E * s.E / (alpha2 * s.X * s.X);
    double two_rho = E_term + 2. * W * s.f * s.f / alpha2 + V(Psi);
    double two_p_r = -E_term + 2. * W * g * g / (s.X * s.X) - V(Psi);

    double dX = s.X * (F1 * two_rho - T1);
    double dPhi = T1 + F1 * two_p_r;
    double E_corr = (radius == 0.) ? (1. / (D - 1.)) : 1.;

    //in the asymptotic region, evolve phi and X normally but do not update f, E (these will be hard-coded to asymptotic expressions)
    if (asymptotic_region)
        return ProcaFieldState{0., dX, dPhi, 0.};

    //t-component of field equation: (r^{D-2} E / (alpha X))' = W f X r^{D-2} / alpha; f' from definition of E
    return ProcaFieldState{static_cast<double>(s.E + frequency * g), dX, dPhi,
                           E_corr * (W * s.f * s.X * s.X - (D - 2.) * T3 + s.E * (dPhi + dX / s.X))};
}

//returns the small-r expansion in X and alpha (2nd order) and f, E (leading order). Used at the second gridpoint for D > 4, as in BosonStar.
ProcaFieldState ProcaStar::state_expansion(const double radius, long double frequency)
{
    const double& r = radius;
    double alpha0 = alpha_central;
    double f0 = A_central * alpha_central;
    double w2 = static_cast<double>(frequency * frequency);

    double Psi0 = -A_central * A_central; // g(0) = 0
    double W0 = dV(Psi0);

    double E1 = W0 * f0 / (D - 1.); // E = E1 r + O(r^3)
    double f2 = E1 * (1. - w2 / (alpha0 * alpha0 * W0)); // f = f0 + f2 r^2 / 2 + O(r^4)

    double two_rho0 = 2. * W0 * f0 * f0 / (alpha0 * alpha0) + V(Psi0);
    double two_p_r0 = -V(Psi0);

    double X2 = 8. * M_PI * G * two_rho0 / ((D - 1.) * (D - 2.)); // X = 1 + X2 r^2 / 2
    double alpha2 = alpha0 * (0.5 * (D - 3.) * X2 + 4. * M_PI * G * two_p_r0 / (D - 2.)); // alpha = alpha0 + alpha2 r^2 / 2

    return ProcaFieldState{f0 + 0.5 * f2 * r * r,
                           1. + 0.5 * X2 * r * r,
                           log(alpha0 + 0.5 * alpha2 * r * r),
                           E1 * r};
}

void ProcaStar::rk4_solve (const long double freq)
{
    state = {ProcaFieldState{A_central * alpha_central, 1.0, log(alpha_central), 0.0 }};
    radius_array = {0.};

    blowup_point = n_gridpoints; //make this 1 larger than max possible value to start, in case solution does not break
    g_solve_failed = 0;
    g_breakdown_r = -1.;

    state.resize(n_gridpoints);
    radius_array.resize(n_gridpoints);

    double dr = R / (n_gridpoints - 1);

    //inter-level state values for RK4 evolution
    ProcaFieldState s1, s2, s3, s4;

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

        //use series expansions for X, Phi near origin for D > 4, as in BosonStar
        if (D > 4. && j < 1.)
        {
            ProcaFieldState small_state = state_expansion(r + dr, freq);

            // replace near-origin values with small-radius series expansion
            state[j + 1].X = small_state.X;
            state[j + 1].phi = small_state.phi;
        }

        if (isnan(state[j].f) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].E) || g_solve_failed)
        {
            blowup_point = j;
            if (g_solve_failed) g_breakdown_r = r;
            return;
        }
    }
}

double ProcaStar::m(int j)
{
    if (j >= n_gridpoints || j < 0)
        {
            cerr << "ERROR: invalid index passed to m(r)" << endl;
            exit(1);
        }

    return (D - 2.) * pow(M_PI, 0.5 * (D - 3.)) * pow(radius_array[j], D - 3.)   * (1. - (1. / (state[j].X * state[j].X))) / (8. * tgamma(0.5 * (D - 1.)));
}

//reconstructs g on the grid from the solved state
void ProcaStar::fill_g_array()
{
    g_array.resize(state.size());
    double W, Psi;

    for (size_t j = 0; j < state.size(); j++)
        g_array[j] = solve_g(omega, state[j], W, Psi);
}

//writes field values to PSdata.dat
void ProcaStar::write_field(string filename)
{
    std::ofstream data_file{filename};

    // If we couldn't open the output file stream for writing
    if (!data_file)
    {
        std::cerr << filename << " could not be opened for writing!\n";
        exit(1);
    }

    if (g_array.size() != state.size())
        fill_g_array();

    // Header
    data_file << "# r f g X phi E m\n";

    for (int j = 0; j < n_gridpoints - 1; ++j)
    {
        data_file << std::setprecision (10) << radius_array[j] << "   " << state[j].f << "    " << g_array[j] << "    " << state[j].X << "    "
                  << state[j].phi << "    " << state[j].E << "    " << m(j) << std::endl;
    }
}

void ProcaStar::double_resolution()
{
    n_gridpoints = 2 * n_gridpoints - 1;
}

//convergence test for the RK4 solver; writes differences in f/alpha (gauge-invariant) between resolutions to conv_proca.dat
void ProcaStar::convergence_test(long double freq)
{
    //ProcaFieldState vectors for low, medium, high resolution, each off by factor of 2
    std::vector<ProcaFieldState> state_l, state_m, state_h;

    if (freq == 0.) freq = omega_pre_rescale;

    solve(1);
    cout << setprecision(10) << static_cast<double>(omega) << endl;
    state_l = state;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1);
    state_m = state;
    cout <<  setprecision(10) << static_cast<double>(omega) << endl;

    n_gridpoints = 2 * n_gridpoints - 1;
    solve(1);
    state_h = state;
    cout << setprecision(10) << static_cast<double>(omega) << endl;

    std::ofstream conv_file{ "conv_proca.dat" };

    if (!conv_file)
    {
        std::cerr << "Error: conv_proca.dat could not be opened for writing\n";
        exit(1);
    }

    auto fa = [](const ProcaFieldState& s) {return s.f * exp(-s.phi);};

    //write change in f/alpha between med/low and high/med resolution to file (rescaled by 2^4 for 4th order convergence)
    for (int j = 0; j < (n_gridpoints +3) / 4; j++)
    {
        conv_file << radius_array[4*j] << "   " << fa(state_m[2*j]) - fa(state_l[j]) << "    " << 16.*(fa(state_h[4*j]) - fa(state_m[2*j])) << endl;
    }
}

//returns number of zero crossings in E (equivalently g) for PS solution. Must call rk4_solve() first!
//E plays the role that A plays for the BS: it is nodeless in the ground state (where f itself has one node, since
//g ~ -f(0) r near the origin but g and f share a sign asymptotically), stays positive below the eigenfrequency and
//crosses zero above it.
int ProcaStar::count_zero_crossings()
{
    int zero_crossings = 0;

    //check for zero crossings up to blowup point (neglecting one extra point as we sometimes get spurious zero crossings on explosion)
    for (int j = 1; j < blowup_point - 1; j++)
    {
        if ((state[j].E == 0.) || (state[j].E > 0. && state[j - 1].E < 0.) || (state[j].E < 0. && state[j - 1].E > 0.) )
            {zero_crossings++;}
    }

    return zero_crossings;
}

//interval bisection algorithm-- initial_guess must be greater than true value for search to work
long double ProcaStar::find_frequency(bool quiet)
{
    //lower bound on frequency
    long double lower_guess = frequency_guess;

    //tolerance for uncertainty in frequency, as well as minimum frequency to try
    long double epsilon = freq_epsilon;

    //at first, compute PS models, halving frequency each time until we find one with fewer/equal zero crossings to the number desired
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

//returns the index of the last minimum of |E|, searching inwards from the blowup point.
int ProcaStar::find_last_minimum()
{
    //E vanishes at the origin, so compare against its maximum magnitude rather than its central value
    double E_max = 0.;
    for (int j = 0; j < blowup_point - 1 && j < n_gridpoints; j++)
        E_max = max(E_max, abs(state[j].E));

    for (int j = blowup_point - 2; j > 1; j--)
    {
     if (abs(state[j].E) <= abs(state[j - 1].E) && abs(state[j].E) <= abs(state[j + 1].E) && abs(state[j].E) < 0.05 * E_max  )
        return j;
    }

    //case where our frequency guess is good enough that the solution does not blow up within the grid
    if (blowup_point == n_gridpoints && abs(state[n_gridpoints - 1].E) <= abs(state[n_gridpoints - 2].E ))
        return (n_gridpoints - 1);

    return -1; //signifies that next function should skip
}

//fill up the solution in the region after the blowup point, using asymptotic matching for f and E and integrating phi and X. Returns 1 if successful
bool ProcaStar:: fill_asymptotic(bool quiet)
{
    int min_index = find_last_minimum();

    if (min_index == -1 && !quiet)
        cerr << "ERROR: No local minimum in |E| found for A_central = " << A_central << endl;

    if (min_index == -1)
        return 0;

    int j_match = round( r_match_fac * min_index); //index at which matching takes place

    double dr = R / (n_gridpoints - 1);
    double r_match = dr * j_match; //matching radius

    double phi_match = state[j_match].phi; //match values in current gauge
    double E_match = state[j_match].E;

    //decay rate; far from the star f obeys f'' + (D-2) f' / r = k^2 f, as for the BS amplitude, and f' = E k^2 / mu^2
    double k = sqrt( mu * mu - static_cast<double>(pow(omega / exp(phi_match), 2)));

    //ratio E / f for the decaying solution f ~ exp(-k r) / r^((D-2)/2)
    auto E_over_f = [&](double r) {return - (mu * mu / (k * k)) * (0.5*(D - 2.) / r + k );};

    //ensures continuity of E (the nodal variable; f may sit near its node at the matching point)
    double f_factor = (E_match / E_over_f(r_match)) * exp(k * r_match) * pow(r_match, 0.5*(D - 2.));

    ProcaFieldState s1, s2, s3, s4;
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

        //fix asymptotic values for f, E
        double r1 = radius_array[j + 1];
        state[j + 1].f = f_factor * exp(-k * r1) / pow(r1, 0.5*(D - 2.));
        state[j + 1].E = E_over_f(r1) * state[j + 1].f;

        if (isnan(state[j].f) || isnan(state[j].X) || isnan(state[j].phi) || isnan(state[j].E))
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

    fill_g_array();

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
        cout << "\nBinding energy is E = " << binding_energy << endl;
     }

     return 1;
}

//rescales lapse by e^(phi_shift) as well as the frequency. Unlike the BS case, f and E (time components) also scale with the lapse
void ProcaStar::rescale_lapse (double phi_shift)
{
    double fac = exp(phi_shift);

    for (ProcaFieldState& s: state)
    {
        s.phi += phi_shift;
        s.f *= fac;
        s.E *= fac;
    }

    omega *= fac;
}

//simply calls find_frequency and fill_asymptotic in one; return 1 iff successful
bool ProcaStar::solve(bool quiet)
{
    find_frequency(quiet);
    bool success = fill_asymptotic(quiet);

    if (!success && g_breakdown_r >= 0. && !quiet)
        cerr << "WARNING: for A_central = " << A_central << ", the algebraic equation V'(Psi) g = -omega E / alpha^2 has no regular root beyond r = "
             << g_breakdown_r << " at the final frequency; no regular solution was found (self-interaction breakdown)" << endl;

    return success;
}

//returns the noether charge associated with the model. Must have computed model + frequency first (polar)
//Q = Omega_{D-2} int r^{D-2} j^t alpha X dr, with j^t = -g E / (alpha^2 X^2)
double ProcaStar::get_noether_charge()
{
    double Q = 0.; //charge to return
    double dr = R / (n_gridpoints - 1);

    if (g_array.size() != state.size())
        fill_g_array();

    //integrate charge over radius
    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = j * dr;
        Q += -pow(r, D - 2.) * dr * g_array[j] * state[j].E * exp(-1. * state[j].phi) / state[j].X;
    }

    Q *= 2. * pow(M_PI, 0.5 * (D - 1.)) / tgamma(0.5 * (D - 1.)); //area of unit (D-2)-sphere
    return Q;
}

//get compactness corresponding to Collodel et. al's C4 definition, defined using the radius within which X is within 0.999 of its Schwarzchild value. Must have computed M!
double ProcaStar::get_c4()
{
    double r4 = R;
    double cutoff_dev = 0.001; //relative deviation from Schwarzschild at which we define PS boundary

    for (int j = n_gridpoints - 3; j > 0; j--)
    {
        double dr = R / (n_gridpoints - 1);
        double r = j * dr;

        double BH_X_sq = 1. /(1. - 2. * M / r);
        double PS_X_sq = state[j].X * state[j].X;

        if (abs((PS_X_sq - BH_X_sq) / BH_X_sq) > cutoff_dev )
        {
            //linearly interpolate to estimate exactly where the relative deviation = cutoff_dev
            double f1 = abs((PS_X_sq - BH_X_sq) / BH_X_sq);

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
double ProcaStar::get_cmax()
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
double ProcaStar::f_RHS(double r, double f)
{
    double dr = R / (n_gridpoints - 1);
    int j_low = floor((r / R) * (n_gridpoints - 1));
    int j0 = bound(j_low - 1, 0, n_gridpoints - 4);

    double X = cubic_interp(r, state[j0].X, state[j0 + 1].X, state[j0 + 2].X, state[j0 + 3].X, j0, dr);

    double r2 = (r == 0.) ? 10e-10 : r;
    return f * (X - 1.) / r2;
}

//returns areal radius corresponding to a given isotropic index
double ProcaStar::r_areal(int j_iso)
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

//fill out arrays for radius, vector components, and conformal factor in isotropic coordinates. Must have solved PS model first!
//f = A_t is unchanged by the radial relabelling; the radial component transforms as A_R = g dr/dR = g psi^2 / X.
void ProcaStar::fill_isotropic_arrays()
{
    double dr = R / (n_gridpoints - 1);

    double f1, f2, f3, f4; //intermediate rk4 values

    if (g_array.size() != state.size())
        fill_g_array();

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

    //fill in f, g, phi, psi arrays
    for (int j = 0; j < n_gridpoints; j++)
    {
        r_areal_array[j] = r_areal(j);

        int j_areal_low = floor(r_areal_array[j] / dr);
        int j_areal_high = ceil (r_areal_array[j] / dr);

        double f, X_val, g_val;

        //cubic interpolation where we won't go OOB in radial array
        if (j_areal_high <= n_gridpoints - 1)
        {
            int j0 = bound(j_areal_low - 1, 0, n_gridpoints - 4);

            phi_iso_array[j] = cubic_interp(r_areal_array[j], state[j0].phi, state[j0 + 1].phi, state[j0 + 2].phi, state[j0 + 3].phi, j0, dr );
            f_iso_array[j] = cubic_interp(r_areal_array[j], state[j0].f, state[j0 + 1].f, state[j0 + 2].f, state[j0 + 3].f, j0, dr );
            g_val = cubic_interp(r_areal_array[j], g_array[j0], g_array[j0 + 1], g_array[j0 + 2], g_array[j0 + 3], j0, dr );
            X_val = cubic_interp(r_areal_array[j], state[j0].X, state[j0 + 1].X, state[j0 + 2].X, state[j0 + 3].X, j0, dr );
            f = cubic_interp(r_areal_array[j], f_array[j0], f_array[j0 + 1], f_array[j0 + 2], f_array[j0 + 3], j0, dr );
        }

        else //extrapolation case; will kick in near boundary
        {
            int N = n_gridpoints;
            double dR = r_areal_array[j] - R;
            phi_iso_array[j] = state[N - 1].phi  + dR * (state[N - 1].phi - state[N - 2].phi ) / dr;
            f_iso_array[j] = state[N - 1].f  + dR * (state[N - 1].f - state[N - 2].f ) / dr;
            g_val = g_array[N - 1]  + dR * (g_array[N - 1] - g_array[N - 2] ) / dr;
            X_val = state[N - 1].X  + dR * (state[N - 1].X - state[N - 2].X ) / dr;
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
        g_iso_array[j] = g_val * psi_iso_array[j] * psi_iso_array[j] / X_val;
    }
}

//writes isotropic values to isotropic_proca.dat
void ProcaStar::write_isotropic()
{
    std::ofstream data_file{ "isotropic_proca.dat" };
    double dr = R / (n_gridpoints - 1);

    if (!data_file)
    {
        std::cerr << "isotropic_proca.dat could not be opened for writing!\n";
        exit(1);
    }

    data_file << "# r r_iso r_areal psi_iso f_iso g_iso phi_iso phi_lap X_lap\n";

    for (int j = 1; j < n_gridpoints - 1; j++)
    {
        double phi_lap = (-2. * phi_iso_array[j] + phi_iso_array[j - 1] + phi_iso_array[j + 1] ) / (dr * dr);
        double X_lap = (-2. * state[j].X + state[j - 1].X + state[j + 1].X ) / (dr * dr);

        data_file << std::setprecision (10) << radius_array[j] << "   " << r_iso_array[j] << "    " << r_areal(j) << "    " << psi_iso_array[j]  <<  "    "
        << f_iso_array[j] << "    " << g_iso_array[j] << "    " << phi_iso_array[j] << "    " << phi_lap <<  "    " << X_lap << std::endl;
    }
}

//constructs a family of PS models with A_central = A_0 + delta_A, A_0 + 2 delta_A, ..., writing to ProcaStars.dat.
//Same columns as BosonStars.dat, plus f0 = f(0) in the gauge alpha(infinity) = 1 (the usual central parameter in the PS literature).
void ProcaStar::cycle_models(int n_stars, double A_0, double delta_A)
{
    read_parameters(0);
    ofstream data_file{"ProcaStars.dat"};
    data_file << "# A_central M r_99 noether_charge binding_energy omega omega_unrescaled phi0 C4 Cmax f0\n";
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
        << "     " << get_cmax() << "     " << state[0].f << endl;

        phi0_values.push_back(state[0].phi);
    }
}

void ProcaStar::broadcast()
{
#ifdef USE_MPI
    // Broadcast scalars that solve() sets.
    double buf[3] = {R, (double)omega, (double)omega_pre_rescale};
    bss_mpi::bcast_doubles(buf, 3);
    R               = buf[0];
    omega           = (long double)buf[1];
    omega_pre_rescale = (long double)buf[2];

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

    // state: vector<ProcaFieldState>, each is 4 contiguous doubles (f, X, phi, E).
    {
        int sz = (int)state.size();
        bss_mpi::bcast_ints(&sz, 1);
        state.resize(sz);
        if (sz > 0)
            bss_mpi::bcast_doubles(reinterpret_cast<double*>(state.data()), sz * 4);
    }

    bcast_vec(g_array);
    bcast_vec(psi_iso_array);
    bcast_vec(f_iso_array);
    bcast_vec(g_iso_array);
    bcast_vec(phi_iso_array);
    bcast_vec(radius_array);
#endif
}

#endif /* PROCASTAR_CPP_ */
