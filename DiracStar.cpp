#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#ifndef DIRACSTAR_CPP_
#define DIRACSTAR_CPP_

#include "DiracStar.h"
#include "mathutils.h"
#include "mpi_helpers.h"
#include <sstream>
#include <algorithm>
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
        fill_parameter(current_line, "perturb = ", perturb, quiet);
        fill_parameter(current_line, "perturb_amp = ", perturb_amp, quiet);
        fill_parameter(current_line, "perturb_spread = ", perturb_spread, quiet);
        fill_parameter(current_line, "perturb_center = ", perturb_center, quiet);

        fill_parameter(current_line, "cycle_only = ", cycle_only, quiet);
        fill_parameter(current_line, "dirac_profile = ", dirac_profile, quiet);
        fill_parameter(current_line, "dirac_profile_file = ", dirac_profile_file, quiet);
        fill_parameter(current_line, "A0 = ", A0, quiet);
        fill_parameter(current_line, "dA = ", dA, quiet);
        fill_parameter(current_line, "n_stars = ", n_stars, quiet);
    }

    check_supported();
}

//the Dirac star solver is implemented in D = 4 only; sets the self-interaction potential
void DiracStar::check_supported()
{
    if (D != 4.)
    {
        cerr << "ERROR: Dirac stars are only implemented for D = 4 (the number of spinors and their angular structure depend on D)" << endl;
        exit(1);
    }
    if (solitonic && sigma <= 0.)
    {
        cerr << "ERROR: the degenerate-vacuum Dirac potential needs sigma > 0" << endl;
        exit(1);
    }
    pot = DiracPotential{mu, lambda, solitonic, sigma};
}

//scalar density S = sum_m psibar_m psi_m = (G^2 - F^2) / (2 pi alpha r^2); at r = 0, S = A_central^2 / (2 pi)
double DiracStar::scalar_density(double r, const DiracFieldState& s)
{
    if (r == 0.)
        return A_central * A_central / (2. * M_PI);
    double Gr = s.G / r, Fr = s.F / r;
    return (Gr * Gr - Fr * Fr) / (2. * M_PI * exp(s.phi));
}

//number density n, energy density rho and radial pressure p_r of the two spinors at radius r (orthonormal frame);
//at r = 0 the regular limits G/r -> g0 sqrt(alpha0), F/r -> 0, F/r^2 -> (mu_eff - omega/alpha0) g0 sqrt(alpha0) / 3 are used.
//With the self-interaction, rho = (omega/alpha) n - P(S) and p_r = [...] - U(S) (see DiracStar.h)
void DiracStar::matter_terms(double r, long double frequency, const DiracFieldState& s, double& rho, double& p_r, double& n)
{
    double alpha = exp(s.phi);
    double w_a = static_cast<double>(frequency) / alpha;
    double S = scalar_density(r, s);
    double mu_eff = pot.mu_eff(S);

    double Gr, Fr, FG_r3; //G/r, F/r and F G / r^3
    if (r == 0.)
    {
        Gr = A_central * sqrt(alpha);
        Fr = 0.;
        FG_r3 = (mu_eff - w_a) * Gr * Gr / 3.;
    }
    else
    {
        Gr = s.G / r;
        Fr = s.F / r;
        FG_r3 = Fr * Gr / r;
    }

    double norm = 1. / (2. * M_PI * alpha);
    n = norm * (Fr * Fr + Gr * Gr);
    rho = w_a * n - pot.P(S);
    p_r = norm * (w_a * (Fr * Fr + Gr * Gr) - mu * (Gr * Gr - Fr * Fr) + 2. * FG_r3) - pot.U(S);
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

    double mu_eff = pot.mu_eff(scalar_density(r, s));
    double dF = s.X * (-Fr + (mu_eff - w_a) * s.G);
    double dG = s.X * (Gr + (mu_eff + w_a) * s.F);

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

        //with the self-interaction a diverging solution runs away within a few grid points (mu_eff grows with the field),
        //and after the metric becomes singular G can change sign, giving spurious zero crossings. Stop as soon as the
        //physical amplitude g = G / (sqrt(alpha) r) or f exceeds 3 A_central (bound states have |f|, |g| <~ A_central).
        //Not applied to the free field, whose shooting is unchanged.
        if (!pot.free_field())
        {
            double bound_amp = 3. * A_central * (j + 1) * dr * exp(0.5 * state[j + 1].phi);
            if (abs(state[j + 1].G) > bound_amp || abs(state[j + 1].F) > bound_amp)
            {
                blowup_point = j + 1;
                return;
            }
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

//simply calls find_frequency and fill_asymptotic in one; return 1 iff successful. With dirac_profile = 1 the star is read
//from a relaxation profile instead (read_profile_solution)
bool DiracStar::solve(bool quiet)
{
    if (dirac_profile)
        return read_profile_solution(quiet);
    find_frequency(quiet);
    return fill_asymptotic(quiet);
}

//Loads a star computed by the relaxation solver (relaxation/export_profile.py) onto this grid (R, n_gridpoints), in place of
//shooting. The file has '#' header lines, among them "# omega = ...", "# sigma = ...", "# lambda = ...", "# mu = ...",
//"# solitonic = ...", followed by columns  r f g m phi  (physical spinor amplitudes f, g = F, G / (sqrt(alpha) r), aspect mass
//m, phi = ln alpha with alpha(infinity) = 1) on an increasing radial grid starting at r = 0. The profile is interpolated with
//4-point Lagrange polynomials; beyond its last radius the exterior is Schwarzschild (m = M) with f, g continued as
//exp(-k r), k = X sqrt(mu^2 - omega^2/alpha^2). The potential parameters in the file must match those of BSParams.par.
bool DiracStar::read_profile_solution(bool quiet)
{
    ifstream in{dirac_profile_file};
    if (!in)
    {
        cerr << "ERROR: could not open Dirac profile " << dirac_profile_file << endl;
        exit(1);
    }
    vector<double> rr, ff, gg, mm, pp;
    double w_file = -1., sig_file = sigma, lam_file = lambda, mu_file = mu, sol_file = solitonic;
    string line;
    while (getline(in, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '#')
        {
            auto grab = [&](const string& key, double& v) {
                size_t p = line.find(key);
                if (p != string::npos) v = stod(line.substr(p + key.size()));
            };
            grab("omega = ", w_file); grab("sigma = ", sig_file); grab("lambda = ", lam_file); grab("mu = ", mu_file);
            grab("solitonic = ", sol_file);
            continue;
        }
        stringstream ss(line);
        double a, b, c, d, e;
        if (ss >> a >> b >> c >> d >> e)
        {
            rr.push_back(a); ff.push_back(b); gg.push_back(c); mm.push_back(d); pp.push_back(e);
        }
    }
    if (w_file <= 0. || rr.size() < 8 || rr[0] != 0.)
    {
        cerr << "ERROR: Dirac profile " << dirac_profile_file << " needs '# omega = ...' and data from r = 0" << endl;
        exit(1);
    }
    if (abs(sig_file - sigma) > 1e-12 * sigma || abs(lam_file - lambda) > 1e-12 * (1. + abs(lambda)) || abs(mu_file - mu) > 1e-12
        || (bool)sol_file != solitonic)
    {
        cerr << "ERROR: the potential of Dirac profile " << dirac_profile_file << " (sigma " << sig_file << ", lambda " << lam_file
             << ", mu " << mu_file << ", solitonic " << sol_file << ") differs from BSParams.par" << endl;
        exit(1);
    }

    const int Nf = rr.size();
    auto interp = [&](const vector<double>& y, double r) {
        int j = upper_bound(rr.begin(), rr.end(), r) - rr.begin() - 2;
        j = max(0, min(j, Nf - 4));
        double res = 0.;
        for (int a = 0; a < 4; a++)
        {
            double L = 1.;
            for (int b = 0; b < 4; b++)
                if (b != a)
                    L *= (r - rr[j + b]) / (rr[j + a] - rr[j + b]);
            res += L * y[j + a];
        }
        return res;
    };

    omega = w_file;
    omega_pre_rescale = w_file;
    A_central = gg[0];
    state.assign(n_gridpoints, DiracFieldState{});
    radius_array.resize(n_gridpoints);
    blowup_point = n_gridpoints;
    const double dr = R / (n_gridpoints - 1);
    const double r_last = rr[Nf - 1], M_file = mm[Nf - 1];
    const double f_last = ff[Nf - 1], g_last = gg[Nf - 1];
    const double X_last = 1. / sqrt(1. - 2. * M_file / r_last);
    const double k = X_last * sqrt(max(mu * mu - pow(w_file / exp(pp[Nf - 1]), 2), 1e-12));
    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = j * dr;
        radius_array[j] = r;
        double f, g, m_, phi;
        if (r <= r_last)
        {
            f = interp(ff, r); g = interp(gg, r); m_ = interp(mm, r); phi = interp(pp, r);
        }
        else
        {
            double dec = exp(-k * (r - r_last));
            f = f_last * dec; g = g_last * dec; m_ = M_file; phi = 0.5 * log(1. - 2. * M_file / r);
        }
        if (r == 0.)
            m_ = 0.;
        double sa = exp(0.5 * phi);
        state[j] = DiracFieldState{sa * r * f, sa * r * g, (r == 0.) ? 1. : 1. / sqrt(1. - 2. * m_ / r), phi};
    }

    M = m(n_gridpoints - 1);
    int j = 0;
    while (j < n_gridpoints - 1 && m(j) < 0.99 * M)
    {
        if (m(j + 1) != m(j))
            r_99 = dr * j + dr * (0.99 * M - m(j)) / (m(j + 1) - m(j));
        j++;
    }
    noether_charge = get_noether_charge();
    binding_energy = M - mu * noether_charge;
    compactness = M / r_99;
    if (!quiet)
    {
        cout << "\nRead Dirac star from " << dirac_profile_file << " (" << Nf << " points to r = " << r_last << "): omega = " << omega
             << ", A_central = " << A_central << "\nMass is M =  " << M << "\nNoether charge (particle number of both spinors) is N = "
             << noether_charge << "\nBinding energy is E = " << binding_energy << "\nr_99 = " << r_99 << ", C = M / r_99 = "
             << compactness << ", Cmax = " << get_cmax() << endl;
    }
    return 1;
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

//Perturbs a solved star by scaling the spinor amplitudes, (u, v) -> (1 + a exp(-(r - center)^2 / k^2)) (u, v), then re-solves
//X from the Hamiltonian constraint and phi from the polar slicing condition with u and v held fixed, as
//ProcaStar::add_perturbation does for the Proca field. At t = 0 u and v are real, so the momentum density
//Im(ubar u' + vbar v') / (2 pi) vanishes and K_ij = 0 remains consistent. Since the perturbed data are not stationary, the
//energy density and radial stress are the 3+1 expressions in terms of u, v and their radial derivatives (ComplexSpinorField):
//    S_rr^ = (v u' - u v') / (2 pi X) + P(S),   S_thth^ = -u v / (2 pi r) + P(S),
//    rho = mu S + U(S) + [kinetic parts of S_rr^ + 2 S_thth^],   S = (u^2 - v^2) / (2 pi),
//which for the unperturbed star agree with rho = omega n / alpha and p_r (a = 0 reproduces the star). a > 0 makes the star
//denser. The lapse is rescaled to alpha(infinity) = 1 afterwards, and M, N, E_b and r_99 are recomputed.
void DiracStar::add_perturbation(double a, double k, double center, bool quiet)
{
    int N = state.size();
    double dr = R / (N - 1);

    //perturbed physical spinor amplitudes on the grid
    vector<double> u(N), v(N);
    for (int j = 0; j < N; j++)
    {
        double r = j * dr;
        double sa = exp(0.5 * state[j].phi);
        double fac = 1. + a * exp(-pow(r - center, 2.) / (k * k));
        u[j] = fac * ((j == 0) ? A_central : state[j].G / (sa * r));
        v[j] = fac * ((j == 0) ? 0. : state[j].F / (sa * r));
    }

    //4th-order radial derivatives, using the parity across r = 0 (u even, v odd) and one-sided stencils at the outer edge
    auto deriv = [&](const vector<double>& y, bool odd)
    {
        vector<double> d(N);
        auto at = [&](int j) { return (j >= 0) ? y[j] : (odd ? -y[-j] : y[-j]); };
        for (int j = 0; j < N; j++)
        {
            if (j <= N - 3)
                d[j] = (at(j - 2) - 8. * at(j - 1) + 8. * at(j + 1) - y[j + 2]) / (12. * dr);
            else
                d[j] = (3. * y[j - 4] - 16. * y[j - 3] + 36. * y[j - 2] - 48. * y[j - 1] + 25. * y[j]) / (12. * dr);
        }
        return d;
    };
    vector<double> du = deriv(u, 0), dv = deriv(v, 1);

    auto interp = [&](const vector<double>& y, double r)
    {
        int j0 = bound((int)floor(r / dr) - 1, 0, N - 4);
        return cubic_interp(r, y[j0], y[j0 + 1], y[j0 + 2], y[j0 + 3], j0, dr);
    };

    //RHS of the Hamiltonian constraint (X') and polar slicing condition (phi') for given u, v, u', v'
    auto rhs = [&](double r, double X, double uu, double vv, double duu, double dvv, double& dX, double& dPhi)
    {
        double S = (uu * uu - vv * vv) / (2. * M_PI);
        double S_rr = (vv * duu - uu * dvv) / (2. * M_PI * X);
        double v_r = (r == 0.) ? dvv : vv / r;
        double S_tt = -uu * v_r / (2. * M_PI);
        double rho = mu * S + pot.U(S) + S_rr + 2. * S_tt;
        S_rr += pot.P(S);
        double T1 = (r == 0.) ? 0. : 0.5 * (X * X - 1.) / r;
        dX = X * (4. * M_PI * G * r * X * X * rho - T1);
        dPhi = T1 + 4. * M_PI * G * r * X * X * S_rr;
    };

    //re-integrate X, phi outwards from the centre (X(0) = 1, phi(0) unchanged)
    state[0].X = 1.;
    for (int j = 0; j < N - 1; j++)
    {
        double r = j * dr, rm = r + 0.5 * dr;
        double X = state[j].X, phi = state[j].phi;
        double um = interp(u, rm), vm = interp(v, rm), dum = interp(du, rm), dvm = interp(dv, rm);
        double k1X, k1p, k2X, k2p, k3X, k3p, k4X, k4p;

        rhs(r, X, u[j], v[j], du[j], dv[j], k1X, k1p);
        rhs(rm, X + 0.5 * dr * k1X, um, vm, dum, dvm, k2X, k2p);
        rhs(rm, X + 0.5 * dr * k2X, um, vm, dum, dvm, k3X, k3p);
        rhs(r + dr, X + dr * k3X, u[j + 1], v[j + 1], du[j + 1], dv[j + 1], k4X, k4p);

        state[j + 1].X = X + (dr / 6.) * (k1X + 2. * k2X + 2. * k3X + k4X);
        state[j + 1].phi = phi + (dr / 6.) * (k1p + 2. * k2p + 2. * k3p + k4p);

        if (isnan(state[j + 1].X) || isnan(state[j + 1].phi))
        {
            cerr << "ERROR: constraint solve for perturbed Dirac star failed at r = " << r + dr << endl;
            exit(1);
        }
    }

    //store the perturbed spinor in the F, G variables (which carry the factor sqrt(alpha) r)
    for (int j = 0; j < N; j++)
    {
        double r = j * dr;
        double sa = exp(0.5 * state[j].phi);
        state[j].G = sa * r * u[j];
        state[j].F = sa * r * v[j];
    }
    A_central = u[0];

    //enforce alpha(infinity) = 1 (the physical u, v are unaffected)
    double phi_shift = -log(state[N - 1].X) - state[N - 1].phi;
    rescale_lapse(phi_shift);

    double M_old = M, N_old = noether_charge;
    M = m(N - 1);
    noether_charge = get_noether_charge();
    binding_energy = M - mu * noether_charge;

    int j = 0;
    while (j < N - 1 && m(j) < 0.99 * M)
    {
        if (m(j + 1) != m(j))
            r_99 = dr * j + dr * (0.99 * M - m(j)) / (m(j + 1) - m(j));
        j++;
    }
    compactness = M / r_99;

    if (!quiet)
        cout << "\nPerturbed Dirac star with amplitude " << a << ", width " << k << ", centre " << center
             << ": M = " << M_old << " -> " << M << ", N = " << N_old << " -> " << noether_charge
             << ", binding energy = " << binding_energy << endl;
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
