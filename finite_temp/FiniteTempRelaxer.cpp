#include "FiniteTempRelaxer.h"
#include "../bose_gas.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include "mathutils.h"
#include "seed_shooter.h"

namespace
{
constexpr double kPi = 3.1415926535897932384626433832795;

double interp_linear_sorted(const std::vector<double>& x, const std::vector<double>& y, double xi)
{
    if (x.empty() || y.empty() || x.size() != y.size())
    {
        return 0.0;
    }

    if (xi <= x.front())
    {
        return y.front();
    }
    if (xi >= x.back())
    {
        return y.back();
    }

    size_t hi = 1;
    while (hi < x.size() && x[hi] < xi)
    {
        ++hi;
    }
    const size_t lo = hi - 1;
    const double dx = x[hi] - x[lo];
    if (dx <= 0.0)
    {
        return y[lo];
    }
    const double t = (xi - x[lo]) / dx;
    return y[lo] * (1.0 - t) + y[hi] * t;
}
}

void FiniteTempRelaxer::read_parameters(bool quiet)
{
    read_parameters("FiniteTempRelaxer.par", quiet);
}

void FiniteTempRelaxer::read_parameters(const std::string& filename, bool quiet)
{
    std::ifstream params{filename};
    if (!params)
    {
        std::cerr << "Could not open " << filename << "\n";
        std::abort();
    }

    std::string current_line;
    while (std::getline(params, current_line))
    {
        fill_parameter(current_line, "n_gridpoints = ", p_.n_gridpoints, quiet);
        fill_parameter(current_line, "rbar_min = ", p_.rbar_min, quiet);
        fill_parameter(current_line, "rbar_max = ", p_.rbar_max, quiet);

        fill_parameter(current_line, "G = ", p_.G, quiet);
        fill_parameter(current_line, "mu = ", p_.mu, quiet);
        fill_parameter(current_line, "lambda = ", p_.lambda, quiet);
        fill_parameter(current_line, "solitonic = ", p_.solitonic, quiet);
        fill_parameter(current_line, "sigma_scale = ", p_.sigma_scale, quiet);

        fill_parameter(current_line, "sigma_central = ", p_.sigma_central, quiet);
        fill_parameter(current_line, "alpha_central = ", p_.alpha_central, quiet);
        fill_parameter(current_line, "mu_chem_infty = ", p_.mu_chem_infty, quiet);
        fill_parameter(current_line, "use_shooting_seed = ", p_.use_shooting_seed, quiet);
        fill_parameter(current_line, "seed_R = ", p_.seed_R, quiet);
        fill_parameter(current_line, "seed_n_gridpoints = ", p_.seed_n_gridpoints, quiet);
        fill_parameter(current_line, "seed_eigen = ", p_.seed_eigen, quiet);
        fill_parameter(current_line, "seed_frequency_guess = ", p_.seed_frequency_guess, quiet);
        fill_parameter(current_line, "seed_freq_epsilon = ", p_.seed_freq_epsilon, quiet);
        fill_parameter(current_line, "seed_r_match_fac = ", p_.seed_r_match_fac, quiet);

        fill_parameter(current_line, "temperature_target = ", p_.temperature_target, quiet);
        fill_parameter(current_line, "temperature_step = ", p_.temperature_step, quiet);
        fill_parameter(current_line, "force_continue = ", p_.force_continue, quiet);
        fill_parameter(current_line, "first_thermalized_only = ", p_.first_thermalized_only, quiet);
        fill_parameter(current_line, "use_suppression_factor = ", p_.use_suppression_factor, quiet);
        fill_parameter(current_line, "decay_exponent = ", p_.decay_exponent, quiet);
        fill_parameter(current_line, "cph_epsilon = ", p_.cph_epsilon, quiet);
        fill_parameter(current_line, "pressure_interp_power = ", p_.pressure_interp_power, quiet);
        fill_parameter(current_line, "crossover_mode = ", p_.crossover_mode, quiet);
        fill_parameter(current_line, "gas_statistics = ", p_.gas_statistics, quiet);
        fill_parameter(current_line, "sigma_phonon = ", p_.sigma_phonon, quiet);
        fill_parameter(current_line, "sigma_gas = ", p_.sigma_gas, quiet);
        fill_parameter(current_line, "crossover_exponent = ", p_.crossover_exponent, quiet);

        fill_parameter(current_line, "omega_target = ", p_.omega_target, quiet);
        fill_parameter(current_line, "omega_step = ", p_.omega_step, quiet);
        fill_parameter(current_line, "charge_step = ", p_.charge_step, quiet);
        fill_parameter(current_line, "sigma_jump_on_fail = ", p_.sigma_jump_on_fail, quiet);
        fill_parameter(current_line, "increment_omega = ", p_.increment_omega, quiet);
        fill_parameter(current_line, "num_models = ", p_.num_models, quiet);

        fill_param_array(current_line, "hi_sigma_vals = ", p_.hi_sigma_vals, quiet);
        fill_param_array(current_line, "lo_sigma_vals = ", p_.lo_sigma_vals, quiet);
        fill_parameter(current_line, "resolution_factor = ", p_.resolution_factor, quiet);

        fill_parameter(current_line, "max_newton_iterations = ", p_.max_newton_iterations, quiet);
        fill_parameter(current_line, "newton_tolerance = ", p_.newton_tolerance, quiet);
        fill_parameter(current_line, "jacobian_fd_step = ", p_.jacobian_fd_step, quiet);
        fill_parameter(current_line, "line_search_shrink = ", p_.line_search_shrink, quiet);
        fill_parameter(current_line, "line_search_min_step = ", p_.line_search_min_step, quiet);

        fill_parameter(current_line, "output_profile = ", p_.output_profile, quiet);
        fill_parameter(current_line, "save_model_profiles = ", p_.save_model_profiles, quiet);
    }
}

void FiniteTempRelaxer::build_grid()
{
    const int n = std::max(3, p_.n_gridpoints);
    p_.n_gridpoints = n;

    state_.rbar.resize(static_cast<size_t>(n));
    state_.sigma.resize(static_cast<size_t>(n));
    state_.eta.resize(static_cast<size_t>(n));
    state_.alpha.resize(static_cast<size_t>(n));
    state_.Xs.resize(static_cast<size_t>(n));
    state_.omega.resize(static_cast<size_t>(n));
    state_.rho.resize(static_cast<size_t>(n));

    const double drbar = (p_.rbar_max - p_.rbar_min) / static_cast<double>(n - 1);
    for (int i = 0; i < n; ++i)
    {
        state_.rbar[static_cast<size_t>(i)] = p_.rbar_min + drbar * static_cast<double>(i);
    }
    
    // Validate the compactification scheme's inverse transformation
    compactification_->validate_inverse(1e-12);
}

double FiniteTempRelaxer::areal_radius(double rbar) const
{
    return compactification_->r_of_rbar(rbar);
}

double FiniteTempRelaxer::potential(double sigma) const
{
    return potential_of_sigma_sq(sigma * sigma);
}

double FiniteTempRelaxer::potential_of_sigma_sq(double sigma_sq) const
{
    const double sigma2 = std::max(0.0, sigma_sq);
    if (!p_.solitonic)
    {
        return p_.mu * p_.mu * sigma2 + p_.lambda * sigma2 * sigma2;
    }

    const double scale2 = std::max(1e-12, p_.sigma_scale * p_.sigma_scale);
    const double shape = 1.0 - 2.0 * sigma2 / scale2;
    return p_.mu * p_.mu * sigma2 * shape * shape + p_.lambda * sigma2 * sigma2;
}

double FiniteTempRelaxer::d_potential_d_sigma2(double sigma) const
{
    const double sigma2 = sigma * sigma;
    if (!p_.solitonic)
    {
        return p_.mu * p_.mu + 2.0 * p_.lambda * sigma2;
    }

    const double inv_scale2 = 1.0 / (std::max(1e-12, p_.sigma_scale) * std::max(1e-12, p_.sigma_scale));
    const double scaled2 = sigma2 * inv_scale2;
    return p_.mu * p_.mu - 8.0 * p_.mu * p_.mu * scaled2 + 12.0 * p_.mu * p_.mu * scaled2 * scaled2 + 2.0 * p_.lambda * sigma2;
}

double FiniteTempRelaxer::ddV(double sigma_sq) const
{
    if (!p_.solitonic)
    {
        return 2.0 * p_.lambda;
    }

    const double scale2 = std::max(1e-12, p_.sigma_scale * p_.sigma_scale);
    const double scale4 = scale2 * scale2;
    return -8.0 * p_.mu * p_.mu / scale2 + 24.0 * p_.mu * p_.mu * sigma_sq / scale4 + 2.0 * p_.lambda;
}

double FiniteTempRelaxer::cph_sq(double sigma_sq, double Xi) const
{
    const double sig2 = std::max(0.0, sigma_sq);
    const double vpp = ddV(sig2);
    const double num = sig2 * vpp + std::max(0.0, p_.cph_epsilon); //exp(-sig2 / std::max(0.0, p_.cph_epsilon)); //+ std::max(0.0, p_.cph_epsilon);
    const double den = sig2 * vpp - 2.0 * Xi;
    const double safe_den = (std::abs(den) < 1e-12) ? ((den >= 0.0) ? 1e-12 : -1e-12) : den;
    return std::max(1e-12, num / safe_den);
}

double FiniteTempRelaxer::scalar_pressure_term(double sigma_sq, double X, double Xi) const
{
    const double sig2 = std::max(0.0, sigma_sq);
    const double Xpos = std::max(0.0, X);
    const double Xi_safe = (std::abs(Xi) < 1e-12) ? ((Xi >= 0.0) ? 1e-12 : -1e-12) : Xi;
    const double V = potential_of_sigma_sq(sig2);
    return -0.5 * (Xpos + sig2 * Xi_safe + V);
}

double FiniteTempRelaxer::thermal_pressure_term(double T, double mu_chem, double sigma_sq, double X, double Xi) const
{
    const double sig2 = std::max(0.0, sigma_sq);
    const double cph2 = cph_sq(sig2, Xi);
    const double cph = std::sqrt(std::max(1e-12, cph2));
    const double safe_Xi = (std::abs(Xi) < 1e-12) ? ((Xi >= 0.0) ? 1e-12 : -1e-12) : Xi;
    const double bracket = 1.0 + (1.0 - cph2) * (mu_chem * mu_chem) / safe_Xi;
    const double bracket_safe = (std::abs(bracket) < 1e-12) ? ((bracket >= 0.0) ? 1e-12 : -1e-12) : bracket;
    (void)X;
    double P_phonon = (kPi * kPi / 90.0) * cph * std::pow(std::max(0.0, T), 4.0)
        * std::pow(bracket_safe, -2.0);

    if (p_.crossover_mode == 1)
    {
        const double kh = 0.5 * p_.crossover_exponent;
        const double x_ph = std::pow(sig2 / std::max(1e-300, p_.sigma_phonon * p_.sigma_phonon), kh);
        const double x_g = std::pow(sig2 / std::max(1e-300, p_.sigma_gas * p_.sigma_gas), kh);
        const double w = -std::expm1(-x_ph);
        double lnF, kap;
        bose_gas::lnF_kappa((mu_chem - p_.mu) / std::max(1e-12, T), p_.gas_statistics, lnF, kap);
        const double lnPg = 1.5 * std::log(p_.mu / (2. * kPi)) + 2.5 * std::log(std::max(1e-300, T)) + lnF;
        return ((w > 0.0) ? w * P_phonon : 0.0) + std::exp(lnPg - x_g);
    }

    double lnF, kap;
    bose_gas::lnF_kappa((mu_chem - p_.mu) / std::max(1e-12, T), p_.gas_statistics, lnF, kap);
    double P_gas = pow(p_.mu / (2. * kPi), 1.5) * pow(std::max(0.0, T), 2.5) * exp(lnF);
    const double q = std::max(1e-12, p_.pressure_interp_power);

    double P_total = pow(pow(P_phonon, -q) + pow(P_gas, -q), -1.0 / q);
    
    return P_total;
}

double FiniteTempRelaxer::compute_suppression_factor(double sigma_sq) const
{
    // Suppress thermal contribution as sigma_sq -> 0 to avoid non-vanishing
    // asymptotic thermal energy density in asymptotically-flat BS solutions.
    const int n = std::max(1, p_.decay_exponent);
    const double sigma0 = state_.sigma.empty() ? 0.0 : state_.sigma.front();
    const double sigma0_sq = std::max(1e-16, sigma0 * sigma0);
    const double x = std::min(std::max(0.0, sigma_sq / sigma0_sq), 1.0);
    const double inner = 1.0 - std::pow(1.0 - x, static_cast<double>(n));
    double suppression_factor = std::pow(std::max(0.0, inner), static_cast<double>(n));
    
    //suppression_factor = pow(x, 4.0);

    double f = 1. / (1. - x) - 1. / x;
    double a = 1.0;
    double k = 1.0; 

    double g = 1.0;

    if (x < 1e-12)
        g = 0.0;

    else if (std::abs(x - 1.0) > 1e-12)
        g = 0.5 * (1.0 + std::tanh(k * (f + a)));

    //suppression_factor = g;

    if (!p_.use_suppression_factor)
    {
        suppression_factor = 1.0;
    }

    return suppression_factor;
}

double FiniteTempRelaxer::pressure_Pthermal(double T, double mu_chem, double sigma_sq, double X, double Xi) const
{
    // Pthermal Eq. (notes): P = scalar piece + suppression * thermal piece
    const double suppression_factor = compute_suppression_factor(sigma_sq);
    if (T < 1e-12)
    {
        return scalar_pressure_term(sigma_sq, X, Xi);
    }

    return scalar_pressure_term(sigma_sq, X, Xi)
        + suppression_factor * thermal_pressure_term(T, mu_chem, sigma_sq, X, Xi);
}

FiniteTempRelaxer::ThermoConjugates FiniteTempRelaxer::conjugates_from_pressure(
    double T,
    double mu_chem,
    double sigma_sq,
    double X,
    double Xi) const
{
    ThermoConjugates out;

    const double fd = std::max(1e-7, p_.jacobian_fd_step);
    const double dT = fd * std::max(1.0, std::abs(T));
    const double dmu = fd * std::max(1.0, std::abs(mu_chem));
    const double ds2 = fd * std::max(1.0, std::abs(sigma_sq));
    const double dXi = fd * std::max(1.0, std::abs(Xi));

    const double pTp = pressure_Pthermal(T + dT, mu_chem, sigma_sq, X, Xi);
    const double pTm = pressure_Pthermal(T - dT, mu_chem, sigma_sq, X, Xi);
    out.s = (pTp - pTm) / (2.0 * dT);

    const double pMup = pressure_Pthermal(T, mu_chem + dmu, sigma_sq, X, Xi);
    const double pMum = pressure_Pthermal(T, mu_chem - dmu, sigma_sq, X, Xi);
    out.n = (pMup - pMum) / (2.0 * dmu);

    const double pS2p = pressure_Pthermal(T, mu_chem, sigma_sq + ds2, X, Xi);
    const double pS2m = pressure_Pthermal(T, mu_chem, std::max(0.0, sigma_sq - ds2), X, Xi);
    out.F = 2.0 * (pS2p - pS2m) / (2.0 * ds2);

    const double pXip = pressure_Pthermal(T, mu_chem, sigma_sq, X, Xi + dXi);
    const double pXim = pressure_Pthermal(T, mu_chem, sigma_sq, X, Xi - dXi);
    out.Upsilon = -2.0 * (pXip - pXim) / (2.0 * dXi);

    return out;
}

double FiniteTempRelaxer::local_temperature(double alpha) const
{
    const double lapse = std::max(1e-12, alpha);
    return current_temperature_ / lapse;
}

void FiniteTempRelaxer::initialize_seed()
{
    const int n = p_.n_gridpoints;

    bool seeded_from_shooting = false;
    if (p_.use_shooting_seed)
    {
        seedbs::BosonStarSeedShooter shooter;
        shooter.D = 4.0;
        shooter.mu = p_.mu;
        shooter.lambda = p_.lambda;
        shooter.solitonic = p_.solitonic;
        shooter.G = p_.G;
        shooter.sigma_scale = p_.sigma_scale;
        shooter.eigen = p_.seed_eigen;
        shooter.alpha_central = p_.alpha_central;
        shooter.A_central = p_.sigma_central;
        shooter.frequency_guess = p_.seed_frequency_guess;
        shooter.r_match_fac = p_.seed_r_match_fac;
        shooter.R = p_.seed_R;
        shooter.n_gridpoints = p_.seed_n_gridpoints;
        shooter.freq_epsilon = p_.seed_freq_epsilon;

        if (shooter.solve(true))
        {
            shooter.write_field("seed_bs.dat");

            std::vector<double> sigma_seed(shooter.state.size());
            std::vector<double> eta_seed(shooter.state.size());
            std::vector<double> alpha_seed(shooter.state.size());
            std::vector<double> X_seed(shooter.state.size());

            for (size_t i = 0; i < shooter.state.size(); ++i)
            {
                sigma_seed[i] = shooter.state[i].A;
                eta_seed[i] = shooter.state[i].eta;
                alpha_seed[i] = std::exp(shooter.state[i].phi);
                X_seed[i] = shooter.state[i].X;
            }

            for (int i = 0; i < n; ++i)
            {
                const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
                state_.sigma[static_cast<size_t>(i)] = interp_linear_sorted(shooter.radius_array, sigma_seed, r);
                state_.eta[static_cast<size_t>(i)] = interp_linear_sorted(shooter.radius_array, eta_seed, r);
                state_.alpha[static_cast<size_t>(i)] = std::max(1e-12, interp_linear_sorted(shooter.radius_array, alpha_seed, r));
                state_.Xs[static_cast<size_t>(i)] = std::max(1e-12, interp_linear_sorted(shooter.radius_array, X_seed, r));
                state_.omega[static_cast<size_t>(i)] = static_cast<double>(shooter.omega);
                state_.rho[static_cast<size_t>(i)] = 1.0;
            }
            p_.mu_chem_infty = static_cast<double>(shooter.omega);
            seeded_from_shooting = true;
        }
    }

    if (!seeded_from_shooting)
    {
        const double r0 = std::max(1e-6, 0.20 * areal_radius(0.5));
        for (int i = 0; i < n; ++i)
        {
            const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
            const double profile = std::exp(-(r * r) / (r0 * r0));
            state_.sigma[static_cast<size_t>(i)] = p_.sigma_central * profile;
            state_.eta[static_cast<size_t>(i)] = 0.0;
            state_.alpha[static_cast<size_t>(i)] = p_.alpha_central;
            state_.Xs[static_cast<size_t>(i)] = 1.0;
            state_.omega[static_cast<size_t>(i)] = p_.mu_chem_infty;
            state_.rho[static_cast<size_t>(i)] = 1.0;
        }
    }

    current_temperature_ = 0.0;
    enforce_boundary_values();
}

bool FiniteTempRelaxer::shoot_for_sigma(double sigma0)
{
    seedbs::BosonStarSeedShooter shooter;
    shooter.D = 4.0;
    shooter.mu = p_.mu;
    shooter.lambda = p_.lambda;
    shooter.solitonic = p_.solitonic;
    shooter.G = p_.G;
    shooter.sigma_scale = p_.sigma_scale;
    shooter.eigen = p_.seed_eigen;
    shooter.alpha_central = p_.alpha_central;
    shooter.A_central = sigma0;
    shooter.frequency_guess = p_.seed_frequency_guess;
    shooter.r_match_fac = p_.seed_r_match_fac;
    shooter.R = p_.seed_R;
    shooter.n_gridpoints = p_.seed_n_gridpoints;
    shooter.freq_epsilon = p_.seed_freq_epsilon;

    if (!shooter.solve(true))
    {
        return false;
    }

    shooter.write_field("seed_bs.dat");

    std::vector<double> sigma_seed(shooter.state.size());
    std::vector<double> eta_seed(shooter.state.size());
    std::vector<double> alpha_seed(shooter.state.size());
    std::vector<double> X_seed(shooter.state.size());

    for (size_t i = 0; i < shooter.state.size(); ++i)
    {
        sigma_seed[i] = shooter.state[i].A;
        eta_seed[i] = shooter.state[i].eta;
        alpha_seed[i] = std::exp(shooter.state[i].phi);
        X_seed[i] = shooter.state[i].X;
    }

    const int n = p_.n_gridpoints;
    for (int i = 0; i < n; ++i)
    {
        const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
        state_.sigma[static_cast<size_t>(i)] = interp_linear_sorted(shooter.radius_array, sigma_seed, r);
        state_.eta[static_cast<size_t>(i)] = interp_linear_sorted(shooter.radius_array, eta_seed, r);
        state_.alpha[static_cast<size_t>(i)] = std::max(1e-12, interp_linear_sorted(shooter.radius_array, alpha_seed, r));
        state_.Xs[static_cast<size_t>(i)] = std::max(1e-12, interp_linear_sorted(shooter.radius_array, X_seed, r));
        state_.omega[static_cast<size_t>(i)] = static_cast<double>(shooter.omega);
        state_.rho[static_cast<size_t>(i)] = 1.0;
    }
    p_.mu_chem_infty = static_cast<double>(shooter.omega);
    enforce_boundary_values();
    return true;
}

void FiniteTempRelaxer::enforce_boundary_values()
{
    const int n = p_.n_gridpoints;
    if (n < 2)
    {
        return;
    }

    state_.eta[0] = 0.0;
    state_.Xs[0] = 1.0;
    state_.omega[0] = state_.omega[1];
    state_.rho[0] = state_.rho[1];
    if (p_.increment_omega)
    {
        state_.sigma[0] = state_.sigma[1];
        state_.omega[static_cast<size_t>(n - 1)] = state_.omega[static_cast<size_t>(n - 2)];
    }

    const int last = n - 1;
    state_.sigma[static_cast<size_t>(last)] = 0.0;
    state_.eta[static_cast<size_t>(last)] = 0.0;
    state_.alpha[static_cast<size_t>(last)] = 1.0;
    state_.Xs[static_cast<size_t>(last)] = 1.0;
    if (!p_.increment_omega)
    {
        state_.omega[static_cast<size_t>(last)] = p_.mu_chem_infty;
        state_.rho[static_cast<size_t>(last)] = 1.0;
    }
}

std::vector<double> FiniteTempRelaxer::pack_state() const
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    const size_t fields = active_field_count();
    std::vector<double> x(fields * n);
    for (size_t i = 0; i < n; ++i)
    {
        x[i] = state_.sigma[i];
        x[n + i] = state_.eta[i];
        x[2 * n + i] = state_.alpha[i];
        x[3 * n + i] = state_.Xs[i];
        x[4 * n + i] = state_.omega[i];
        if (p_.increment_omega)
        {
            x[5 * n + i] = state_.rho[i];
        }
    }
    return x;
}

void FiniteTempRelaxer::unpack_state(const std::vector<double>& x)
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    for (size_t i = 0; i < n; ++i)
    {
        state_.sigma[i] = x[i];
        state_.eta[i] = x[n + i];
        state_.alpha[i] = std::max(1e-12, x[2 * n + i]);
        state_.Xs[i] = std::max(1e-12, x[3 * n + i]);
        state_.omega[i] = x[4 * n + i];
        if (p_.increment_omega)
        {
            state_.rho[i] = x[5 * n + i];
        }
    }
}

FiniteTempRelaxer::Residuals FiniteTempRelaxer::compute_residuals() const
{
    const int n = p_.n_gridpoints;
    const double drbar = (p_.rbar_max - p_.rbar_min) / static_cast<double>(n - 1);

    Residuals R;
    R.sigma.assign(static_cast<size_t>(n), 0.0);
    R.eta.assign(static_cast<size_t>(n), 0.0);
    R.alpha.assign(static_cast<size_t>(n), 0.0);
    R.Xs.assign(static_cast<size_t>(n), 0.0);
    R.omega.assign(static_cast<size_t>(n), 0.0);
    R.rho.assign(static_cast<size_t>(n), 0.0);

    R.sigma[0] = state_.sigma[1] - state_.sigma[0];
    R.eta[0] = state_.eta[0];
    R.alpha[0] = state_.alpha[1] - state_.alpha[0];
    R.Xs[0] = state_.Xs[0] - 1.0;
    if (p_.increment_omega)
    {
        R.rho[0] = state_.rho[1] - state_.rho[0];
    }

    R.omega[0] = state_.omega[0] - p_.mu_chem_infty;

    const int last = n - 1;
    R.sigma[static_cast<size_t>(last)] = state_.sigma[static_cast<size_t>(last)] - 0.0;
    R.eta[static_cast<size_t>(last)] = state_.eta[static_cast<size_t>(last)];
    R.alpha[static_cast<size_t>(last)] = state_.alpha[static_cast<size_t>(last)] - 1.0;
    R.Xs[static_cast<size_t>(last)] = state_.Xs[static_cast<size_t>(last)] - 1.0;
    if (p_.increment_omega)
    {
        R.rho[static_cast<size_t>(last)] = particle_number_from_rho() - target_particle_number_;
    }

    if (!p_.increment_omega)
    {
        R.omega[static_cast<size_t>(last)] = state_.omega[static_cast<size_t>(last)] - p_.mu_chem_infty;
    }
    else
    {
        R.omega[static_cast<size_t>(last)] = state_.omega[static_cast<size_t>(last)] - state_.omega[static_cast<size_t>(last - 1)];
    }

    for (int i = 1; i < last; ++i)
    {
        const size_t j = static_cast<size_t>(i);
        const double rbar = state_.rbar[j];
        const double r = std::max(1e-12, areal_radius(rbar));
        const double ds_drb = d1_rbar(state_.sigma, i, drbar);
        const double de_drb = d1_rbar(state_.eta, i, drbar);
        const double da_drb = d1_rbar(state_.alpha, i, drbar);
        const double dx_drb = d1_rbar(state_.Xs, i, drbar);
        const double ds_dr = d1_r_from_rbar(ds_drb, rbar);
        const double de_dr = d1_r_from_rbar(de_drb, rbar);
        const double da_dr = d1_r_from_rbar(da_drb, rbar);
        const double dx_dr = d1_r_from_rbar(dx_drb, rbar);

        const double alpha = std::max(1e-12, state_.alpha[j]);
        const double Xs = std::max(1e-12, state_.Xs[j]);
        const double sigma = state_.sigma[j];
        const double sigma_sq = sigma * sigma;
        const double eta = state_.eta[j];
        const double X = eta * eta;
        const double temp_local = local_temperature(alpha);

        const double omega = state_.omega[j];
        const double mu_chem = omega / alpha;
        const double Xi = -mu_chem * mu_chem;
        const double P = pressure_Pthermal(temp_local, mu_chem, sigma_sq, X, Xi);
        const ThermoConjugates tc = conjugates_from_pressure(temp_local, mu_chem, sigma_sq, X, Xi);

        const double F = tc.F;
        const double source_X = mu_chem * mu_chem * tc.Upsilon + mu_chem * tc.n + tc.s * temp_local - P;

        R.sigma[j] = ds_dr - Xs * eta;
        R.eta[j] = de_dr + 2.0 * eta / r + (da_dr / alpha) * eta + Xs * sigma * F;
        R.alpha[j] = da_dr / alpha - ((Xs * Xs - 1.0) / (2.0 * r) + 4.0 * kPi * r * Xs * Xs * (eta * eta + P));
        R.Xs[j] = dx_dr / Xs - ((1.0 - Xs * Xs) / (2.0 * r) + 4.0 * kPi * r * Xs * Xs * source_X);
        if (p_.increment_omega)
        {
            const double drho_drb = d1_rbar(state_.rho, i, drbar);
            const double d2rho_drb2 = d2_rbar(state_.rho, i, drbar);
            const double drho_dr = d1_r_from_rbar(drho_drb, rbar);
            const double d2rho_dr2 = d2_r_from_rbar(drho_drb, d2rho_drb2, rbar);
            const double omega_safe = std::max(1e-12, std::abs(omega));
            const double source_rho = tc.Upsilon / (alpha * alpha) + tc.n / (alpha * omega_safe);
            R.rho[j] = d2rho_dr2 + (2.0 / r + (da_dr / alpha) - (dx_dr / Xs)) * drho_dr - source_rho;
        }

        if (!p_.increment_omega)
        {
            R.omega[j] = state_.omega[j] - p_.mu_chem_infty;
        }
        else
        {
            R.omega[j] = state_.omega[j] - state_.omega[j - 1];
        }
    }

    return R;
}

std::vector<double> FiniteTempRelaxer::pack_residual_vector(const Residuals& R) const
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    const size_t fields = active_field_count();
    std::vector<double> v(fields * n);
    for (size_t i = 0; i < n; ++i)
    {
        v[i] = R.sigma[i];
        v[n + i] = R.eta[i];
        v[2 * n + i] = R.alpha[i];
        v[3 * n + i] = R.Xs[i];
        v[4 * n + i] = R.omega[i];
        if (p_.increment_omega)
        {
            v[5 * n + i] = R.rho[i];
        }
    }
    return v;
}

double FiniteTempRelaxer::residual_norm(const Residuals& R) const
{
    return vector_norm_inf(pack_residual_vector(R));
}

double FiniteTempRelaxer::vector_norm_inf(const std::vector<double>& v) const
{
    double nrm = 0.0;
    for (double value : v)
    {
        nrm = std::max(nrm, std::abs(value));
    }
    return nrm;
}

void FiniteTempRelaxer::build_jacobian_fd(const std::vector<double>& x_base,
                                          const std::vector<double>& r_base,
                                          std::vector<double>& J) const
{
    const size_t nvar = x_base.size();
    const double eps = std::max(1e-12, p_.jacobian_fd_step);
    J.assign(nvar * nvar, 0.0);

    for (size_t col = 0; col < nvar; ++col)
    {
        std::vector<double> x_pert = x_base;
        const double step = eps * std::max(1.0, std::abs(x_base[col]));
        x_pert[col] += step;
        const_cast<FiniteTempRelaxer*>(this)->unpack_state(x_pert);
        const Residuals Rp = compute_residuals();
        const std::vector<double> r_pert = pack_residual_vector(Rp);

        for (size_t row = 0; row < nvar; ++row)
        {
            J[row * nvar + col] = (r_pert[row] - r_base[row]) / step;
        }
    }

    const_cast<FiniteTempRelaxer*>(this)->unpack_state(x_base);
}

bool FiniteTempRelaxer::solve_linear_system(std::vector<double> A,
                                            std::vector<double> b,
                                            std::vector<double>& x) const
{
    const size_t n = b.size();
    x.assign(n, 0.0);

    for (size_t k = 0; k < n; ++k)
    {
        size_t pivot = k;
        double pivot_value = std::abs(A[k * n + k]);
        for (size_t i = k + 1; i < n; ++i)
        {
            const double candidate = std::abs(A[i * n + k]);
            if (candidate > pivot_value)
            {
                pivot_value = candidate;
                pivot = i;
            }
        }

        if (pivot_value < 1e-18 || !std::isfinite(pivot_value))
        {
            if (p_.increment_omega)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[newton][charge] linear solve failed at pivot=" << k
                          << " pivot_value=" << pivot_value << '\n';
            }
            return false;
        }

        if (pivot != k)
        {
            for (size_t j = k; j < n; ++j)
            {
                std::swap(A[k * n + j], A[pivot * n + j]);
            }
            std::swap(b[k], b[pivot]);
        }

        const double diag = A[k * n + k];
        for (size_t i = k + 1; i < n; ++i)
        {
            const double factor = A[i * n + k] / diag;
            A[i * n + k] = 0.0;
            for (size_t j = k + 1; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = b[static_cast<size_t>(i)];
        for (size_t j = static_cast<size_t>(i) + 1; j < n; ++j)
        {
            sum -= A[static_cast<size_t>(i) * n + j] * x[j];
        }
        const double diag = A[static_cast<size_t>(i) * n + static_cast<size_t>(i)];
        if (std::abs(diag) < 1e-18 || !std::isfinite(diag))
        {
            if (p_.increment_omega)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[newton][charge] back-substitution failed at row=" << i
                          << " diag=" << diag << '\n';
            }
            return false;
        }
        x[static_cast<size_t>(i)] = sum / diag;
    }

    return true;
}

bool FiniteTempRelaxer::newton_solve_stage()
{
    const size_t nvar = active_field_count() * static_cast<size_t>(p_.n_gridpoints);
    std::vector<double> x = pack_state();
    std::vector<double> best_x = x;
    double best_norm = std::numeric_limits<double>::infinity();

    auto max_abs_with_index = [](const std::vector<double>& values) -> std::pair<double, int>
    {
        double max_value = 0.0;
        int max_index = 0;
        for (int i = 0; i < static_cast<int>(values.size()); ++i)
        {
            const double abs_value = std::abs(values[static_cast<size_t>(i)]);
            if (abs_value > max_value)
            {
                max_value = abs_value;
                max_index = i;
            }
        }
        return {max_value, max_index};
    };

    for (int iter = 0; iter < p_.max_newton_iterations; ++iter)
    {
        const Residuals R = compute_residuals();
        const std::vector<double> r_base = pack_residual_vector(R);
        const double base_norm = vector_norm_inf(r_base);
        if (std::isfinite(base_norm) && base_norm < best_norm)
        {
            best_norm = base_norm;
            best_x = x;
        }
        last_residual_norm_ = base_norm;
        last_newton_iterations_ = iter + 1;

        const auto [sig_max, sig_idx] = max_abs_with_index(R.sigma);
        const auto [eta_max, eta_idx] = max_abs_with_index(R.eta);
        const auto [alp_max, alp_idx] = max_abs_with_index(R.alpha);
        const auto [xs_max, xs_idx] = max_abs_with_index(R.Xs);
        const auto [om_max, om_idx] = max_abs_with_index(R.omega);

        if (p_.increment_omega)
        {
            const auto [rho_max, rho_idx] = max_abs_with_index(R.rho);
            std::cout << std::scientific << std::setprecision(6)
                      << "[newton][charge] iter=" << iter
                      << " |res|_inf=" << base_norm
                      << " sigma=" << sig_max << "@" << sig_idx
                      << " eta=" << eta_max << "@" << eta_idx
                      << " alpha=" << alp_max << "@" << alp_idx
                      << " Xs=" << xs_max << "@" << xs_idx
                      << " omega=" << om_max << "@" << om_idx
                      << " rho=" << rho_max << "@" << rho_idx
                      << " |bdy sigma0=" << R.sigma[0]
                      << " sigmaN=" << R.sigma[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << " eta0=" << R.eta[0]
                      << " etaN=" << R.eta[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << " alphaN=" << R.alpha[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << " XsN=" << R.Xs[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << " omegaN=" << R.omega[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << " rhoN=" << R.rho[static_cast<size_t>(p_.n_gridpoints - 1)]
                      << '\n';
        }

        if (base_norm < p_.newton_tolerance)
        {
            return true;
        }

        std::vector<double> J;
        build_jacobian_fd(x, r_base, J);

        std::vector<double> rhs = r_base;
        for (double& value : rhs)
        {
            value = -value;
        }

        std::vector<double> delta;
        if (!solve_linear_system(J, rhs, delta))
        {
            unpack_state(best_x);
            last_residual_norm_ = best_norm;
            return false;
        }

        double step = 1.0;
        bool accepted = false;
        double last_trial_norm = std::numeric_limits<double>::infinity();
        while (step >= p_.line_search_min_step)
        {
            std::vector<double> trial = x;
            for (size_t i = 0; i < nvar; ++i)
            {
                trial[i] += step * delta[i];
            }

            unpack_state(trial);
            const double trial_norm = residual_norm(compute_residuals());
            last_trial_norm = trial_norm;
            if (std::isfinite(trial_norm) && trial_norm < base_norm)
            {
                x = trial;
                accepted = true;
                break;
            }

            step *= p_.line_search_shrink;
        }

        if (!accepted)
        {
            if (p_.increment_omega)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[newton][charge] line search failed: base=" << base_norm
                          << " last_trial=" << last_trial_norm
                          << " min_step=" << p_.line_search_min_step << '\n';
            }
            unpack_state(best_x);
            last_residual_norm_ = best_norm;
            return false;
        }

        unpack_state(x);
    }

    const double final_norm = residual_norm(compute_residuals());
    if (std::isfinite(final_norm) && final_norm < best_norm)
    {
        best_norm = final_norm;
        best_x = x;
    }
    if (best_norm < p_.newton_tolerance)
    {
        unpack_state(best_x);
        last_residual_norm_ = best_norm;
        return true;
    }

    unpack_state(best_x);
    last_residual_norm_ = best_norm;
    return false;
}

bool FiniteTempRelaxer::relax_current_state()
{
    enforce_boundary_values();
    return newton_solve_stage();
}

size_t FiniteTempRelaxer::active_field_count() const
{
    return p_.increment_omega ? 6U : 5U;
}

std::vector<double> FiniteTempRelaxer::stage_targets(bool charge_mode, double start_value) const
{
    std::vector<double> stages;
    double step = charge_mode ? p_.charge_step : p_.omega_step;
    if (std::abs(step) < 1e-12)
    {
        step = 1e-12;
    }
    const int count = std::max(1, p_.num_models);
    stages.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i)
    {
        if (charge_mode)
        {
            stages.push_back(start_value - step * static_cast<double>(i + 1));
        }
        else
        {
            stages.push_back(start_value - step * static_cast<double>(i));
        }
    }

    if (stages.empty())
    {
        stages.push_back(p_.mu_chem_infty);
    }

    return stages;
}

bool FiniteTempRelaxer::run(bool quiet)
{
    build_grid();
    initialize_seed();

    // initialize sigma trackers used by shooting fallback
    last_saved_sigma0_ = state_.sigma.front();
    last_sigma_step_ = 0.0;

    {
        std::ofstream models_out{"bs_models.dat"};
        if (models_out)
        {
            models_out << "# sigma0 M omega r99 N\n";
        }
    }

    auto max_abs_with_index = [](const std::vector<double>& values) -> std::pair<double, int>
    {
        double max_value = 0.0;
        int max_index = 0;
        for (int i = 0; i < static_cast<int>(values.size()); ++i)
        {
            const double abs_value = std::abs(values[static_cast<size_t>(i)]);
            if (abs_value > max_value)
            {
                max_value = abs_value;
                max_index = i;
            }
        }
        return {max_value, max_index};
    };

    auto print_residual_breakdown = [&](const std::string& label)
    {
        const Residuals R = compute_residuals();
        const auto [sig_max, sig_idx] = max_abs_with_index(R.sigma);
        const auto [eta_max, eta_idx] = max_abs_with_index(R.eta);
        const auto [alp_max, alp_idx] = max_abs_with_index(R.alpha);
        const auto [xs_max, xs_idx] = max_abs_with_index(R.Xs);
        const auto [om_max, om_idx] = max_abs_with_index(R.omega);
        const std::vector<double> packed = pack_residual_vector(R);
        const double inf_norm = vector_norm_inf(packed);

        std::cout << std::scientific << std::setprecision(6)
                  << "[residuals][" << label << "]"
                  << " inf=" << inf_norm
                  << " sigma=" << sig_max << "@" << sig_idx
                  << " eta=" << eta_max << "@" << eta_idx
                  << " alpha=" << alp_max << "@" << alp_idx
                  << " Xs=" << xs_max << "@" << xs_idx
                  << " omega=" << om_max << "@" << om_idx;

        if (p_.increment_omega)
        {
            const auto [rho_max, rho_idx] = max_abs_with_index(R.rho);
            std::cout << " rho=" << rho_max << "@" << rho_idx;
        }

        const int n_last = p_.n_gridpoints - 1;
        std::cout << " |bdy sigma0=" << R.sigma[0]
                  << " sigmaN=" << R.sigma[static_cast<size_t>(n_last)]
                  << " eta0=" << R.eta[0]
                  << " etaN=" << R.eta[static_cast<size_t>(n_last)]
                  << " alphaN=" << R.alpha[static_cast<size_t>(n_last)]
                  << " XsN=" << R.Xs[static_cast<size_t>(n_last)]
                  << " omegaN=" << R.omega[static_cast<size_t>(n_last)];

        if (p_.increment_omega)
        {
            std::cout << " rhoN=" << R.rho[static_cast<size_t>(n_last)];
        }

        std::cout << '\n';
    };

    // Temperature evolution stage: if temperature_target > 0, evolve from T=0 to T_target
    std::cout << "[DEBUG] temperature_target = " << p_.temperature_target
              << ", temperature_step = " << p_.temperature_step << '\n';
    
    if (p_.temperature_target > 0.0)
    {
        std::cout << "[INFO] === BEGINNING TEMPERATURE EVOLUTION ===" << '\n';
        std::cout << "[INFO] Evolving from T=0 to T=" << p_.temperature_target
                  << " in steps of " << p_.temperature_step << '\n';

        double current_T = 0.0;
        int temp_step_count = 0;

        while (current_T < p_.temperature_target)
        {
            // Increment temperature
            current_T = std::min(p_.temperature_target, current_T + p_.temperature_step);
            current_temperature_ = current_T;
            ++temp_step_count;

            std::cout << "[TEMP] Step " << temp_step_count << ": T = " << current_T << " ... ";
            std::cout.flush();

            enforce_boundary_values();
            const bool ok = newton_solve_stage();
            // print_residual_breakdown("temp_step_post_newton");
            if (!ok)
            {
                std::cout << "FAILED (residual = " << last_residual_norm_ << ")\n";
                if (p_.force_continue)
                {
                    std::cout << "[WARNING] Temperature step did not reach tolerance at T = " << current_T
                              << ", using best attained residual = " << last_residual_norm_
                              << " and continuing.\n";
                    print_residual_breakdown("temp_step_forced_continue");

                    // Keep a snapshot of the best state retained by Newton fallback.
                    write_profile(p_.output_profile);
                    continue;
                }

                std::cout << "[ERROR] Temperature evolution failed at T = " << current_T << '\n';
                print_residual_breakdown("temp_step_failure");
                
                // Print detailed diagnostics at the point of failure
                print_field_diagnostics("at_failure");
                
                // Find and print the grid point with max eta residual
                const Residuals R_fail = compute_residuals();
                double max_eta_res = 0.0;
                int max_eta_idx = 0;
                for (int ii = 0; ii < p_.n_gridpoints; ++ii)
                {
                    if (std::abs(R_fail.eta[static_cast<size_t>(ii)]) > max_eta_res)
                    {
                        max_eta_res = std::abs(R_fail.eta[static_cast<size_t>(ii)]);
                        max_eta_idx = ii;
                    }
                }
                std::cout << "[ERROR] Max eta residual at index " << max_eta_idx << " with value " << max_eta_res << "\n";
                print_residual_components(max_eta_idx);
                
                write_profile(p_.output_profile);
                return false;
            }

            std::cout << "OK (" << last_newton_iterations_ << " iterations, residual = "
                      << last_residual_norm_ << ")\n";

            // Sanity checks on constructed model
            if (p_.mu_chem_infty > 1.0)
            {
                std::cout << "[ERROR] Sanity check failed: omega = " << p_.mu_chem_infty 
                          << " > 1.0 at T = " << current_T << "\n";
                return false;
            }
            
            const double sigma_min = *std::min_element(state_.sigma.begin(), state_.sigma.end());
            if (sigma_min < -1e-6)  // Allow small numerical errors
            {
                std::cout << "[WARNING] Sanity check failed: sigma_min = " << sigma_min 
                          << " < -1e-6 at T = " << current_T << "\n";
                // return false;
            }

            // Update sigma trackers for possible fallback stepping
            const double new_sigma0 = state_.sigma.front();
            last_sigma_step_ = new_sigma0 - last_saved_sigma0_;
            last_saved_sigma0_ = new_sigma0;
            // Write the current successful profile so progress can be monitored
            write_profile(p_.output_profile);
        }

        std::cout << "[INFO] === TEMPERATURE EVOLUTION COMPLETE ===" << '\n';
        std::cout << "[INFO] Total steps: " << temp_step_count << ", final T = " << current_T << '\n';
        
        if (p_.first_thermalized_only)
        {
            std::cout << "[INFO] first_thermalized_only is true; writing profile and exiting.\n";
            write_profile("profile.dat");
            std::cout << "[INFO] Profile written to profile.dat.\n";
            return true;
        }
        
        std::cout << "[INFO] Starting omega continuation at T = " << p_.temperature_target << '\n';
    }
    else
    {
        std::cout << "[INFO] No temperature evolution (temperature_target = 0)\n";
    }

    if (p_.first_thermalized_only && p_.temperature_target <= 0.0)
    {
        std::cout << "[INFO] first_thermalized_only is true and temperature_target is 0; writing zero-temperature profile and exiting.\n";
        write_profile("profile.dat");
        std::cout << "[INFO] Profile written to profile.dat.\n";
        return true;
    }

    std::cout << "[INFO] === STARTING OMEGA CONTINUATION ===" << '\n';

    // Clear threshold tracking sets for this continuation phase
    crossed_hi_sigma_thresholds_.clear();
    crossed_lo_sigma_thresholds_.clear();

    auto run_continuation = [&](bool charge_mode, double start_value, bool preserve_initial_freq = false) -> bool
    {
        const std::vector<double> stages = stage_targets(charge_mode, start_value);
        for (size_t stage = 0; stage < stages.size(); ++stage)
        {
            // Only reset temperature to 0 if we didn't do temperature evolution.
            // If temperature_target > 0, keep current_temperature_ at the target value.
            if (p_.temperature_target <= 0.0)
            {
                current_temperature_ = 0.0;
            }
            if (charge_mode)
            {
                target_particle_number_ = stages[stage];
            }
            else
            {
                // If preserve_initial_freq is set, keep the current p_.mu_chem_infty
                // for the first stage so that a freshly-shot model's frequency
                // is not overwritten by the continuation stage value.
                if (!(preserve_initial_freq && stage == 0))
                {
                    p_.mu_chem_infty = stages[stage];
                }
            }
            if (!quiet)
            {
                std::cout << "[info] Stage " << (stage + 1) << "/" << stages.size()
                          << (charge_mode ? ": N_target = " : ": omega_infty = ")
                          << (charge_mode ? target_particle_number_ : p_.mu_chem_infty) << '\n';
            }

            enforce_boundary_values();
            const bool ok = newton_solve_stage();
            if (!ok)
            {
                if (!quiet)
                {
                    std::cout << "[warn] Newton stage failed at "
                              << (charge_mode ? "N_target = " : "omega_infty = ")
                              << (charge_mode ? target_particle_number_ : p_.mu_chem_infty)
                              << ", residual = " << last_residual_norm_ << '\n';
                }
                write_profile(p_.output_profile);
                return false;
            }

            // Sanity checks on constructed model
            if (p_.mu_chem_infty > 1.0)
            {
                if (!quiet)
                {
                    std::cout << "[ERROR] Sanity check failed: omega = " << p_.mu_chem_infty 
                              << " > 1.0\n";
                }
                return false;
            }

            const double sigma_min = *std::min_element(state_.sigma.begin(), state_.sigma.end());
            if (sigma_min < -1e-6)  // Allow small numerical errors
            {
                if (!quiet)
                {
                    std::cout << "[WARNING] Sanity check failed: sigma_min = " << sigma_min 
                              << " < -1e-6\n";
                }
                //return false;
            }

            // Check for sigma threshold crossings and refine grid if needed
            const double sigma0 = state_.sigma.front();
            for (const double hi_sig : p_.hi_sigma_vals)
            {
                if (sigma0 > hi_sig && crossed_hi_sigma_thresholds_.find(hi_sig) == crossed_hi_sigma_thresholds_.end())
                {
                    crossed_hi_sigma_thresholds_.insert(hi_sig);
                    if (!quiet)
                    {
                        std::cout << "[info] Central amplitude sigma0 = " << sigma0 
                                  << " exceeded high threshold " << hi_sig 
                                  << "; refining grid by factor " << p_.resolution_factor << '\n';
                    }
                    refine_grid(p_.resolution_factor);
                    if (!quiet)
                    {
                        std::cout << "[info] Grid refinement complete; new n_gridpoints = " 
                                  << p_.n_gridpoints << '\n';
                    }
                    break;  // Only refine once per stage
                }
            }

            for (const double lo_sig : p_.lo_sigma_vals)
            {
                if (sigma0 < lo_sig && crossed_lo_sigma_thresholds_.find(lo_sig) == crossed_lo_sigma_thresholds_.end())
                {
                    crossed_lo_sigma_thresholds_.insert(lo_sig);
                    if (!quiet)
                    {
                        std::cout << "[info] Central amplitude sigma0 = " << sigma0 
                                  << " went below low threshold " << lo_sig 
                                  << "; refining grid by factor " << p_.resolution_factor << '\n';
                    }
                    refine_grid(p_.resolution_factor);
                    if (!quiet)
                    {
                        std::cout << "[info] Grid refinement complete; new n_gridpoints = " 
                                  << p_.n_gridpoints << '\n';
                    }
                    break;  // Only refine once per stage
                }
            }

            // Update sigma trackers for possible fallback stepping
            {
                const double new_sigma0 = state_.sigma.front();
                last_sigma_step_ = new_sigma0 - last_saved_sigma0_;
                last_saved_sigma0_ = new_sigma0;
                (void)0; // no-op: do not track particle number for charge fallback
            }

            if (!quiet)
            {
                std::cout << "[info] Converged in " << last_newton_iterations_
                          << " Newton iterations, residual = " << last_residual_norm_ << '\n';
            }

            append_model_summary("bs_models.dat");
            if (p_.save_model_profiles)
            {
                // one profile per family member, e.g. for picking SBSE backgrounds along a branch
                static int model_counter = 0;
                std::ostringstream name;
                name << "model_" << std::setw(4) << std::setfill('0') << model_counter++ << ".dat";
                write_profile(name.str());
            }
        }
        return true;
    };

    if (!p_.increment_omega)
    {
        const bool omega_ok = run_continuation(false, p_.mu_chem_infty);
        if (!omega_ok)
        {
            if (!quiet)
            {
                std::cout << "[info] Omega continuation failed; attempting shooting fallback by stepping central amplitude\n";
            }

            // compute proposed new central amplitude using previous step
            double step = last_sigma_step_;
            if (std::abs(step) < 1e-14)
            {
                step = 1e-3 * std::max(1.0, std::abs(last_saved_sigma0_));
            }
            const double new_sigma0 = last_saved_sigma0_ + p_.sigma_jump_on_fail * step;

            if (!quiet)
            {
                std::cout << "[info] Shooting for sigma0 = " << new_sigma0 << '\n';
            }

            const bool shot_ok = shoot_for_sigma(new_sigma0);
            if (shot_ok)
            {
                if (!quiet)
                {
                    std::cout << "[info] Shooting produced a seed\n";
                }

                // If we did temperature evolution initially, re-thermalize the new seed
                if (p_.temperature_target > 0.0)
                {
                    if (!quiet)
                    {
                        std::cout << "[info] Re-thermalizing new seed from T=0 to T=" << p_.temperature_target << "\n";
                    }

                    double current_T = 0.0;
                    int temp_step_count = 0;
                    bool rethermal_ok = true;

                    while (current_T < p_.temperature_target)
                    {
                        current_T = std::min(p_.temperature_target, current_T + p_.temperature_step);
                        current_temperature_ = current_T;
                        ++temp_step_count;

                        enforce_boundary_values();
                        const bool temp_ok = newton_solve_stage();
                        if (!temp_ok)
                        {
                            if (!quiet)
                            {
                                std::cout << "[ERROR] Re-thermalization failed at T = " << current_T << '\n';
                            }
                            rethermal_ok = false;
                            break;
                        }

                        // Sanity checks on constructed model
                        if (p_.mu_chem_infty > 1.0)
                        {
                            if (!quiet)
                            {
                                std::cout << "[ERROR] Sanity check failed during re-thermalization: omega = " 
                                          << p_.mu_chem_infty << " > 1.0 at T = " << current_T << "\n";
                            }
                            rethermal_ok = false;
                            break;
                        }
                        
                        const double sigma_min = *std::min_element(state_.sigma.begin(), state_.sigma.end());
                        if (sigma_min < -1e-6)
                        {
                            if (!quiet)
                            {
                                std::cout << "[WARNING] Sanity check failed during re-thermalization: sigma_min = " 
                                          << sigma_min << " < -1e-6 at T = " << current_T << "\n";
                            }
                            //rethermal_ok = false;
                            //break;
                        }

                        // Update sigma trackers
                        const double new_sigma0_temp = state_.sigma.front();
                        last_sigma_step_ = new_sigma0_temp - last_saved_sigma0_;
                        last_saved_sigma0_ = new_sigma0_temp;
                    }

                    if (!rethermal_ok)
                    {
                        if (!quiet)
                        {
                            std::cout << "[warn] Re-thermalization did not succeed; aborting (no charge fallback)\n";
                        }
                        return false;
                    }

                    if (!quiet)
                    {
                        std::cout << "[info] Re-thermalization complete; reversing omega step direction and retrying omega continuation\n";
                    }
                }
                else
                {
                    if (!quiet)
                    {
                        std::cout << "[info] Reversing omega step direction and retrying omega continuation\n";
                    }
                }

                // Reverse omega stepping direction because we passed a turning point
                p_.omega_step = -p_.omega_step;

                // reset trackers to reflect the new shot seed (or re-thermalized seed)
                last_sigma_step_ = state_.sigma.front() - last_saved_sigma0_;
                last_saved_sigma0_ = state_.sigma.front();

                // Try omega continuation again from current omega_infty
                // Preserve the shooter's frequency for the first stage so it isn't overwritten.
                const bool retry_ok = run_continuation(false, p_.mu_chem_infty, true);
                if (!retry_ok)
                {
                    if (!quiet)
                    {
                        std::cout << "[warn] Shooting fallback did not recover omega continuation; aborting (no charge fallback)\n";
                    }
                    write_profile(p_.output_profile);
                    return false;
                }
            }
            else
            {
                if (!quiet)
                {
                    std::cout << "[info] Shooting fallback failed; aborting (no charge fallback)\n";
                }
                write_profile(p_.output_profile);
                return false;
            }
        }
    }
    else
    {
        if (!run_continuation(true, p_.sigma_central))
        {
            write_profile(p_.output_profile);
            return false;
        }
    }

    write_profile(p_.output_profile);
    return true;
}

void FiniteTempRelaxer::write_profile(const std::string& filename) const
{
    std::ofstream out{filename};
    if (!out)
    {
        std::cerr << "Could not open " << filename << " for writing\n";
        return;
    }

    out << std::scientific << std::setprecision(12);
    // Metadata lines (read back by load_profile; older readers skip '#' lines).
    out << "# T_inf = " << current_temperature_ << '\n';
    out << "# omega_inf = " << omega_infinity() << '\n';
    out << "# mu = " << p_.mu << " lambda = " << p_.lambda << " solitonic = " << p_.solitonic
        << " sigma_scale = " << p_.sigma_scale << " G = " << p_.G << '\n';
    out << "# n_gridpoints = " << p_.n_gridpoints << " rbar_min = " << p_.rbar_min << " rbar_max = " << p_.rbar_max << '\n';
    // Columns 1-11 are unchanged from the original relaxer output; 12-17 were appended for SBSE.
    out << "# rbar r sigma eta alpha Xs omega Tlocal P_scalar P_thermal_suppressed P_thermal_raw"
        << " mu_local P_total s n Upsilon rho\n";

    for (int i = 0; i < p_.n_gridpoints - 1; ++i)
    {
        const double alpha = std::max(1e-12, state_.alpha[static_cast<size_t>(i)]);
        const double Tlocal = local_temperature(alpha);
        const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
        const double sigma_sq = state_.sigma[static_cast<size_t>(i)] * state_.sigma[static_cast<size_t>(i)];
        const double X = state_.eta[static_cast<size_t>(i)] * state_.eta[static_cast<size_t>(i)];
        const double Xi = -std::pow(state_.omega[static_cast<size_t>(i)] / alpha, 2.0);
        const double mu_chem = state_.omega[static_cast<size_t>(i)] / alpha;
        const double P_scalar = scalar_pressure_term(sigma_sq, X, Xi);
        const double P_thermal = thermal_pressure_term(Tlocal, mu_chem, sigma_sq, X, Xi);
        const double suppression_factor = compute_suppression_factor(sigma_sq);
        
        const double P_thermal_suppressed = suppression_factor * P_thermal;

        out << state_.rbar[static_cast<size_t>(i)] << ' '
            << r << ' '
            << state_.sigma[static_cast<size_t>(i)] << ' '
            << state_.eta[static_cast<size_t>(i)] << ' '
            << state_.alpha[static_cast<size_t>(i)] << ' '
            << state_.Xs[static_cast<size_t>(i)] << ' '
            << state_.omega[static_cast<size_t>(i)] << ' '
            << Tlocal << ' '
            << P_scalar << ' '
            << P_thermal_suppressed << ' '
            << P_thermal;

        const PointThermo th = thermo_at(state_.sigma[static_cast<size_t>(i)], state_.eta[static_cast<size_t>(i)],
                                         alpha, state_.omega[static_cast<size_t>(i)]);
        out << ' ' << th.mu << ' ' << th.P << ' ' << th.s << ' ' << th.n << ' ' << th.Upsilon << ' ' << th.rho << '\n';
    }
}

FiniteTempRelaxer::PointThermo FiniteTempRelaxer::thermo_at(double sigma, double eta, double alpha, double omega) const
{
    PointThermo th;
    const double a = std::max(1e-12, alpha);
    th.T = local_temperature(a);
    th.mu = omega / a;
    th.X = eta * eta;
    th.Xi = -th.mu * th.mu;
    const double sigma_sq = sigma * sigma;

    th.P = pressure_Pthermal(th.T, th.mu, sigma_sq, th.X, th.Xi);
    const ThermoConjugates tc = conjugates_from_pressure(th.T, th.mu, sigma_sq, th.X, th.Xi);
    th.s = tc.s;
    th.n = tc.n;
    th.F = tc.F;
    th.Upsilon = tc.Upsilon;
    th.w = th.T * th.s + th.mu * th.n;
    // Same combination as the source of the Xs equation in compute_residuals().
    th.rho = th.mu * th.mu * th.Upsilon + th.mu * th.n + th.s * th.T - th.P;
    return th;
}

bool FiniteTempRelaxer::load_profile(const std::string& filename, bool quiet)
{
    std::ifstream in{filename};
    if (!in)
    {
        std::cerr << "FiniteTempRelaxer::load_profile: could not open " << filename << '\n';
        return false;
    }

    std::vector<double> rbar, sigma, eta, alpha, Xs, omega, Tloc;
    double T_inf_meta = -1.0;
    double rbar_max_meta = -1.0;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (line[0] == '#')
        {
            std::string meta = line.substr(1);
            fill_parameter(meta, "T_inf = ", T_inf_meta, true);
            const size_t pos = line.find("rbar_max = ");
            if (pos != std::string::npos)
            {
                std::stringstream ss(line.substr(pos + 11));
                ss >> rbar_max_meta;
            }
            continue;
        }
        std::stringstream ss(line);
        double v[8];
        bool ok = true;
        for (double& x : v)
        {
            if (!(ss >> x))
            {
                ok = false;
                break;
            }
        }
        if (!ok)
        {
            continue;
        }
        rbar.push_back(v[0]);
        sigma.push_back(v[2]);
        eta.push_back(v[3]);
        alpha.push_back(v[4]);
        Xs.push_back(v[5]);
        omega.push_back(v[6]);
        Tloc.push_back(v[7]);
    }

    if (rbar.size() < 3)
    {
        std::cerr << "FiniteTempRelaxer::load_profile: too few rows in " << filename << '\n';
        return false;
    }

    // write_profile() omits the last (rbar = rbar_max) point; the grid is uniform in rbar.
    const double drbar = (rbar.back() - rbar.front()) / static_cast<double>(rbar.size() - 1);
    for (size_t i = 1; i < rbar.size(); ++i)
    {
        if (std::abs((rbar[i] - rbar[i - 1]) - drbar) > 1e-6 * drbar)
        {
            std::cerr << "FiniteTempRelaxer::load_profile: rbar grid in " << filename << " is not uniform\n";
            return false;
        }
    }
    const double rbar_end = rbar.back() + drbar;
    if (rbar_max_meta > 0.0 && std::abs(rbar_end - rbar_max_meta) > 1e-6)
    {
        std::cerr << "FiniteTempRelaxer::load_profile: WARNING rbar_max from header (" << rbar_max_meta
                  << ") differs from the grid end (" << rbar_end << ")\n";
    }

    p_.n_gridpoints = static_cast<int>(rbar.size()) + 1;
    p_.rbar_min = rbar.front();
    p_.rbar_max = rbar_end;
    build_grid();

    const size_t n = rbar.size();
    for (size_t i = 0; i < n; ++i)
    {
        state_.sigma[i] = sigma[i];
        state_.eta[i] = eta[i];
        state_.alpha[i] = alpha[i];
        state_.Xs[i] = Xs[i];
        state_.omega[i] = omega[i];
        state_.rho[i] = 1.0;
    }
    p_.mu_chem_infty = omega.back();
    state_.omega[n] = p_.mu_chem_infty;
    state_.rho[n] = 1.0;
    {
        // boundary values at r = infinity (normally set by enforce_boundary_values)
        state_.sigma[n] = 0.0;
        state_.eta[n] = 0.0;
        state_.alpha[n] = 1.0;
        state_.Xs[n] = 1.0;
    }

    // Tolman: T_local * alpha = T_inf. Prefer the header value; otherwise recover it from the columns.
    double T_inf = T_inf_meta;
    if (T_inf < 0.0)
    {
        T_inf = Tloc.front() * alpha.front();
    }
    double max_dev = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        max_dev = std::max(max_dev, std::abs(Tloc[i] * alpha[i] - T_inf));
    }
    current_temperature_ = T_inf;

    if (!quiet)
    {
        std::cout << "[FiniteTempRelaxer] loaded " << filename << ": " << p_.n_gridpoints << " points, rbar in ["
                  << p_.rbar_min << ", " << p_.rbar_max << "], T_inf = " << T_inf
                  << ", omega_inf = " << p_.mu_chem_infty << ", sigma0 = " << state_.sigma.front()
                  << ", max |T_local alpha - T_inf| = " << max_dev << '\n';
    }
    return true;
}

void FiniteTempRelaxer::refine_grid(double res_factor)
{
    const int old_n = p_.n_gridpoints;
    const int new_n = std::max(3, static_cast<int>(std::round(old_n * res_factor)));

    if (new_n == old_n)
    {
        return;  // No refinement needed
    }

    // Save old state
    std::vector<double> old_rbar = state_.rbar;
    std::vector<double> old_sigma = state_.sigma;
    std::vector<double> old_eta = state_.eta;
    std::vector<double> old_alpha = state_.alpha;
    std::vector<double> old_Xs = state_.Xs;
    std::vector<double> old_omega = state_.omega;
    std::vector<double> old_rho = state_.rho;

    // Update n_gridpoints and build new grid
    p_.n_gridpoints = new_n;
    build_grid();  // This recreates state_.rbar with new size

    // Interpolate all fields onto the new grid
    for (int i = 0; i < new_n; ++i)
    {
        const double new_rbar_i = state_.rbar[static_cast<size_t>(i)];
        state_.sigma[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_sigma, new_rbar_i);
        state_.eta[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_eta, new_rbar_i);
        state_.alpha[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_alpha, new_rbar_i);
        state_.Xs[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_Xs, new_rbar_i);
        state_.omega[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_omega, new_rbar_i);
        state_.rho[static_cast<size_t>(i)] = interp_linear_sorted(old_rbar, old_rho, new_rbar_i);
    }

    // Ensure boundary conditions are satisfied
    enforce_boundary_values();
}

double FiniteTempRelaxer::metric_mass(double r, double Xs) const
{
    const double safe_Xs = std::max(1e-12, Xs);
    return 0.5 * r * (1.0 - 1.0 / (safe_Xs * safe_Xs));
}

double FiniteTempRelaxer::asymptotic_mass() const
{
    const int last_interior = std::max(0, p_.n_gridpoints - 2);
    const double r = areal_radius(state_.rbar[static_cast<size_t>(last_interior)]);
    return metric_mass(r, state_.Xs[static_cast<size_t>(last_interior)]);
}

double FiniteTempRelaxer::radius_at_mass_fraction(double mass_fraction) const
{
    const int last_interior = std::max(1, p_.n_gridpoints - 2);
    const double mass_inf = asymptotic_mass();
    const double mf = std::min(1.0, std::max(0.0, mass_fraction));
    const double target_mass = mf * mass_inf;

    if (!std::isfinite(mass_inf) || mass_inf <= 0.0)
    {
        return areal_radius(state_.rbar[static_cast<size_t>(last_interior)]);
    }

    // Build vectors of rbar and masses for interior points
    std::vector<double> rbar_vals;
    std::vector<double> masses;
    
    for (int i = 0; i <= last_interior; ++i)
    {
        const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
        const double m = metric_mass(r, state_.Xs[static_cast<size_t>(i)]);
        rbar_vals.push_back(state_.rbar[static_cast<size_t>(i)]);
        masses.push_back(m);
    }

    if (rbar_vals.size() < 2)
    {
        return areal_radius(rbar_vals.back());
    }

    // Handle out-of-bounds cases
    if (target_mass <= masses.front())
    {
        return areal_radius(rbar_vals.front());
    }

    if (target_mass >= masses.back())
    {
        return areal_radius(rbar_vals.back());
    }

    // Linear interpolation in (mass, rbar) space to find rbar_target
    for (size_t i = 0; i < masses.size() - 1; ++i)
    {
        const double m_prev = masses[i];
        const double m_curr = masses[i + 1];
        const bool bracketed = (m_prev <= target_mass && target_mass <= m_curr) ||
                               (m_curr <= target_mass && target_mass <= m_prev);
        if (bracketed)
        {
            const double denom = m_curr - m_prev;
            if (std::abs(denom) < 1e-18)
            {
                return areal_radius(rbar_vals[i + 1]);
            }
            const double t = (target_mass - m_prev) / denom;
            const double rbar_target = rbar_vals[i] * (1.0 - t) + rbar_vals[i + 1] * t;
            return areal_radius(rbar_target);
        }
    }

    return areal_radius(rbar_vals.back());
}

double FiniteTempRelaxer::particle_number() const
{
    const int n = p_.n_gridpoints;
    if (n < 2)
    {
        return 0.0;
    }

    const int last_interior = std::max(1, n - 2);
    auto integrand = [this](int i) -> double
    {
        const size_t j = static_cast<size_t>(i);
        const double r = areal_radius(state_.rbar[j]);
        const double alpha = std::max(1e-12, state_.alpha[j]);
        const double Xs = std::max(1e-12, state_.Xs[j]);
        const double omega = state_.omega[j];
        const double sigma_sq = state_.sigma[j] * state_.sigma[j];
        const double X = state_.eta[j] * state_.eta[j];
        const double mu_chem = omega / alpha;
        const double Xi = -mu_chem * mu_chem;
        const double Tlocal = local_temperature(alpha);
        const ThermoConjugates tc = conjugates_from_pressure(Tlocal, mu_chem, sigma_sq, X, Xi);

        return 4.0 * kPi * Xs * r * r * (tc.Upsilon * omega / alpha + tc.n);
    };

    double integral = 0.0;
    double prev_r = areal_radius(state_.rbar[0]);
    double prev_f = integrand(0);
    for (int i = 1; i <= last_interior; ++i)
    {
        const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
        const double f = integrand(i);
        integral += 0.5 * (prev_f + f) * (r - prev_r);
        prev_r = r;
        prev_f = f;
    }

    return integral;
}

double FiniteTempRelaxer::particle_number_from_rho() const
{
    const int n = p_.n_gridpoints;
    if (n < 2)
    {
        return 0.0;
    }

    const int i = std::max(1, n - 2);
    const double rbar = state_.rbar[static_cast<size_t>(i)];
    const double r = std::max(1e-12, areal_radius(rbar));
    const double drho_drb = d1_rbar(state_.rho, i, (p_.rbar_max - p_.rbar_min) / static_cast<double>(n - 1));
    const double drho_dr = d1_r_from_rbar(drho_drb, rbar);
    return 8.0 * kPi * r * r * drho_dr;
}

void FiniteTempRelaxer::append_model_summary(const std::string& filename) const
{
    std::ofstream out{filename, std::ios::app};
    if (!out)
    {
        std::cerr << "Could not open " << filename << " for appending\n";
        return;
    }

    const double sigma0 = state_.sigma.front();
    const double M = asymptotic_mass();
    const double omega = state_.omega.front();
    const double r99 = radius_at_mass_fraction(0.99);
    const double r90 = radius_at_mass_fraction(0.90);
    const double N = particle_number();

    out << std::scientific << std::setprecision(12)
        << sigma0 << ' '
        << M << ' '
        << omega << ' '
        << r99 << ' '
        << r90 << ' '
        << N << '\n';
}

double FiniteTempRelaxer::d1_rbar(const std::vector<double>& u, int i, double h) const
{
    const int N = static_cast<int>(u.size());
    if (i <= 0)
    {
        return 0.0;
    }
    if (i >= 2 && i <= N - 3)
    {
        return fivePointDeriv(h, 1,
                              u[i - 2], u[i - 1], u[i], u[i + 1], u[i + 2]);
    }
    if (i >= N - 1)
    {
        return (u[N - 1] - u[N - 2]) / h;
    }
    return (u[i + 1] - u[i - 1]) / (2.0 * h);
}

double FiniteTempRelaxer::d2_rbar(const std::vector<double>& u, int i, double h) const
{
    const int N = static_cast<int>(u.size());
    if (i <= 0)
    {
        return 2.0 * (u[1] - u[0]) / (h * h);
    }
    if (i >= 2 && i <= N - 3)
    {
        return fivePointDeriv(h, 2,
                              u[i - 2], u[i - 1], u[i], u[i + 1], u[i + 2]);
    }
    if (i >= N - 1)
    {
        return 0.0;
    }
    return (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (h * h);
}

double FiniteTempRelaxer::d1_r_from_rbar(double u_b, double rbar) const
{
    return compactification_->d1_r_factor(rbar) * u_b;
}

double FiniteTempRelaxer::d2_r_from_rbar(double u_b, double u_bb, double rbar) const
{
    return compactification_->d2_r_factor_first(rbar) * u_bb 
         + compactification_->d2_r_factor_second(rbar) * u_b;
}

void FiniteTempRelaxer::print_field_diagnostics(const std::string& label) const
{
    std::cout << "\n[DIAG] Field diagnostics at T = " << current_temperature_ << " (" << label << ")\n";
    std::cout << "[DIAG] Grid points: " << p_.n_gridpoints << "\n";
    std::cout << "[DIAG] Index |    rbar    |      r     |    sigma   |     eta    |    alpha   |     Xs     |    omega\n";
    std::cout << "[DIAG] ------|------------|------------|------------|------------|------------|------------|----------\n";
    
    for (int i = 0; i < p_.n_gridpoints; i += std::max(1, p_.n_gridpoints / 10))
    {
        const size_t j = static_cast<size_t>(i);
        const double r = areal_radius(state_.rbar[j]);
        std::cout << std::scientific << std::setprecision(4)
                  << "[DIAG] " << std::setw(5) << i << " | "
                  << state_.rbar[j] << " | "
                  << r << " | "
                  << state_.sigma[j] << " | "
                  << state_.eta[j] << " | "
                  << state_.alpha[j] << " | "
                  << state_.Xs[j] << " | "
                  << state_.omega[j] << "\n";
    }
    
    // Also print the last few points (near boundary)
    std::cout << "[DIAG] ------|------------|------------|------------|------------|------------|------------|----------\n";
    for (int i = std::max(1, p_.n_gridpoints - 5); i < p_.n_gridpoints; ++i)
    {
        const size_t j = static_cast<size_t>(i);
        const double r = areal_radius(state_.rbar[j]);
        std::cout << std::scientific << std::setprecision(4)
                  << "[DIAG] " << std::setw(5) << i << " | "
                  << state_.rbar[j] << " | "
                  << r << " | "
                  << state_.sigma[j] << " | "
                  << state_.eta[j] << " | "
                  << state_.alpha[j] << " | "
                  << state_.Xs[j] << " | "
                  << state_.omega[j] << "\n";
    }
}

void FiniteTempRelaxer::print_residual_components(int index) const
{
    if (index < 0 || index >= p_.n_gridpoints)
    {
        return;
    }

    const int n = p_.n_gridpoints;
    const double drbar = (p_.rbar_max - p_.rbar_min) / static_cast<double>(n - 1);
    const int i = index;
    const size_t j = static_cast<size_t>(i);
    const double rbar = state_.rbar[j];
    const double r = std::max(1e-12, areal_radius(rbar));
    
    const double ds_drb = d1_rbar(state_.sigma, i, drbar);
    const double de_drb = d1_rbar(state_.eta, i, drbar);
    const double da_drb = d1_rbar(state_.alpha, i, drbar);
    const double dx_drb = d1_rbar(state_.Xs, i, drbar);
    const double ds_dr = d1_r_from_rbar(ds_drb, rbar);
    const double de_dr = d1_r_from_rbar(de_drb, rbar);
    const double da_dr = d1_r_from_rbar(da_drb, rbar);
    const double dx_dr = d1_r_from_rbar(dx_drb, rbar);

    const double alpha = std::max(1e-12, state_.alpha[j]);
    const double Xs = std::max(1e-12, state_.Xs[j]);
    const double sigma = state_.sigma[j];
    const double sigma_sq = sigma * sigma;
    const double eta = state_.eta[j];
    const double X = eta * eta;
    const double temp_local = local_temperature(alpha);

    const double omega = state_.omega[j];
    const double mu_chem = omega / alpha;
    const double Xi = -mu_chem * mu_chem;
    const double P = pressure_Pthermal(temp_local, mu_chem, sigma_sq, X, Xi);
    const ThermoConjugates tc = conjugates_from_pressure(temp_local, mu_chem, sigma_sq, X, Xi);

    const double F = tc.F;

    std::cout << "\n[DIAG] Residual component analysis at grid point " << index << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[DIAG]   rbar = " << rbar << ", r = " << r << ", T_local = " << temp_local << "\n";
    std::cout << "[DIAG]   sigma = " << sigma << ", eta = " << eta << ", alpha = " << alpha << ", Xs = " << Xs << "\n";
    std::cout << "[DIAG]   mu_chem = " << mu_chem << ", Xi = " << Xi << "\n";
    std::cout << "[DIAG]   P_total = " << P << "\n";
    std::cout << "[DIAG]   Thermodynamic conjugates: s=" << tc.s << ", n=" << tc.n << ", F=" << F << ", Upsilon=" << tc.Upsilon << "\n";
    
    std::cout << "\n[DIAG] ETA EQUATION COMPONENTS:\n";
    std::cout << "[DIAG]   de_dr = " << de_dr << "\n";
    std::cout << "[DIAG]   2*eta/r = " << (2.0 * eta / r) << "\n";
    std::cout << "[DIAG]   (da_dr/alpha)*eta = " << ((da_dr / alpha) * eta) << "\n";
    std::cout << "[DIAG]   Xs*sigma*F = " << (Xs * sigma * F) << "\n";
    std::cout << "[DIAG]   Total R.eta = " << (de_dr + 2.0 * eta / r + (da_dr / alpha) * eta + Xs * sigma * F) << "\n";
}