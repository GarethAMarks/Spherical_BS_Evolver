#include "FiniteTempRelaxer.h"
#include "../bose_gas.h"
#include "../bogoliubov.h"

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
        fill_parameter(current_line, "cutoff_percent = ", p_.cutoff_percent, quiet);
        fill_parameter(current_line, "decay_exponent = ", p_.decay_exponent, quiet);
        fill_parameter(current_line, "cph_epsilon = ", p_.cph_epsilon, quiet);
        fill_parameter(current_line, "pressure_interp_power = ", p_.pressure_interp_power, quiet);
        fill_parameter(current_line, "crossover_mode = ", p_.crossover_mode, quiet);
        fill_parameter(current_line, "bogoliubov_gap_width = ", p_.bogoliubov_gap_width, quiet);
        fill_parameter(current_line, "bogoliubov_sigma_gap = ", p_.bogoliubov_sigma_gap, quiet);
        fill_parameter(current_line, "gas_statistics = ", p_.gas_statistics, quiet);
        fill_parameter(current_line, "sigma_phonon = ", p_.sigma_phonon, quiet);
        fill_parameter(current_line, "sigma_gas = ", p_.sigma_gas, quiet);
        fill_parameter(current_line, "crossover_exponent = ", p_.crossover_exponent, quiet);
        fill_parameter(current_line, "sigma_decouple = ", p_.sigma_decouple, quiet);
        fill_parameter(current_line, "decouple_exponent = ", p_.decouple_exponent, quiet);

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
        fill_parameter(current_line, "solver_mode = ", p_.solver_mode, quiet);
        fill_parameter(current_line, "qball_mode = ", p_.qball_mode, quiet);

        fill_parameter(current_line, "output_profile = ", p_.output_profile, quiet);
        fill_parameter(current_line, "save_model_profiles = ", p_.save_model_profiles, quiet);
    }

    if (p_.cutoff_percent < 0.0 || p_.cutoff_percent > 100.0)
    {
        std::cerr << "[WARNING] cutoff_percent = " << p_.cutoff_percent
                  << " is out of range [0, 100]; clamping.\n";
        p_.cutoff_percent = std::min(100.0, std::max(0.0, p_.cutoff_percent));
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
    if (p_.crossover_mode == 2)
        return bogoliubov_pressure(T, mu_chem, sigma_sq, Xi);
    if (p_.crossover_mode == 3)
        return (Xi < 0.0) ? gapped_qp_pressure(T, mu_chem, std::sqrt(-Xi), sigma_sq, X) : 0.0;
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

double FiniteTempRelaxer::bogoliubov_pressure(double T, double mu_chem, double sigma_sq, double Xi) const
{
    // mu_c = sqrt(-Xi) (timelike phase gradient), a = sigma^2 V''(sigma^2), mass m = p_.mu.
    if (T <= 0.0 || Xi >= 0.0)
        return 0.0;
    const double sig2 = std::max(0.0, sigma_sq);
    const double P = bogoliubov::pressure(T, mu_chem, std::sqrt(-Xi), sig2 * ddV(sig2), p_.mu, p_.bogoliubov_gap_width,
                                                bogoliubov::gap_factor(sig2, p_.bogoliubov_sigma_gap));
    // Beyond the Landau limit P is +inf; return a large finite value so line searches reject the trial.
    return std::isfinite(P) ? P : 1e30;
}

double FiniteTempRelaxer::gapped_qp_pressure(double T, double mu_chem, double mu_c, double sigma_sq, double X) const
{
    if (T <= 0.0 || mu_c <= 0.0)
        return 0.0;
    const double sig2 = std::max(0.0, sigma_sq);
    const double m = p_.mu;
    const double sd = std::max(1e-300, p_.sigma_decouple);
    const double xw = std::pow(sig2 / (sd * sd), 0.5 * p_.decouple_exponent);
    const double w = -std::expm1(-xw);
    double P = 0.0;
    if (w > 1e-14)
    {
        const double delta = std::max(0.0, X) / (2.0 * m * sig2 * T);
        const double gam = sig2 * ddV(sig2) / (2.0 * mu_c * T);
        const double eta = (mu_chem - mu_c) / T;
        const double pre = std::pow(mu_c * T / (2.0 * kPi), 1.5) * T;
        double Pqp;
        if (!std::isfinite(delta) || delta > 700.0)
            Pqp = 0.0;
        else
        {
            const double lnf = bogoliubov::lnf_quad(delta, gam, eta);
            if (std::isfinite(lnf))
                Pqp = pre * std::exp(lnf);
            else
            {
                const double fv = bogoliubov::f(delta, gam, eta);
                if (!std::isfinite(fv)) return 1e30; // beyond the Landau limit: reject in line searches
                Pqp = pre * std::max(0.0, fv);
            }
        }
        P += w * Pqp;
    }
    if (1.0 - w > 1e-14)
    {
        double lnF, kap;
        const double x = (mu_chem - m) / std::max(1e-12, T);
        bose_gas::lnF_kappa(std::min(700.0, x), 1, lnF, kap);
        P += (1.0 - w) * std::pow(m / (2.0 * kPi), 1.5) * std::pow(T, 2.5) * std::exp(std::min(700.0, lnF));
    }
    return P;
}

double FiniteTempRelaxer::thermal_Y(double T, double mu_chem, double sigma_sq, double X, double Xi, double r) const
{
    if (p_.crossover_mode != 3 || T < 1e-12 || Xi >= 0.0)
        return 1.0;
    const double mu_c = std::sqrt(-Xi);
    const double sig2 = std::max(1e-300, sigma_sq);
    // step: relative in X, or a change of 1e-4 in delta = X/(2 m sigma^2 T) near X = 0
    const double hX = std::max(1e-4 * X, 1e-4 * 2.0 * p_.mu * sig2 * T);
    double dPdX;
    if (X > hX)
        dPdX = (gapped_qp_pressure(T, mu_chem, mu_c, sig2, X + hX) - gapped_qp_pressure(T, mu_chem, mu_c, sig2, X - hX)) / (2.0 * hX);
    else
        dPdX = (-3.0 * gapped_qp_pressure(T, mu_chem, mu_c, sig2, X) + 4.0 * gapped_qp_pressure(T, mu_chem, mu_c, sig2, X + hX)
                - gapped_qp_pressure(T, mu_chem, mu_c, sig2, X + 2.0 * hX)) / (2.0 * hX);
    return 1.0 - 2.0 * compute_suppression_factor(r) * dPdX;
}

double FiniteTempRelaxer::eta_dlnY_dr(double sigma, double eta, double alpha, double omega, double r,
                                      double ds_dr, double de_dr, double da_dr, double Y0) const
{
    if (p_.crossover_mode != 3 || current_temperature_ < 1e-12)
        return 0.0;
    // Directional derivative of Y along the profile: evaluate Y at the linearly extrapolated local state r +- h.
    const double h = 1e-4 * std::max(1.0, r);
    auto Yat = [&](double sgn) -> double
    {
        const double s = sigma + sgn * h * ds_dr;
        const double e = eta + sgn * h * de_dr;
        const double a = std::max(1e-12, alpha + sgn * h * da_dr);
        const double T = current_temperature_ / a;
        const double mu = omega / a;
        return thermal_Y(T, mu, s * s, e * e, -mu * mu, r + sgn * h);
    };
    const double dY = (Yat(1.0) - Yat(-1.0)) / (2.0 * h);
    return eta * dY / std::max(1e-12, Y0);
}

double FiniteTempRelaxer::compute_suppression_factor(double r) const
{
    if (!p_.use_suppression_factor)
        return 1.0;

    if (r_cutoff_ <= 0.0)
        return 1.0;

    if (r <= r_cutoff_)
        return 1.0;

    const double dr = r - r_cutoff_;
    return std::exp(-16. * (dr * dr) / (r_cutoff_ * r_cutoff_));
}

double FiniteTempRelaxer::pressure_Pthermal(double T, double mu_chem, double sigma_sq, double X, double Xi, double r) const
{
    // Pthermal Eq. (notes): P = scalar piece + suppression * thermal piece
    const double suppression_factor = compute_suppression_factor(r);
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
    double Xi,
    double r) const
{
    ThermoConjugates out;

    const double fd = std::max(1e-7, p_.jacobian_fd_step);
    const double dT = fd * std::max(1.0, std::abs(T));
    const double dmu = fd * std::max(1.0, std::abs(mu_chem));
    const double ds2 = fd * std::max(1.0, std::abs(sigma_sq));
    const double dXi = fd * std::max(1.0, std::abs(Xi));

    const double pTp = pressure_Pthermal(T + dT, mu_chem, sigma_sq, X, Xi, r);
    const double pTm = pressure_Pthermal(T - dT, mu_chem, sigma_sq, X, Xi, r);
    out.s = (pTp - pTm) / (2.0 * dT);

    const double pMup = pressure_Pthermal(T, mu_chem + dmu, sigma_sq, X, Xi, r);
    const double pMum = pressure_Pthermal(T, mu_chem - dmu, sigma_sq, X, Xi, r);
    out.n = (pMup - pMum) / (2.0 * dmu);

    const double pS2p = pressure_Pthermal(T, mu_chem, sigma_sq + ds2, X, Xi, r);
    const double pS2m = pressure_Pthermal(T, mu_chem, std::max(0.0, sigma_sq - ds2), X, Xi, r);
    out.F = 2.0 * (pS2p - pS2m) / (2.0 * ds2);

    const double pXip = pressure_Pthermal(T, mu_chem, sigma_sq, X, Xi + dXi, r);
    const double pXim = pressure_Pthermal(T, mu_chem, sigma_sq, X, Xi - dXi, r);
    out.Upsilon = -2.0 * (pXip - pXim) / (2.0 * dXi);

    return out;
}

FiniteTempRelaxer::ThermoConjugates FiniteTempRelaxer::conjugates_analytic(
    double T,
    double mu_chem,
    double sigma_sq,
    double X,
    double Xi,
    double r) const
{
    ThermoConjugates out;

    const double sig2 = std::max(0.0, sigma_sq);
    // Scalar part: P_scalar = -0.5*(X + sig2*Xi + V(sig2))
    //   dP_scalar/dsig2 = -0.5*(Xi + V'(sig2))
    //   dP_scalar/dXi   = -0.5*sig2
    // d_potential_d_sigma2 expects the field amplitude sigma (not sig2); callers pass sqrt(sig2).
    const double Vp = d_potential_d_sigma2(std::sqrt(sig2));
    const double dPsc_ds2 = -0.5 * (Xi + Vp);
    const double dPsc_dXi = -0.5 * sig2;

    const double sup = compute_suppression_factor(r);

    if (p_.crossover_mode == 3 && T >= 1e-12)
    {
        // Recommended EoS R: central differences of the smooth quadrature (lnf_quad), one-sided near sigma^2 = 0.
        const double mu_c = (Xi < 0.0) ? std::sqrt(-Xi) : 0.0;
        auto Pq = [&](double TT, double mm, double mc, double ss) -> double
        {
            return gapped_qp_pressure(TT, mm, mc, ss, X);
        };
        const double hT = 1e-4 * T, hm = 1e-4 * T, hc = 1e-4 * T;
        const double hs = 1e-4 * std::max(sig2, 1e-300);
        const double dP_dT = (Pq(T + hT, mu_chem, mu_c, sig2) - Pq(T - hT, mu_chem, mu_c, sig2)) / (2.0 * hT);
        const double dP_dmu = (Pq(T, mu_chem + hm, mu_c, sig2) - Pq(T, mu_chem - hm, mu_c, sig2)) / (2.0 * hm);
        const double dP_dmuc = (mu_c > 2.0 * hc)
            ? (Pq(T, mu_chem, mu_c + hc, sig2) - Pq(T, mu_chem, mu_c - hc, sig2)) / (2.0 * hc) : 0.0;
        const double dP_ds2 = (sig2 > 0.0)
            ? (Pq(T, mu_chem, mu_c, sig2 + hs) - Pq(T, mu_chem, mu_c, sig2 - hs)) / (2.0 * hs) : 0.0;
        const double dP_dXi = (mu_c > 0.0) ? -dP_dmuc / (2.0 * mu_c) : 0.0;
        out.s = sup * dP_dT;
        out.n = sup * dP_dmu;
        out.F = 2.0 * (dPsc_ds2 + sup * dP_ds2);
        out.Upsilon = -2.0 * (dPsc_dXi + sup * dP_dXi);
        out.Y = thermal_Y(T, mu_chem, sig2, X, Xi, r);
        return out;
    }

    if (p_.crossover_mode == 2 && T >= 1e-12)
    {
        // Bogoliubov EoS: central differences of the (smooth) quadrature. sigma^2 enters only through
        // a = sigma^2 V'', so dP/dsigma^2 = (dP/da) (V'' + sigma^2 V'''); the a-step is relative to a, one-sided near a = 0.
        const double mu_c = (Xi < 0.0) ? std::sqrt(-Xi) : 0.0;
        const double Vpp = ddV(sig2);
        double Vppp = 0.0;
        if (p_.solitonic)
        {
            const double scale2 = std::max(1e-12, p_.sigma_scale * p_.sigma_scale);
            Vppp = 24.0 * p_.mu * p_.mu / (scale2 * scale2);
        }
        const double a = sig2 * Vpp;
        double dg_ds2 = 0.0;
        const double gf = bogoliubov::gap_factor(sig2, p_.bogoliubov_sigma_gap, &dg_ds2);
        auto Pb = [&](double TT, double mm, double aa, double mc, double gg) -> double
        {
            if (TT <= 0.0 || mc <= 0.0) return 0.0;
            const double P = bogoliubov::pressure(TT, mm, mc, std::max(0.0, aa), p_.mu, p_.bogoliubov_gap_width, gg);
            return std::isfinite(P) ? P : 1e30;
        };
        const double hT = 1e-4 * T;
        const double hm = 1e-4 * T;
        const double ha = 1e-4 * std::max(a, 1e-3 * mu_c * T);
        const double hc = 1e-4 * T;
        const double dP_dT = (Pb(T + hT, mu_chem, a, mu_c, gf) - Pb(T - hT, mu_chem, a, mu_c, gf)) / (2.0 * hT);
        const double dP_dmu = (Pb(T, mu_chem + hm, a, mu_c, gf) - Pb(T, mu_chem - hm, a, mu_c, gf)) / (2.0 * hm);
        const double dP_da = (a > ha)
            ? (Pb(T, mu_chem, a + ha, mu_c, gf) - Pb(T, mu_chem, a - ha, mu_c, gf)) / (2.0 * ha)
            : (-3.0 * Pb(T, mu_chem, a, mu_c, gf) + 4.0 * Pb(T, mu_chem, a + ha, mu_c, gf) - Pb(T, mu_chem, a + 2.0 * ha, mu_c, gf)) / (2.0 * ha);
        const double dP_dmuc = (mu_c > 2.0 * hc)
            ? (Pb(T, mu_chem, a, mu_c + hc, gf) - Pb(T, mu_chem, a, mu_c - hc, gf)) / (2.0 * hc)
            : 0.0;
        const double hg = 1e-4;
        const double dP_dg = (dg_ds2 != 0.0)
            ? (Pb(T, mu_chem, a, mu_c, std::min(1.0, gf + hg)) - Pb(T, mu_chem, a, mu_c, std::max(0.0, gf - hg)))
                / (std::min(1.0, gf + hg) - std::max(0.0, gf - hg))
            : 0.0;
        // dmu_c/dXi = -1/(2 mu_c)
        const double dP_dXi = (mu_c > 0.0) ? -dP_dmuc / (2.0 * mu_c) : 0.0;
        out.s = sup * dP_dT;
        out.n = sup * dP_dmu;
        out.F = 2.0 * (dPsc_ds2 + sup * (dP_da * (Vpp + sig2 * Vppp) + dP_dg * dg_ds2));
        out.Upsilon = -2.0 * (dPsc_dXi + sup * dP_dXi);
        (void)X;
        return out;
    }

    if (T < 1e-12)
    {
        // At T = 0 the thermal piece is defined to be identically zero.
        out.s = 0.0;
        out.n = 0.0;
        out.F = 2.0 * dPsc_ds2;
        out.Upsilon = -2.0 * dPsc_dXi;
        return out;
    }

    // ---------- phonon sector ----------
    const double eps_cph = std::max(0.0, p_.cph_epsilon);
    const double Vpp = ddV(sig2);
    // V''' is zero for the plain quartic; nonzero for the solitonic shape.
    double Vppp = 0.0;
    if (p_.solitonic)
    {
        const double scale2 = std::max(1e-12, p_.sigma_scale * p_.sigma_scale);
        Vppp = 24.0 * p_.mu * p_.mu / (scale2 * scale2);
    }
    const double g = sig2 * Vpp;
    const double num = g + eps_cph;
    const double den_raw = g - 2.0 * Xi;
    const double den = (std::abs(den_raw) < 1e-12) ? ((den_raw >= 0.0) ? 1e-12 : -1e-12) : den_raw;
    const bool den_clamped = (std::abs(den_raw) < 1e-12);
    const double cph2_raw = num / den;
    const bool cph2_floored = (cph2_raw < 1e-12);
    const double cph2 = cph2_floored ? 1e-12 : cph2_raw;
    const double cph = std::sqrt(cph2);

    const double Xi_safe = (std::abs(Xi) < 1e-12) ? ((Xi >= 0.0) ? 1e-12 : -1e-12) : Xi;
    const bool Xi_clamped = (std::abs(Xi) < 1e-12);
    const double mu2 = mu_chem * mu_chem;
    const double bracket_raw = 1.0 + (1.0 - cph2) * mu2 / Xi_safe;
    const double bracket = (std::abs(bracket_raw) < 1e-12) ? ((bracket_raw >= 0.0) ? 1e-12 : -1e-12) : bracket_raw;

    const double T2 = T * T;
    const double T4 = T2 * T2;
    const double inv_bracket = 1.0 / bracket;
    const double inv_bracket2 = inv_bracket * inv_bracket;
    constexpr double kPhonCoef = 3.1415926535897932384626433832795 * 3.1415926535897932384626433832795 / 90.0;
    const double A = kPhonCoef * cph * T4 * inv_bracket2;

    // Partial derivatives of cph^2 wrt sig2 and Xi (zero if floored or den clamped).
    double dcph2_ds2 = 0.0;
    double dcph2_dXi = 0.0;
    if (!cph2_floored && !den_clamped)
    {
        const double dg_ds2 = Vpp + sig2 * Vppp;
        const double inv_den = 1.0 / den;
        const double inv_den2 = inv_den * inv_den;
        dcph2_ds2 = dg_ds2 * (-2.0 * Xi - eps_cph) * inv_den2;
        dcph2_dXi = 2.0 * num * inv_den2;
    }
    // d(cph)/dx = d(cph^2)/dx / (2 cph).
    const double two_cph = 2.0 * cph;
    const double dcph_ds2 = dcph2_ds2 / two_cph;
    const double dcph_dXi = dcph2_dXi / two_cph;

    // Partial derivatives of bracket.
    const double dbr_dmu = 2.0 * (1.0 - cph2) * mu_chem / Xi_safe;
    const double dbr_ds2 = -dcph2_ds2 * mu2 / Xi_safe;
    double dbr_dXi = -dcph2_dXi * mu2 / Xi_safe;
    if (!Xi_clamped)
    {
        dbr_dXi -= (1.0 - cph2) * mu2 / (Xi_safe * Xi_safe);
    }

    // Partial derivatives of A = kPhonCoef * cph * T^4 / bracket^2.
    const double dA_dT = 4.0 * A / T;
    const double dA_dmu = A * (-2.0 * dbr_dmu * inv_bracket);
    const double dA_ds2 = (cph > 1e-20)
        ? A * (dcph_ds2 / cph - 2.0 * dbr_ds2 * inv_bracket)
        : 0.0;
    const double dA_dXi = (cph > 1e-20)
        ? A * (dcph_dXi / cph - 2.0 * dbr_dXi * inv_bracket)
        : 0.0;

    // ---------- gas sector ----------
    // B = (mu_mass/(2pi))^(3/2) * T^(5/2) * exp((mu_chem - mu_mass)/T)
    constexpr double kTwoPi = 2.0 * 3.1415926535897932384626433832795;
    const double mass = p_.mu;
    const double gas_prefac = std::pow(mass / kTwoPi, 1.5);
    const double T_sqrt = std::sqrt(T);
    const double T52 = T2 * T_sqrt;             // T^(5/2)
    // Boltzmann (gas_statistics = 0) or Bose-Einstein (1): B = prefac T^{5/2} F(x), x = (mu - m)/T, kappa = F'/F.
    // Cap ln B growth to avoid overflow when a line-search trial drives alpha toward zero (mu_chem = omega/alpha
    // -> infinity). The downstream residual will still be huge and the line search will reject the trial.
    const double gas_x = (mu_chem - mass) / std::max(1e-12, T);
    double gas_lnF, gas_kappa;
    bose_gas::lnF_kappa(std::min(700.0, gas_x), p_.gas_statistics, gas_lnF, gas_kappa);
    const double B = gas_prefac * T52 * std::exp(std::min(700.0, gas_lnF));

    const double dB_dT = B * (2.5 - gas_x * gas_kappa) / T;
    const double dB_dmu = B * gas_kappa / T;
    // dB/ds2 = 0, dB/dXi = 0.

    if (p_.crossover_mode == 1)
    {
        // Convex blend P_th = w(sig2) A + G(sig2) B, w = 1 - exp(-x_ph), G = exp(-x_g), x = (sig2/sigma_*^2)^{k/2}
        // (same as thermal_pressure_term).
        const double kh = 0.5 * p_.crossover_exponent;
        const double x_ph = std::pow(sig2 / std::max(1e-300, p_.sigma_phonon * p_.sigma_phonon), kh);
        const double x_g = std::pow(sig2 / std::max(1e-300, p_.sigma_gas * p_.sigma_gas), kh);
        const double wph = -std::expm1(-x_ph);
        const double G = std::exp(-x_g);
        const double inv_s2 = (sig2 > 1e-300) ? 1.0 / sig2 : 0.0;
        const double dw_ds2 = std::exp(-x_ph) * kh * x_ph * inv_s2;
        const double dG_ds2 = -G * kh * x_g * inv_s2;
        const double Aw = (wph > 0.0) ? wph : 0.0;
        out.s = sup * (Aw * dA_dT + G * dB_dT);
        out.n = sup * (Aw * dA_dmu + G * dB_dmu);
        out.F = 2.0 * (dPsc_ds2 + sup * (Aw * dA_ds2 + ((wph > 0.0) ? dw_ds2 * A : 0.0) + dG_ds2 * B));
        out.Upsilon = -2.0 * (dPsc_dXi + sup * Aw * dA_dXi);
        (void)X;
        return out;
    }

    // ---------- interpolation: P_th = (A^-q + B^-q)^(-1/q) ----------
    // Factor out the smaller of A, B to avoid overflow in the intermediates
    // when either phonon or gas pressure underflows (very common when alpha
    // -> 0 during a line-search probe, which makes mu_chem -> infinity).
    //
    // Let m = min(A, B), M = max(A, B), r = m/M <= 1, w = r^q.
    // Then P_th = m * (1 + w)^(-1/q), and
    //   dP_th/dx = u/(1+w) * (dm/dx + m * w/(1+w) * dM/dx/M)
    //            = u/(1+w) * dm/dx + u * m * w / ((1+w)*(1+w)) * dM/dx / M
    // where u = (1 + w)^(-1/q).
    const double q = std::max(1e-12, p_.pressure_interp_power);
    const double A_pos = std::max(1e-300, A);
    const double B_pos = std::max(1e-300, B);
    const bool A_smaller = A_pos <= B_pos;
    const double m_val = A_smaller ? A_pos : B_pos;
    const double M_val = A_smaller ? B_pos : A_pos;
    const double ratio = m_val / M_val;         // in (0, 1]
    const double w = std::pow(ratio, q);        // in (0, 1]
    const double one_plus_w = 1.0 + w;
    const double u = std::pow(one_plus_w, -1.0 / q);
    const double P_th = m_val * u;

    const double coef_m = u / one_plus_w;             // for dm/dx
    const double coef_M = m_val * u * w / one_plus_w; // for dM/dx / M

    auto dPth = [&](double dA_dx, double dB_dx) -> double
    {
        const double dm_dx = A_smaller ? dA_dx : dB_dx;
        const double dM_dx = A_smaller ? dB_dx : dA_dx;
        return coef_m * dm_dx + coef_M * dM_dx / M_val;
    };

    const double dPth_dT   = dPth(dA_dT,  dB_dT);
    const double dPth_dmu  = dPth(dA_dmu, dB_dmu);
    const double dPth_ds2  = dPth(dA_ds2, 0.0);
    const double dPth_dXi  = dPth(dA_dXi, 0.0);

    // Apply suppression (independent of thermodynamic variables).
    out.s        = sup * dPth_dT;
    out.n        = sup * dPth_dmu;
    out.F        = 2.0 * (dPsc_ds2 + sup * dPth_ds2);
    out.Upsilon  = -2.0 * (dPsc_dXi + sup * dPth_dXi);
    (void)X;    // scalar partial wrt X is -0.5, but X is a dependent variable (eta^2), not needed.
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
        shooter.qball_mode = p_.qball_mode;
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
    shooter.qball_mode = p_.qball_mode;
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

    if (p_.qball_mode)
    {
        for (int i = 0; i < n; ++i)
        {
            state_.alpha[static_cast<size_t>(i)] = 1.0;
            state_.Xs[static_cast<size_t>(i)] = 1.0;
        }
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

std::vector<double> FiniteTempRelaxer::pack_state_interleaved() const
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    const size_t fields = active_field_count();
    std::vector<double> x(fields * n);
    for (size_t i = 0; i < n; ++i)
    {
        const size_t base = fields * i;
        x[base] = state_.sigma[i];
        x[base + 1] = state_.eta[i];
        x[base + 2] = state_.alpha[i];
        x[base + 3] = state_.Xs[i];
        x[base + 4] = state_.omega[i];
        if (p_.increment_omega)
        {
            x[base + 5] = state_.rho[i];
        }
    }
    return x;
}

void FiniteTempRelaxer::unpack_state_interleaved(const std::vector<double>& x)
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    const size_t fields = active_field_count();
    for (size_t i = 0; i < n; ++i)
    {
        const size_t base = fields * i;
        state_.sigma[i] = x[base];
        state_.eta[i] = x[base + 1];
        state_.alpha[i] = std::max(1e-12, x[base + 2]);
        state_.Xs[i] = std::max(1e-12, x[base + 3]);
        state_.omega[i] = x[base + 4];
        if (p_.increment_omega)
        {
            state_.rho[i] = x[base + 5];
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
        const double P = pressure_Pthermal(temp_local, mu_chem, sigma_sq, X, Xi, r);
        const ThermoConjugates tc = conjugates_analytic(temp_local, mu_chem, sigma_sq, X, Xi, r);

        const double F = tc.F;
        const double source_X = mu_chem * mu_chem * tc.Upsilon + mu_chem * tc.n + tc.s * temp_local - P;

        const double Ysig = tc.Y;
        const double yterm = eta_dlnY_dr(sigma, eta, alpha, omega, r, ds_dr, de_dr, da_dr, Ysig);
        R.sigma[j] = ds_dr - Xs * eta;
        R.eta[j] = de_dr + 2.0 * eta / r + (da_dr / alpha) * eta + yterm + Xs * sigma * F / Ysig;
        R.alpha[j] = da_dr / alpha - ((Xs * Xs - 1.0) / (2.0 * r) + 4.0 * kPi * r * Xs * Xs * (Ysig * eta * eta + P));
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

    if (p_.qball_mode)
    {
        for (int i = 0; i < n; ++i)
        {
            const size_t j = static_cast<size_t>(i);
            R.alpha[j] = state_.alpha[j] - 1.0;
            R.Xs[j]    = state_.Xs[j]    - 1.0;
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

std::vector<double> FiniteTempRelaxer::pack_residual_vector_interleaved(const Residuals& R) const
{
    const size_t n = static_cast<size_t>(p_.n_gridpoints);
    const size_t fields = active_field_count();
    std::vector<double> v(fields * n);
    for (size_t i = 0; i < n; ++i)
    {
        const size_t base = fields * i;
        v[base] = R.sigma[i];
        v[base + 1] = R.eta[i];
        v[base + 2] = R.alpha[i];
        v[base + 3] = R.Xs[i];
        v[base + 4] = R.omega[i];
        if (p_.increment_omega)
        {
            v[base + 5] = R.rho[i];
        }
    }
    return v;
}

double FiniteTempRelaxer::residual_norm(const Residuals& R) const
{
    if (p_.solver_mode == static_cast<int>(SolverMode::BlockBanded))
    {
        return vector_norm_inf(pack_residual_vector_interleaved(R));
    }
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

    const bool banded_mode = (p_.solver_mode == static_cast<int>(SolverMode::BlockBanded));
    const size_t fields = active_field_count();
    const size_t half_bandwidth = banded_mode ? (3U * fields - 1U) : 0U;

    if (banded_mode)
    {
        const size_t stride = 2U * half_bandwidth + 1U;
        J.assign(stride * nvar, 0.0);
    }
    else
    {
        J.assign(nvar * nvar, 0.0);
    }

    auto band_index = [nvar, half_bandwidth](size_t row, size_t col) -> size_t
    {
        return (half_bandwidth + row - col) * nvar + col;
    };

    if (!banded_mode)
    {
        // Dense mode: perturb each column, recompute all residuals, store full column.
        for (size_t col = 0; col < nvar; ++col)
        {
            std::vector<double> x_pert = x_base;
            const double step = eps * std::max(1.0, std::abs(x_base[col]));
            x_pert[col] += step;
            const_cast<FiniteTempRelaxer*>(this)->unpack_state(x_pert);
            const Residuals Rp = compute_residuals();
            const std::vector<double> r_pert = pack_residual_vector(Rp);
            for (size_t row = 0; row < nvar; ++row)
                J[row * nvar + col] = (r_pert[row] - r_base[row]) / step;
        }
        const_cast<FiniteTempRelaxer*>(this)->unpack_state(x_base);
        return;
    }

    // Banded mode: local-window optimisation.
    //
    // The residual at interior grid point k depends on field values at k-2..k+2 (5-point stencil).
    // Boundary residuals at k=0 depend on points 0 and 1; at k=last on points last and last-1.
    // Therefore perturbing the unknown at grid point i_col can only change residuals at grid
    // points in [i_col-2, i_col+2] (plus boundary handling), so we recompute only that window
    // instead of calling the full compute_residuals() over all n points.
    //
    // Exception: in increment_omega mode, R.rho[last] = particle_number_from_rho() - target
    // integrates over all rho values, so any rho perturbation also affects that one entry.

    const int n    = p_.n_gridpoints;
    const int last = n - 1;
    const double drbar = (p_.rbar_max - p_.rbar_min) / static_cast<double>(n - 1);

    // Parallelism: each thread perturbs its own copy of the five (or six) field vectors.
    // state_ itself is untouched — safe to read concurrently from residual helpers.
    // Each column writes to a disjoint slice of J[], so no reduction is needed.
    #pragma omp parallel
    {
        RelaxerState ls = state_;                 // thread-local scratch copy of state
        std::vector<double> r_window = r_base;    // thread-local window residual buffer

        #pragma omp for schedule(static)
        for (long col_signed = 0; col_signed < static_cast<long>(nvar); ++col_signed)
        {
            const size_t col   = static_cast<size_t>(col_signed);
            const size_t i_col = col / fields;
            const size_t f_col = col % fields;
            const double orig  = x_base[col];
            const double step  = eps * std::max(1.0, std::abs(orig));

            // Perturb one state variable in the thread-local copy.
            double* perturbed_ptr = nullptr;
            switch (f_col)
            {
                case 0: ls.sigma[i_col] = orig + step; perturbed_ptr = &ls.sigma[i_col]; break;
                case 1: ls.eta[i_col]   = orig + step; perturbed_ptr = &ls.eta[i_col];   break;
                case 2: ls.alpha[i_col] = std::max(1e-12, orig + step); perturbed_ptr = &ls.alpha[i_col]; break;
                case 3: ls.Xs[i_col]    = std::max(1e-12, orig + step); perturbed_ptr = &ls.Xs[i_col];    break;
                case 4: ls.omega[i_col] = orig + step; perturbed_ptr = &ls.omega[i_col]; break;
                default: ls.rho[i_col]  = orig + step; perturbed_ptr = &ls.rho[i_col];   break;
            }

            const int win_min = std::max(0,    (int)i_col - 2);
            const int win_max = std::min(last, (int)i_col + 2);

            for (int k = win_min; k <= win_max; ++k)
            {
                const size_t kj   = static_cast<size_t>(k);
                const size_t base = fields * kj;

                if (k == 0)
                {
                    r_window[base + 0] = ls.sigma[1] - ls.sigma[0];
                    r_window[base + 1] = ls.eta[0];
                    r_window[base + 2] = ls.alpha[1] - ls.alpha[0];
                    r_window[base + 3] = ls.Xs[0] - 1.0;
                    r_window[base + 4] = ls.omega[0] - p_.mu_chem_infty;
                    if (fields > 5)
                        r_window[base + 5] = ls.rho[1] - ls.rho[0];
                }
                else if (k == last)
                {
                    r_window[base + 0] = ls.sigma[kj];
                    r_window[base + 1] = ls.eta[kj];
                    r_window[base + 2] = ls.alpha[kj] - 1.0;
                    r_window[base + 3] = ls.Xs[kj] - 1.0;
                    if (!p_.increment_omega)
                    {
                        r_window[base + 4] = ls.omega[kj] - p_.mu_chem_infty;
                    }
                    else
                    {
                        r_window[base + 4] = ls.omega[kj] - ls.omega[kj - 1];
                        // R.rho[last] uses the local rho view, consistent with the perturbation.
                        if (fields > 5)
                        {
                            const int rn = p_.n_gridpoints;
                            const int li = std::max(1, rn - 2);
                            const double rbar_last = ls.rbar[static_cast<size_t>(li)];
                            const double r_last = std::max(1e-12, areal_radius(rbar_last));
                            const double drho_drb_last = d1_rbar(ls.rho, li, drbar);
                            const double drho_dr_last  = d1_r_from_rbar(drho_drb_last, rbar_last);
                            r_window[base + 5] = 8.0 * 3.1415926535897932384626433832795 * r_last * r_last * drho_dr_last - target_particle_number_;
                        }
                    }
                }
                else
                {
                    const double rbar   = ls.rbar[kj];
                    const double r      = std::max(1e-12, areal_radius(rbar));
                    const double ds_drb = d1_rbar(ls.sigma, k, drbar);
                    const double de_drb = d1_rbar(ls.eta,   k, drbar);
                    const double da_drb = d1_rbar(ls.alpha, k, drbar);
                    const double dx_drb = d1_rbar(ls.Xs,    k, drbar);
                    const double ds_dr  = d1_r_from_rbar(ds_drb, rbar);
                    const double de_dr  = d1_r_from_rbar(de_drb, rbar);
                    const double da_dr  = d1_r_from_rbar(da_drb, rbar);
                    const double dx_dr  = d1_r_from_rbar(dx_drb, rbar);

                    const double alpha    = std::max(1e-12, ls.alpha[kj]);
                    const double Xs       = std::max(1e-12, ls.Xs[kj]);
                    const double sigma    = ls.sigma[kj];
                    const double sigma_sq = sigma * sigma;
                    const double eta      = ls.eta[kj];
                    const double X        = eta * eta;
                    const double temp_local = local_temperature(alpha);
                    const double omega    = ls.omega[kj];
                    const double mu_chem  = omega / alpha;
                    const double Xi       = -mu_chem * mu_chem;
                    const double P        = pressure_Pthermal(temp_local, mu_chem, sigma_sq, X, Xi, r);
                    const ThermoConjugates tc = conjugates_analytic(temp_local, mu_chem, sigma_sq, X, Xi, r);

                    const double F        = tc.F;
                    const double source_X = mu_chem * mu_chem * tc.Upsilon + mu_chem * tc.n + tc.s * temp_local - P;

                    const double Ysig = tc.Y;
                    const double yterm = eta_dlnY_dr(sigma, eta, alpha, omega, r, ds_dr, de_dr, da_dr, Ysig);
                    r_window[base + 0] = ds_dr - Xs * eta;
                    r_window[base + 1] = de_dr + 2.0 * eta / r + (da_dr / alpha) * eta + yterm + Xs * sigma * F / Ysig;
                    r_window[base + 2] = da_dr / alpha - ((Xs * Xs - 1.0) / (2.0 * r) + 4.0 * kPi * r * Xs * Xs * (Ysig * eta * eta + P));
                    r_window[base + 3] = dx_dr / Xs - ((1.0 - Xs * Xs) / (2.0 * r) + 4.0 * kPi * r * Xs * Xs * source_X);

                    if (!p_.increment_omega)
                    {
                        r_window[base + 4] = ls.omega[kj] - p_.mu_chem_infty;
                    }
                    else
                    {
                        r_window[base + 4] = ls.omega[kj] - ls.omega[kj - 1];
                        if (fields > 5)
                        {
                            const double drho_drb   = d1_rbar(ls.rho, k, drbar);
                            const double d2rho_drb2 = d2_rbar(ls.rho, k, drbar);
                            const double drho_dr    = d1_r_from_rbar(drho_drb, rbar);
                            const double d2rho_dr2  = d2_r_from_rbar(drho_drb, d2rho_drb2, rbar);
                            const double omega_safe = std::max(1e-12, std::abs(omega));
                            const double source_rho = tc.Upsilon / (alpha * alpha) + tc.n / (alpha * omega_safe);
                            r_window[base + 5] = d2rho_dr2 + (2.0 / r + (da_dr / alpha) - (dx_dr / Xs)) * drho_dr - source_rho;
                        }
                    }
                }

                if (p_.qball_mode)
                {
                    r_window[base + 2] = ls.alpha[kj] - 1.0;
                    r_window[base + 3] = ls.Xs[kj]    - 1.0;
                }
            }

            // In increment_omega mode, R.rho[last] is a global integral of rho and depends on every
            // interior rho sample. If the current column perturbs rho at i_col and the window does
            // not already cover `last`, recompute the global entry using the local rho view.
            if (fields > 5 && f_col == 5 && (int)i_col + 2 < last)
            {
                const size_t base_last = fields * static_cast<size_t>(last);
                const int rn = p_.n_gridpoints;
                const int li = std::max(1, rn - 2);
                const double rbar_last = ls.rbar[static_cast<size_t>(li)];
                const double r_last = std::max(1e-12, areal_radius(rbar_last));
                const double drho_drb_last = d1_rbar(ls.rho, li, drbar);
                const double drho_dr_last  = d1_r_from_rbar(drho_drb_last, rbar_last);
                r_window[base_last + 5] = 8.0 * 3.1415926535897932384626433832795 * r_last * r_last * drho_dr_last - target_particle_number_;
            }

            const size_t row_begin = (col > half_bandwidth) ? (col - half_bandwidth) : 0U;
            const size_t row_end   = std::min(nvar - 1U, col + half_bandwidth);
            for (size_t row = row_begin; row <= row_end; ++row)
                J[band_index(row, col)] = (r_window[row] - r_base[row]) / step;

            // Restore the one perturbed state variable to the base value for the next iteration.
            *perturbed_ptr = orig;

            // Restore r_window entries for the touched window.
            for (int k = win_min; k <= win_max; ++k)
            {
                const size_t base = fields * static_cast<size_t>(k);
                for (size_t f = 0; f < fields; ++f)
                    r_window[base + f] = r_base[base + f];
            }
            if (fields > 5 && f_col == 5 && (int)i_col + 2 < last)
            {
                const size_t base_last = fields * static_cast<size_t>(last);
                r_window[base_last + 5] = r_base[base_last + 5];
            }
        }
    }
    // state_ was never mutated; nothing to restore.
}

bool FiniteTempRelaxer::solve_banded_system(std::vector<double> A,
                                            std::vector<double> b,
                                            std::vector<double>& x,
                                            size_t half_bandwidth) const
{
    return solve_compact_banded_system(std::move(A), std::move(b), x, half_bandwidth);
}

bool FiniteTempRelaxer::solve_compact_banded_system(std::vector<double> A,
                                                    std::vector<double> b,
                                                    std::vector<double>& x,
                                                    size_t half_bandwidth) const
{
    const size_t n = b.size();
    const size_t stride = 2U * half_bandwidth + 1U;
    x.assign(n, 0.0);

    auto band_index = [stride, half_bandwidth, n](size_t row, size_t col) -> size_t
    {
        return (half_bandwidth + row - col) * n + col;
    };

    auto band_get = [&](size_t row, size_t col) -> double
    {
        if (col > row + half_bandwidth || row > col + half_bandwidth)
        {
            return 0.0;
        }
        return A[band_index(row, col)];
    };

    auto band_set = [&](size_t row, size_t col, double value)
    {
        if (col > row + half_bandwidth || row > col + half_bandwidth)
        {
            return;
        }
        A[band_index(row, col)] = value;
    };

    for (size_t k = 0; k < n; ++k)
    {
        const double diag = band_get(k, k);
        if (std::abs(diag) < 1e-18 || !std::isfinite(diag))
        {
            return false;
        }

        const size_t row_end = std::min(n - 1, k + half_bandwidth);
        const size_t col_end = std::min(n - 1, k + half_bandwidth);

        for (size_t i = k + 1; i <= row_end; ++i)
        {
            const double entry = band_get(i, k);
            if (std::abs(entry) < 1e-30)
            {
                continue;
            }

            const double factor = entry / diag;
            band_set(i, k, 0.0);
            const size_t local_end = std::min(col_end, i + half_bandwidth);
            for (size_t j = k + 1; j <= local_end; ++j)
            {
                band_set(i, j, band_get(i, j) - factor * band_get(k, j));
            }
            b[i] -= factor * b[k];
        }
    }

    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = b[static_cast<size_t>(i)];
        const size_t col_end = std::min(n - 1, static_cast<size_t>(i) + half_bandwidth);
        for (size_t j = static_cast<size_t>(i) + 1; j <= col_end; ++j)
        {
            sum -= band_get(static_cast<size_t>(i), j) * x[j];
        }

        const double diag = band_get(static_cast<size_t>(i), static_cast<size_t>(i));
        if (std::abs(diag) < 1e-18 || !std::isfinite(diag))
        {
            return false;
        }

        x[static_cast<size_t>(i)] = sum / diag;
    }

    return true;
}

bool FiniteTempRelaxer::solve_bordered_banded_system(const std::vector<double>& J,
                                                     const std::vector<double>& rhs,
                                                     std::vector<double>& x) const
{
    const size_t n = rhs.size();
    if (n < 2)
    {
        return false;
    }

    const size_t m = n - 1;
    const size_t fields = active_field_count();
    const size_t half_bandwidth = 3U * fields - 1U;

    std::vector<double> B((2U * half_bandwidth + 1U) * m, 0.0);
    std::vector<double> rhs_y(m, 0.0);
    std::vector<double> bvec(m, 0.0);
    std::vector<double> cvec(m, 0.0);

    auto band_index = [half_bandwidth, m](size_t row, size_t col) -> size_t
    {
        return (half_bandwidth + row - col) * m + col;
    };

    auto J_get = [half_bandwidth, n, &J](size_t row, size_t col) -> double
    {
        if (col > row + half_bandwidth || row > col + half_bandwidth)
        {
            return 0.0;
        }
        return J[(half_bandwidth + row - col) * n + col];
    };

    auto band_set = [&](size_t row, size_t col, double value)
    {
        if (col > row + half_bandwidth || row > col + half_bandwidth)
        {
            return;
        }
        B[band_index(row, col)] = value;
    };

    for (size_t i = 0; i < m; ++i)
    {
        rhs_y[i] = rhs[i + 1];
        bvec[i] = J_get(i + 1, 0);
        cvec[i] = J_get(0, i + 1);

        const size_t row_begin = (i > half_bandwidth) ? (i - half_bandwidth) : 0U;
        const size_t row_end = std::min(m - 1U, i + half_bandwidth);
        for (size_t j = row_begin; j <= row_end; ++j)
        {
            band_set(i, j, J_get(i + 1, j + 1));
        }
    }

    std::vector<double> y0;
    std::vector<double> y1;
    if (!solve_compact_banded_system(B, rhs_y, y0, half_bandwidth))
    {
        return false;
    }
    if (!solve_compact_banded_system(B, bvec, y1, half_bandwidth))
    {
        return false;
    }

    const double a = J[0];
    double cTy0 = 0.0;
    double cTy1 = 0.0;
    for (size_t i = 0; i < m; ++i)
    {
        cTy0 += cvec[i] * y0[i];
        cTy1 += cvec[i] * y1[i];
    }

    const double schur = a - cTy1;
    if (std::abs(schur) < 1e-18 || !std::isfinite(schur))
    {
        return false;
    }

    x.assign(n, 0.0);
    x[0] = (rhs[0] - cTy0) / schur;
    for (size_t i = 0; i < m; ++i)
    {
        x[i + 1] = y0[i] - y1[i] * x[0];
    }

    return true;
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
    const bool banded_mode = (p_.solver_mode == static_cast<int>(SolverMode::BlockBanded));
    const size_t nvar = active_field_count() * static_cast<size_t>(p_.n_gridpoints);
    std::vector<double> x = banded_mode ? pack_state_interleaved() : pack_state();
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
        const std::vector<double> r_base = banded_mode ? pack_residual_vector_interleaved(R)
                                                       : pack_residual_vector(R);
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
        bool banded_ok = false;
        if (banded_mode)
        {
            // Direct banded solve on the full nvar×nvar matrix. The previous
            // bordered-banded scheme would often fail on well-conditioned
            // Jacobians and silently fall back to a dense O(N^3) solve, which
            // dominated runtime. The matrix is strictly banded (half-bandwidth
            // = 3*fields - 1), so direct Gaussian elimination on the band is
            // the right primitive.
            const size_t half_bandwidth = 3U * active_field_count() - 1U;
            banded_ok = solve_compact_banded_system(J, rhs, delta, half_bandwidth);
        }
        if (!banded_ok)
        {
            if (banded_mode)
            {
                const size_t half_bandwidth = 3U * active_field_count() - 1U;
                std::vector<double> dense_J(nvar * nvar, 0.0);
                for (size_t col = 0; col < nvar; ++col)
                {
                    const size_t row_begin = (col > half_bandwidth) ? (col - half_bandwidth) : 0U;
                    const size_t row_end = std::min(nvar - 1U, col + half_bandwidth);
                    for (size_t row = row_begin; row <= row_end; ++row)
                    {
                        dense_J[row * nvar + col] = J[(half_bandwidth + row - col) * nvar + col];
                    }
                }

                if (!solve_linear_system(std::move(dense_J), rhs, delta))
                {
                    unpack_state(best_x);
                    last_residual_norm_ = best_norm;
                    return false;
                }
            }
            else if (!solve_linear_system(J, rhs, delta))
            {
                unpack_state(best_x);
                last_residual_norm_ = best_norm;
                return false;
            }
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

            if (banded_mode)
            {
                unpack_state_interleaved(trial);
            }
            else
            {
                unpack_state(trial);
            }
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
            if (banded_mode)
            {
                unpack_state_interleaved(best_x);
            }
            else
            {
                unpack_state(best_x);
            }
            last_residual_norm_ = best_norm;
            return false;
        }

        if (banded_mode)
        {
            unpack_state_interleaved(x);
        }
        else
        {
            unpack_state(x);
        }
    }

    const double final_norm = residual_norm(compute_residuals());
    if (std::isfinite(final_norm) && final_norm < best_norm)
    {
        best_norm = final_norm;
        best_x = x;
    }
    if (best_norm < p_.newton_tolerance)
    {
        if (banded_mode)
        {
            unpack_state_interleaved(best_x);
        }
        else
        {
            unpack_state(best_x);
        }
        last_residual_norm_ = best_norm;
        return true;
    }

    if (banded_mode)
    {
        unpack_state_interleaved(best_x);
    }
    else
    {
        unpack_state(best_x);
    }
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

    if (p_.use_suppression_factor)
        r_cutoff_ = radius_at_mass_fraction(p_.cutoff_percent / 100.0);

    if (p_.solver_mode != static_cast<int>(SolverMode::DenseLegacy) &&
        p_.solver_mode != static_cast<int>(SolverMode::BlockBanded))
    {
        if (!quiet)
        {
            std::cout << "[warn] Unknown solver_mode = " << p_.solver_mode
                      << ", falling back to dense legacy solver\n";
        }
        p_.solver_mode = static_cast<int>(SolverMode::DenseLegacy);
    }

    if (!quiet)
    {
        std::cout << "[INFO] solver_mode = "
                  << (p_.solver_mode == static_cast<int>(SolverMode::BlockBanded) ? "block_banded" : "dense_legacy")
                  << '\n';
    }

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
            if (p_.use_suppression_factor)
                r_cutoff_ = radius_at_mass_fraction(p_.cutoff_percent / 100.0);
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
            write_profile(p_.output_profile);
            if (p_.use_suppression_factor)
                r_cutoff_ = radius_at_mass_fraction(p_.cutoff_percent / 100.0);
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
    // Columns 1-11 are unchanged from the original relaxer output; 12-17 were appended for SBSE (FiniteTempBackground
    // reads them), and 18-19 are the particle-number integrands that the optimized upstream relaxer writes as 14-15.
    out << "# rbar r sigma eta alpha Xs omega Tlocal P_scalar P_thermal_suppressed P_thermal_raw"
        << " mu_local P_total s n Upsilon rho n_integrand Noether_integrand\n";

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
        const double suppression_factor = compute_suppression_factor(r);
        const ThermoConjugates tc = conjugates_analytic(Tlocal, mu_chem, sigma_sq, X, Xi, r);
        const double omega_i = state_.omega[static_cast<size_t>(i)];
        const double Xs_i = std::max(1e-12, state_.Xs[static_cast<size_t>(i)]);
        const double N_integrand = 4.0 * kPi * Xs_i * r * r * (tc.Upsilon * omega_i / alpha + tc.n);
        const double Noether_integrand = 4.0 * kPi * Xs_i * r * r * (omega_i * sigma_sq / alpha);
        
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
                                         alpha, state_.omega[static_cast<size_t>(i)], r);
        out << ' ' << th.mu << ' ' << th.P << ' ' << th.s << ' ' << th.n << ' ' << th.Upsilon << ' ' << th.rho
            << ' ' << N_integrand << ' ' << Noether_integrand << '\n';
    }
}

FiniteTempRelaxer::PointThermo FiniteTempRelaxer::thermo_at(double sigma, double eta, double alpha, double omega, double r) const
{
    PointThermo th;
    const double a = std::max(1e-12, alpha);
    th.T = local_temperature(a);
    th.mu = omega / a;
    th.X = eta * eta;
    th.Xi = -th.mu * th.mu;
    const double sigma_sq = sigma * sigma;

    th.P = pressure_Pthermal(th.T, th.mu, sigma_sq, th.X, th.Xi, r);
    const ThermoConjugates tc = conjugates_analytic(th.T, th.mu, sigma_sq, th.X, th.Xi, r);
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

double FiniteTempRelaxer::qball_total_energy() const
{
    const int n = p_.n_gridpoints;
    if (n < 2)
        return 0.0;

    const int last_interior = std::max(1, n - 2);
    auto integrand = [this](int i) -> double
    {
        const size_t j        = static_cast<size_t>(i);
        const double r        = areal_radius(state_.rbar[j]);
        const double alpha    = std::max(1e-12, state_.alpha[j]);
        const double sigma_sq = state_.sigma[j] * state_.sigma[j];
        const double X        = state_.eta[j] * state_.eta[j];
        const double omega    = state_.omega[j];
        const double mu_chem  = omega / alpha;
        const double Xi       = -mu_chem * mu_chem;
        const double Tlocal   = local_temperature(alpha);
        const double P        = pressure_Pthermal(Tlocal, mu_chem, sigma_sq, X, Xi, r);
        const ThermoConjugates tc = conjugates_analytic(Tlocal, mu_chem, sigma_sq, X, Xi, r);
        const double energy_density = tc.s * Tlocal + mu_chem * tc.n - P + tc.Upsilon * omega * omega;
        return 4.0 * kPi * r * r * energy_density;
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

double FiniteTempRelaxer::asymptotic_mass() const
{
    if (p_.qball_mode)
        return qball_total_energy();

    const int last_interior = std::max(0, p_.n_gridpoints - 2);
    const double r = areal_radius(state_.rbar[static_cast<size_t>(last_interior)]);
    return metric_mass(r, state_.Xs[static_cast<size_t>(last_interior)]);
}

double FiniteTempRelaxer::radius_at_mass_fraction(double mass_fraction) const
{
    const int last_interior = std::max(1, p_.n_gridpoints - 2);
    const double mf = std::min(1.0, std::max(0.0, mass_fraction));

    std::vector<double> rbar_vals;
    std::vector<double> masses;

    if (p_.qball_mode)
    {
        // Build cumulative energy integral using trapezoidal rule in r
        auto integrand = [this](int i) -> double
        {
            const size_t j        = static_cast<size_t>(i);
            const double r        = areal_radius(state_.rbar[j]);
            const double alpha    = std::max(1e-12, state_.alpha[j]);
            const double sigma_sq = state_.sigma[j] * state_.sigma[j];
            const double X        = state_.eta[j] * state_.eta[j];
            const double omega    = state_.omega[j];
            const double mu_chem  = omega / alpha;
            const double Xi       = -mu_chem * mu_chem;
            const double Tlocal   = local_temperature(alpha);
            const double P        = pressure_Pthermal(Tlocal, mu_chem, sigma_sq, X, Xi, r);
            const ThermoConjugates tc = conjugates_analytic(Tlocal, mu_chem, sigma_sq, X, Xi, r);
            const double energy_density = tc.s * Tlocal + mu_chem * tc.n - P + tc.Upsilon * omega * omega;
            return 4.0 * kPi * r * r * energy_density;
        };

        double cumulative = 0.0;
        double prev_r = areal_radius(state_.rbar[0]);
        double prev_f = integrand(0);
        rbar_vals.push_back(state_.rbar[0]);
        masses.push_back(0.0);
        for (int i = 1; i <= last_interior; ++i)
        {
            const size_t j = static_cast<size_t>(i);
            const double r = areal_radius(state_.rbar[j]);
            const double f = integrand(i);
            cumulative += 0.5 * (prev_f + f) * (r - prev_r);
            rbar_vals.push_back(state_.rbar[j]);
            masses.push_back(cumulative);
            prev_r = r;
            prev_f = f;
        }
    }
    else
    {
        for (int i = 0; i <= last_interior; ++i)
        {
            const double r = areal_radius(state_.rbar[static_cast<size_t>(i)]);
            const double m = metric_mass(r, state_.Xs[static_cast<size_t>(i)]);
            rbar_vals.push_back(state_.rbar[static_cast<size_t>(i)]);
            masses.push_back(m);
        }
    }

    const double mass_inf = masses.back();
    const double target_mass = mf * mass_inf;

    if (!std::isfinite(mass_inf) || mass_inf <= 0.0)
        return areal_radius(rbar_vals.back());

    if (rbar_vals.size() < 2)
        return areal_radius(rbar_vals.back());

    if (target_mass <= masses.front())
        return areal_radius(rbar_vals.front());

    if (target_mass >= masses.back())
        return areal_radius(rbar_vals.back());

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
                return areal_radius(rbar_vals[i + 1]);
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
        const ThermoConjugates tc = conjugates_analytic(Tlocal, mu_chem, sigma_sq, X, Xi, r);

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
    const double P = pressure_Pthermal(temp_local, mu_chem, sigma_sq, X, Xi, r);
    const ThermoConjugates tc = conjugates_analytic(temp_local, mu_chem, sigma_sq, X, Xi, r);

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