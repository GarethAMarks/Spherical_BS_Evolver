#include "seed_shooter.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace seedbs
{

SeedState operator+(const SeedState& s1, const SeedState& s2)
{
    return SeedState{s1.A + s2.A, s1.X + s2.X, s1.phi + s2.phi, s1.eta + s2.eta};
}

SeedState operator-(const SeedState& s1, const SeedState& s2)
{
    return SeedState{s1.A - s2.A, s1.X - s2.X, s1.phi - s2.phi, s1.eta - s2.eta};
}

SeedState operator*(double c, const SeedState& s)
{
    return SeedState{c * s.A, c * s.X, c * s.phi, c * s.eta};
}

double BosonStarSeedShooter::V(const double A) const
{
    if (!solitonic)
    {
        return mu * mu * A * A + lambda * std::pow(A, 4);
    }

    return mu * mu * A * A * std::pow((1.0 - 2.0 * std::pow(A / sigma_scale, 2)), 2) + lambda * std::pow(A, 4);
}

double BosonStarSeedShooter::dV(const double A) const
{
    if (!solitonic)
    {
        return mu * mu + 2.0 * A * A * lambda;
    }

    return mu * mu - 8.0 * mu * mu * std::pow(A / sigma_scale, 2)
        + 12.0 * mu * mu * std::pow(A / sigma_scale, 4)
        + 2.0 * A * A * lambda;
}

double BosonStarSeedShooter::ddV(const double A) const
{
    if (!solitonic)
    {
        return 2.0 * lambda;
    }

    return -8.0 * mu * mu * std::pow(1.0 / sigma_scale, 2)
        + 24.0 * mu * mu * A * A * std::pow(1.0 / sigma_scale, 4)
        + 2.0 * lambda;
}

SeedState BosonStarSeedShooter::state_RHS(const double radius, long double frequency, SeedState s, bool asymptotic_region) const
{
    if (qball_mode)
    {
        // Flat Minkowski: X = 1, alpha = 1 (phi = 0) throughout.
        // Metric fields do not evolve; only the scalar sector is integrated.
        const double eta_corr = (radius == 0.0) ? (1.0 / (D - 1.0)) : 1.0;
        const double T3       = (radius == 0.0) ? 0.0 : s.eta / radius;
        return SeedState{
            s.eta,
            0.0,
            0.0,
            static_cast<double>(eta_corr * (-(D - 2.0) * T3
                + s.A * (dV(s.A) - static_cast<double>(frequency * frequency))))};
    }

    const double epsilon = 1e-7;
    const double r = (radius == 0.0) ? epsilon : radius;

    double T1 = 0.5 * (D - 3.0) * (s.X * s.X - 1.0) / r;
    double T2 = s.eta * s.eta + frequency * frequency * s.A * s.A / std::exp(2.0 * s.phi);
    double T3 = s.eta / r;

    const double F1 = 4.0 * M_PI * G * r * s.X * s.X / (D - 2.0);

    if (radius == 0.0)
    {
        T1 = 0.0;
        T3 = 0.0;
    }

    const double dPhi = T1 + F1 * (T2 - V(s.A));
    const double eta_corr = (radius == 0.0) ? (1.0 / (D - 1.0)) : 1.0;

    if (asymptotic_region)
    {
        return SeedState{0.0, s.X * (F1 * (T2 + V(s.A)) - T1), dPhi, 0.0};
    }

    return SeedState{
        s.X * s.eta,
        s.X * (F1 * (T2 + V(s.A)) - T1),
        dPhi,
        static_cast<double>(eta_corr * (-(D - 2.0) * T3 - s.eta * dPhi + s.X * s.A * (dV(s.A) - frequency * frequency / std::exp(2.0 * s.phi))))};
}

SeedState BosonStarSeedShooter::state_expansion(const double radius, long double frequency) const
{
    const double& r = radius;
    const double A0 = A_central;
    const double alpha0 = alpha_central;
    const double fac_1 = frequency * frequency * A_central * A_central / (alpha_central * alpha_central);

    const double X2 = 8.0 * M_PI * fac_1 / ((D - 2.0) * (D - 1.0));
    const double alpha2 = alpha_central * (0.5 * (D - 3.0) * X2 + 4.0 * M_PI * fac_1) / (D - 2.0);

    const double A2 = (fac_1 - A_central * dV(A_central)) / (1.0 - D);

    const double X4 = (144 * M_PI * alpha0 * alpha0 * alpha0 * (X2 * r * r + 2.0 / 3.0) * V(A0)
        + 96 * A0 * A2 * dV(A0) * M_PI * r * r * alpha0 * alpha0 * alpha0
        + (((-9 * D * D + 45 * D - 54) * X2 * X2 + 96 * M_PI * A2 * A2) * r * r - 12 * X2 * (-1 + D) * (D - 2)) * alpha0 * alpha0 * alpha0
        + 96 * frequency * frequency * M_PI * ((3.0 / 2.0 * X2 * A0 + A2) * r * r + A0) * A0 * alpha0
        - 96 * A0 * A0 * M_PI * frequency * frequency * r * r * alpha2)
        / (D - 2) / (alpha0 * alpha0 * alpha0) / (r * r) / (1 + D);

    const double alpha4 = 0.25 * (-48 * (2 * (X2 * r * r + 1) * alpha0 + alpha2 * r * r) * M_PI * alpha0 * alpha0 * V(A0)
        - 96 * A2 * A0 * dV(A0) * M_PI * r * r * alpha0 * alpha0 * alpha0
        + ((3 * (D * D - 5 * D + 6) * X2 * X2 + 96 * M_PI * A2 * A2 + D * D * X4 - 5 * D * X4 + 6 * X4) * r * r + 12 * X2 * (D - 2) * (D - 3)) * alpha0 * alpha0 * alpha0
        + 6 * (-4 + (D - 3) * X2 * r * r) * alpha2 * (D - 2) * alpha0 * alpha0
        + 96 * A0 * M_PI * frequency * frequency * ((A0 * X2 + A2) * r * r + A0) * alpha0
        - 48 * A0 * A0 * M_PI * frequency * frequency * r * r * alpha2)
        / (alpha0 * alpha0) / (D - 2) / (r * r);

    return SeedState{
        A_central + 0.5 * A2 * r * r,
        1.0 + 0.5 * X2 * r * r + X4 * std::pow(r, 4.0) / 24.0,
        std::log(alpha_central + 0.5 * alpha2 * r * r + alpha4 * std::pow(r, 4.0) / 24.0),
        A2 * r};
}

void BosonStarSeedShooter::rk4_solve(long double freq)
{
    state = {SeedState{A_central, 1.0, std::log(alpha_central), 0.0}};
    radius_array = {0.0};

    blowup_point = n_gridpoints;
    state.resize(static_cast<size_t>(n_gridpoints));
    radius_array.resize(static_cast<size_t>(n_gridpoints));

    const double dr = R / (n_gridpoints - 1);

    SeedState s1, s2, s3, s4;
    for (int j = 0; j < n_gridpoints - 1; ++j)
    {
        const double r = j * dr;

        if (D > 100.0 && j < 1)
        {
            s1 = state_RHS(r, freq, state[static_cast<size_t>(j)], false);

            SeedState t1 = state[static_cast<size_t>(j)] + 0.5 * dr * s1;
            SeedState r2_state = state_expansion(r + 0.5 * dr, freq);
            t1.X = r2_state.X;
            t1.phi = r2_state.phi;

            s2 = state_RHS(r + dr / 2.0, freq, t1, false);

            SeedState t2 = state[static_cast<size_t>(j)] + 0.5 * dr * s2;
            t2.X = r2_state.X;
            t2.phi = r2_state.phi;

            s3 = state_RHS(r + dr / 2.0, freq, t2, false);

            SeedState t3 = state[static_cast<size_t>(j)] + 0.5 * dr * s3;
            SeedState r_state = state_expansion(r + dr, freq);
            t3.X = r_state.X;
            t3.phi = r_state.phi;

            s4 = state_RHS(r + dr, freq, t3, false);
        }
        else
        {
            s1 = state_RHS(r, freq, state[static_cast<size_t>(j)], false);
            s2 = state_RHS(r + dr / 2.0, freq, state[static_cast<size_t>(j)] + 0.5 * dr * s1, false);
            s3 = state_RHS(r + dr / 2.0, freq, state[static_cast<size_t>(j)] + 0.5 * dr * s2, false);
            s4 = state_RHS(r + dr, freq, state[static_cast<size_t>(j)] + dr * s3, false);
        }

        state[static_cast<size_t>(j + 1)] = state[static_cast<size_t>(j)] + (dr / 6.0) * (s1 + 2 * s2 + 2 * s3 + s4);
        radius_array[static_cast<size_t>(j + 1)] = (j + 1) * dr;

        if (D > 4.0 && j < 1)
        {
            SeedState small_state = state_expansion(r + dr, freq);
            state[static_cast<size_t>(j + 1)].X = small_state.X;
            state[static_cast<size_t>(j + 1)].phi = small_state.phi;
        }

        if (std::isnan(state[static_cast<size_t>(j)].A)
            || std::isnan(state[static_cast<size_t>(j)].X)
            || std::isnan(state[static_cast<size_t>(j)].phi)
            || std::isnan(state[static_cast<size_t>(j)].eta))
        {
            blowup_point = j;
            return;
        }
    }
}

double BosonStarSeedShooter::m(int j) const
{
    if (j >= n_gridpoints || j < 0)
    {
        std::cerr << "ERROR: invalid index passed to m(r)" << std::endl;
        std::exit(1);
    }

    return (D - 2.0) * std::pow(M_PI, 0.5 * (D - 3.0)) * std::pow(radius_array[static_cast<size_t>(j)], D - 3.0)
        * (1.0 - (1.0 / (state[static_cast<size_t>(j)].X * state[static_cast<size_t>(j)].X)))
        / (8.0 * tgamma(0.5 * (D - 1.0)));
}

int BosonStarSeedShooter::count_zero_crossings() const
{
    int zero_crossings = 0;
    for (int j = 1; j < blowup_point - 1; ++j)
    {
        const double Aj = state[static_cast<size_t>(j)].A;
        const double Ajm = state[static_cast<size_t>(j - 1)].A;
        if ((Aj == 0.0) || (Aj > 0.0 && Ajm < 0.0) || (Aj < 0.0 && Ajm > 0.0))
        {
            ++zero_crossings;
        }
    }
    return zero_crossings;
}

long double BosonStarSeedShooter::find_frequency(bool quiet)
{
    long double lower_guess = frequency_guess;
    const long double epsilon = freq_epsilon;

    rk4_solve(lower_guess);
    if (count_zero_crossings() <= eigen && !quiet)
    {
        std::cerr << "WARNING: initial guess of " << frequency_guess << " too small for A_central = " << A_central << std::endl;
    }

    while (count_zero_crossings() > eigen && lower_guess > epsilon)
    {
        lower_guess /= 2.0;
        rk4_solve(lower_guess);
    }

    if (lower_guess <= epsilon)
    {
        std::cerr << "ERROR: could not find suitable lower frequency bound for A_central = " << A_central << std::endl;
    }

    long double upper_guess = lower_guess * 2.0;
    long double midpoint = 0.5 * (upper_guess + lower_guess);
    int num_solves = 0;

    while ((upper_guess - lower_guess) > epsilon && num_solves < 100)
    {
        midpoint = 0.5 * (upper_guess + lower_guess);
        rk4_solve(midpoint);

        if (count_zero_crossings() > eigen)
        {
            upper_guess = midpoint;
        }
        else
        {
            lower_guess = midpoint;
        }

        ++num_solves;
    }

    omega = lower_guess;
    rk4_solve(omega);

    if (count_zero_crossings() != eigen)
    {
        std::cerr << "WARNING: zero crossing count may be incorrect for A_central = " << A_central << ", found " << count_zero_crossings() << std::endl;
    }

    if (!quiet)
    {
        std::cout << "\nFound solution in eigenstate " << eigen << " with frequency " << static_cast<double>(omega) << std::endl;
    }

    return lower_guess;
}

int BosonStarSeedShooter::find_last_minimum() const
{
    for (int j = blowup_point - 2; j > 1; --j)
    {
        const double Aj = std::abs(state[static_cast<size_t>(j)].A);
        const double Ajm = std::abs(state[static_cast<size_t>(j - 1)].A);
        const double Ajp = std::abs(state[static_cast<size_t>(j + 1)].A);
        if (Aj <= Ajm && Aj <= Ajp && state[static_cast<size_t>(j)].A < 0.05 * A_central)
        {
            return j;
        }
    }

    if (blowup_point == n_gridpoints
        && std::abs(state[static_cast<size_t>(n_gridpoints - 1)].A) <= std::abs(state[static_cast<size_t>(n_gridpoints - 2)].A))
    {
        return n_gridpoints - 1;
    }

    return -1;
}

bool BosonStarSeedShooter::fill_asymptotic(bool quiet)
{
    const int min_index = find_last_minimum();
    if (min_index == -1 && !quiet)
    {
        std::cerr << "ERROR: No local minimum in |A| found for A_central = " << A_central << std::endl;
    }
    if (min_index == -1)
    {
        return false;
    }

    const int j_match = static_cast<int>(std::round(r_match_fac * min_index));

    const double dr = R / (n_gridpoints - 1);
    const double r_match = dr * j_match;

    const double phi_match = state[static_cast<size_t>(j_match)].phi;
    const double A_match = state[static_cast<size_t>(j_match)].A;

    const double kappa = std::sqrt(std::max(1e-14, mu * mu - std::pow(static_cast<double>(omega) / std::exp(phi_match), 2)));
    const double A_factor = A_match * std::exp(kappa * r_match) * std::pow(r_match, 0.5 * (D - 2.0));

    SeedState s1, s2, s3, s4;
    for (int j = j_match; j < n_gridpoints - 1; ++j)
    {
        const long double r = j * dr;

        s1 = state_RHS(r, omega, state[static_cast<size_t>(j)], true);
        s2 = state_RHS(r + dr / 2.0, omega, state[static_cast<size_t>(j)] + 0.5 * dr * s1, true);
        s3 = state_RHS(r + dr / 2.0, omega, state[static_cast<size_t>(j)] + 0.5 * dr * s2, true);
        s4 = state_RHS(r + dr, omega, state[static_cast<size_t>(j)] + dr * s3, true);

        state[static_cast<size_t>(j + 1)] = state[static_cast<size_t>(j)] + (dr / 6.0) * (s1 + 2 * s2 + 2 * s3 + s4);
        radius_array[static_cast<size_t>(j + 1)] = (j + 1) * dr;

        const double rp = radius_array[static_cast<size_t>(j + 1)];
        state[static_cast<size_t>(j + 1)].A = A_factor * std::exp(-kappa * rp) / std::pow(rp, 0.5 * (D - 2.0));
        state[static_cast<size_t>(j + 1)].eta = -(0.5 * (D - 2.0) / rp + kappa) * state[static_cast<size_t>(j + 1)].A;

        if (std::isnan(state[static_cast<size_t>(j)].A)
            || std::isnan(state[static_cast<size_t>(j)].X)
            || std::isnan(state[static_cast<size_t>(j)].phi)
            || std::isnan(state[static_cast<size_t>(j)].eta))
        {
            std::cerr << "WARNING: During asymptotic phase, state values became nan on step " << j << std::endl;
            return false;
        }
    }

    const double phi_shift = -std::log(state[static_cast<size_t>(n_gridpoints - 1)].X) - state[static_cast<size_t>(n_gridpoints - 1)].phi;
    rescale_lapse(phi_shift);

    return true;
}

void BosonStarSeedShooter::rescale_lapse(double phi_shift)
{
    for (SeedState& s : state)
    {
        s.phi += phi_shift;
    }

    omega *= std::exp(phi_shift);
}

bool BosonStarSeedShooter::solve(bool quiet)
{
    find_frequency(quiet);
    return fill_asymptotic(quiet);
}

void BosonStarSeedShooter::write_field(const std::string& filename) const
{
    std::ofstream data_file{filename};
    if (!data_file)
    {
        std::cerr << "BS seed file could not be opened for writing!\n";
        return;
    }

    data_file << "# r sigma Xs alpha eta omega m\n";
    data_file << std::scientific << std::setprecision(12);

    const int last = std::max(0, n_gridpoints - 1);
    for (int j = 0; j < last; ++j)
    {
        const double r = radius_array[static_cast<size_t>(j)];
        const double sigma = state[static_cast<size_t>(j)].A;
        const double Xs = state[static_cast<size_t>(j)].X;
        const double alpha = std::exp(state[static_cast<size_t>(j)].phi);
        const double eta = state[static_cast<size_t>(j)].eta;
        const double om = static_cast<double>(omega);

        data_file << r << ' '
                  << sigma << ' '
                  << Xs << ' '
                  << alpha << ' '
                  << eta << ' '
                  << om << ' '
                  << m(j) << '\n';
    }
}

} // namespace seedbs
