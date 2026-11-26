#include "FiniteTempMatter.h"
#include "EvolutionVariables.h"
#include "mathutils.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using std::max;
using std::min;

// ---------------------------------------------------------------------------------------------------------------------
// FTM arithmetic
FTM operator+(const FTM& a, const FTM& b) { return FTM{a.K_sigma + b.K_sigma, a.N_S + b.N_S, a.N_N + b.N_N, a.Sigma + b.Sigma, a.P_z + b.P_z}; }
FTM operator-(const FTM& a, const FTM& b) { return FTM{a.K_sigma - b.K_sigma, a.N_S - b.N_S, a.N_N - b.N_N, a.Sigma - b.Sigma, a.P_z - b.P_z}; }
FTM operator*(const FTM& a, const FTM& b) { return FTM{a.K_sigma * b.K_sigma, a.N_S * b.N_S, a.N_N * b.N_N, a.Sigma * b.Sigma, a.P_z * b.P_z}; }
FTM operator*(double c, const FTM& a) { return FTM{c * a.K_sigma, c * a.N_S, c * a.N_N, c * a.Sigma, c * a.P_z}; }

// ---------------------------------------------------------------------------------------------------------------------
// Equation of state (same function as FiniteTempRelaxer::pressure_Pthermal with use_suppression_factor = 0, but with
// analytic derivatives instead of finite differences)
double FiniteTempEoS::V(double s2) const
{
    s2 = max(0., s2);
    if (!solitonic) return m * m * s2 + lambda * s2 * s2;
    const double sc2 = max(1e-12, sigma_scale * sigma_scale);
    const double sh = 1. - 2. * s2 / sc2;
    return m * m * s2 * sh * sh + lambda * s2 * s2;
}
double FiniteTempEoS::dV(double s2) const
{
    if (!solitonic) return m * m + 2. * lambda * s2;
    const double x = s2 / max(1e-24, sigma_scale * sigma_scale);
    return m * m * (1. - 8. * x + 12. * x * x) + 2. * lambda * s2;
}
double FiniteTempEoS::ddV(double s2) const
{
    if (!solitonic) return 2. * lambda;
    const double sc2 = max(1e-12, sigma_scale * sigma_scale);
    return -8. * m * m / sc2 + 24. * m * m * s2 / (sc2 * sc2) + 2. * lambda;
}
double FiniteTempEoS::dddV(double /*s2*/) const
{
    if (!solitonic) return 0.;
    const double sc2 = max(1e-12, sigma_scale * sigma_scale);
    return 24. * m * m / (sc2 * sc2);
}

FiniteTempEoS::Thermo FiniteTempEoS::thermal(double T, double mu, double sigma_sq, double Xi) const
{
    Thermo th;
    if (!(T >= 1e-12)) return th;

    const double s2 = max(0., sigma_sq);
    const double Xs = (std::abs(Xi) < 1e-12) ? ((Xi >= 0.) ? 1e-12 : -1e-12) : Xi;

    // phonon speed (same regularisation as the relaxer)
    const double vpp = ddV(s2);
    const double a = s2 * vpp;
    const double a_p = vpp + s2 * dddV(s2); // da / d(sigma^2)
    const double num = a + max(0., cph_epsilon);
    double den = a - 2. * Xs;
    if (std::abs(den) < 1e-12) den = (den >= 0.) ? 1e-12 : -1e-12;
    double c2 = num / den;
    double dc2_dXi = 2. * num / (den * den);
    double dc2_ds2 = a_p * (den - num) / (den * den);
    if (c2 < 1e-12) { c2 = 1e-12; dc2_dXi = 0.; dc2_ds2 = 0.; }
    const double c = std::sqrt(c2);

    const double B0 = 1. + (1. - c2) * mu * mu / Xs;
    const double B_mu = 2. * (1. - c2) * mu / Xs, B_c2 = -mu * mu / Xs, B_Xs = -(1. - c2) * mu * mu / (Xs * Xs);
    // optional convex floor on B (B_floor = beta > 0): B_eff = B for B >= B1 = beta c^2, B_eff = B1^2/(2 B1 - B) below.
    // Equilibrium has B = c^2, so the background is untouched for beta < 1; the floor removes the pole of P_ph at B = 0
    // (and the non-convex B < 0 branch). B_eff is C^1, monotone and keeps P_ph convex in mu.
    double B = B0, Be_B = 1., Be_c2 = 0.;
    const double B1 = B_floor * c2;
    if (B_floor > 0. && B0 < B1)
    {
        const double d = 2. * B1 - B0; // > B1 > 0
        B = B1 * B1 / d; Be_B = B * B / (B1 * B1); Be_c2 = B_floor * 2. * B1 * (B1 - B0) / (d * d);
    }
    if (std::abs(B) < 1e-12) B = (B >= 0.) ? 1e-12 : -1e-12;

    const double kph = M_PI * M_PI / 90.;
    const double T4 = T * T * T * T;
    const double Pph = kph * c * T4 / (B * B);
    // partial derivatives of P_ph
    const double Pph_T = 4. * Pph / T;
    const double Pph_mu = -2. * Pph / B * Be_B * B_mu;
    const double Pph_c2 = Pph * (0.5 / c2 - 2. / B * (Be_B * B_c2 + Be_c2));
    const double Pph_Xi = -2. * Pph / B * Be_B * B_Xs + Pph_c2 * dc2_dXi;
    const double Pph_s2 = Pph_c2 * dc2_ds2;

    // free Boltzmann gas
    const double lnPg = 1.5 * std::log(m / (2. * M_PI)) + 2.5 * std::log(T) - (m - mu) / T;
    const double Pg_T_over_Pg = 2.5 / T + (m - mu) / (T * T);
    const double Pg_mu_over_Pg = 1. / T;

    if (crossover_mode == 1)
    {
        // convex blend P_th = w(s2) P_ph + g(s2) P_gas
        const double kh = 0.5 * crossover_exponent;
        if (!(s2 > 0.))
        {
            const double Pg = std::exp(lnPg);
            th.P_th = Pg; th.s = Pg * Pg_T_over_Pg; th.n = Pg * Pg_mu_over_Pg;
            return th;
        }
        const double x_ph = std::pow(s2 / (sigma_phonon * sigma_phonon), kh);
        const double x_g = std::pow(s2 / (sigma_gas * sigma_gas), kh);
        const double w = -std::expm1(-x_ph);
        const double w_s2 = std::exp(-x_ph) * x_ph * kh / s2;
        const double Pg = std::exp(lnPg - x_g);  // g P_gas
        const double Pg_s2 = -Pg * x_g * kh / s2; // d(g P_gas)/d(s2)
        const bool ph = (w > 0.) && std::isfinite(Pph);
        th.P_th = (ph ? w * Pph : 0.) + Pg;
        th.s = (ph ? w * Pph_T : 0.) + Pg * Pg_T_over_Pg;
        th.n = (ph ? w * Pph_mu : 0.) + Pg * Pg_mu_over_Pg;
        th.dF = 2. * ((ph ? w * Pph_s2 + w_s2 * Pph : 0.) + Pg_s2);
        th.dUps = ph ? -2. * w * Pph_Xi : 0.;
        if (!std::isfinite(th.s) || !std::isfinite(th.n) || !std::isfinite(th.dF) || !std::isfinite(th.dUps)) { th = Thermo{}; th.P_th = Pg; }
        return th;
    }

    // q-norm interpolation, in logs: P_th = (P_ph^-q + P_g^-q)^(-1/q)
    const double qq = max(1e-12, q);
    if (!(Pph > 0.))
    {
        // degenerate phonon branch (B^2 -> inf etc.): fall back to the gas
        const double Pg = std::exp(lnPg);
        th.P_th = Pg; th.s = Pg * Pg_T_over_Pg; th.n = Pg * Pg_mu_over_Pg;
        return th;
    }
    const double lnPph = std::log(Pph);
    const double u1 = -qq * lnPph, u2 = -qq * lnPg;
    const double umax = max(u1, u2);
    const double lse = umax + std::log1p(std::exp(min(u1, u2) - umax));
    const double lnPth = -lse / qq;
    const double Pth = std::exp(lnPth);
    const double w_ph = std::exp((1. + qq) * (lnPth - lnPph)); // (P_th / P_ph)^{1+q}
    const double w_g = std::exp((1. + qq) * (lnPth - lnPg));
    const double Pg = std::exp(lnPg);

    th.P_th = Pth;
    th.s = w_ph * Pph_T + w_g * Pg * Pg_T_over_Pg;
    th.n = w_ph * Pph_mu + w_g * Pg * Pg_mu_over_Pg;
    th.dF = 2. * w_ph * Pph_s2;
    th.dUps = -2. * w_ph * Pph_Xi;
    if (!std::isfinite(th.s) || !std::isfinite(th.n)) { th = Thermo{}; th.P_th = Pth; }
    return th;
}

// ---------------------------------------------------------------------------------------------------------------------
// helpers

double FiniteTempMatter::deriv(const BSSNSlice* slice, const std::vector<double>& f, int j, int order, bool odd) const
{
    BSSNSlice::DerivStencil stn;
    slice->build_stencil(j, stn);
    double F[5];
    for (int k = 0; k < 5; k++) F[k] = f[stn.J[k]];
    if (odd && stn.flip0) F[0] = -F[0];
    if (odd && stn.flip1) F[1] = -F[1];
    const double h = stn.h;
    if (stn.outer_pos != 0)
    {
        if (order == 1)
            return (stn.outer_pos == 5) ? p5_stencil(h, F[0], F[1], F[2], F[3], F[4]) : p4_stencil(h, F[0], F[1], F[2], F[3], F[4]);
        return (stn.outer_pos == 5) ? (11. * F[0] - 56. * F[1] + 114. * F[2] - 104. * F[3] + 35. * F[4]) / (12. * h * h)
                                    : (-F[0] + 4. * F[1] + 6. * F[2] - 20. * F[3] + 11. * F[4]) / (12. * h * h);
    }
    return fivePointDeriv(h, order, F[0], F[1], F[2], F[3], F[4]);
}

namespace
{
bool solve4(double A[4][4], double b[4], double x[4])
{
    int idx[4] = {0, 1, 2, 3};
    double M[4][5];
    for (int i = 0; i < 4; i++) { for (int k = 0; k < 4; k++) M[i][k] = A[i][k]; M[i][4] = b[i]; }
    for (int c = 0; c < 4; c++)
    {
        int p = c;
        for (int r = c + 1; r < 4; r++) if (std::abs(M[r][c]) > std::abs(M[p][c])) p = r;
        if (!(std::abs(M[p][c]) > 1e-300)) return false;
        if (p != c) for (int k = 0; k < 5; k++) std::swap(M[p][k], M[c][k]);
        for (int r = c + 1; r < 4; r++)
        {
            const double f = M[r][c] / M[c][c];
            for (int k = c; k < 5; k++) M[r][k] -= f * M[c][k];
        }
    }
    for (int i = 3; i >= 0; i--)
    {
        double s = M[i][4];
        for (int k = i + 1; k < 4; k++) s -= M[i][k] * x[k];
        x[i] = s / M[i][i];
    }
    (void)idx;
    return true;
}
} // namespace

// Primitive recovery: find (T, mu, v, Pi_psi) with W n = N_N, W s = Sigma, w W^2 gamma_zz v = P_z, -Upsilon Pi = N_S.
bool FiniteTempMatter::recover(const FTM& U, double floor, const BSSNState& b, double sigma, double dz_sigma, double dz_psi,
                               double K_sigma, Point& p, double& res) const
{
    (void)dz_sigma; (void)K_sigma; // the thermal EoS has no X dependence (delta Y = 0)
    const double chi = max(st->min_chi, b.chi);
    const double g_zz = b.h_zz / chi;
    const double g_ZZ = chi / b.h_zz;
    const double s2 = sigma * sigma;
    const double grad_psi2 = g_ZZ * dz_psi * dz_psi;

    const bool no_fluid = !(U.N_N > floor) || !(U.Sigma > floor);
    if (no_fluid)
    {
        // T = 0: Upsilon = sigma^2 exactly
        p.T = 0.; p.v = 0.;
        if (s2 > 1e-300) p.Pi = -U.N_S / s2;
        p.mu = -p.Pi;
        res = 0.;
        return true;
    }

    // unknowns y = (ln T, mu, sqrt(g_zz) v, Pi)
    double y[4] = {std::log(max(1e-300, p.T > 0. ? p.T : 1e-3)), p.mu, std::sqrt(g_zz) * p.v, p.Pi};
    const double scale_S = max(std::abs(U.N_S), 1e-300);

    auto residual = [&](const double* yy, double* r) -> bool
    {
        const double T = std::exp(yy[0]), mu = yy[1], u = yy[2], Pi = yy[3];
        if (!(u * u < 1.)) return false;
        const double W = 1. / std::sqrt(1. - u * u);
        const double Xi = -Pi * Pi + grad_psi2;
        const FiniteTempEoS::Thermo th = eos.thermal(T, mu, s2, Xi);
        const double w = T * th.s + mu * th.n;
        const double Ups = s2 + th.dUps;
        r[0] = W * th.n / U.N_N - 1.;
        r[1] = W * th.s / U.Sigma - 1.;
        // P_z = w W^2 g_zz v = w W^2 sqrt(g_zz) u, scaled by the enthalpy density W w ~ T Sigma + mu N_N
        r[2] = (w * W * W * std::sqrt(g_zz) * u - U.P_z) / (std::sqrt(g_zz) * max(1e-300, T * U.Sigma + mu * U.N_N));
        r[3] = (-Ups * Pi - U.N_S) / scale_S;
        return std::isfinite(r[0]) && std::isfinite(r[1]) && std::isfinite(r[2]) && std::isfinite(r[3]);
    };
    auto norm = [](const double* r) { return max(max(std::abs(r[0]), std::abs(r[1])), max(std::abs(r[2]), std::abs(r[3]))); };

    double r[4];
    if (!residual(y, r))
    {
        // reset to the T -> 0 guess
        y[0] = std::log(1e-3); y[2] = 0.; if (s2 > 1e-300) y[3] = -U.N_S / s2; y[1] = -y[3];
        if (!residual(y, r)) { res = 1e300; return false; }
    }
    double nr = norm(r);
    int it = 0;
    for (; it < c2p_max_iter && nr > c2p_tol; it++)
    {
        double J[4][4];
        for (int k = 0; k < 4; k++)
        {
            double yp[4] = {y[0], y[1], y[2], y[3]};
            const double hk = 1e-7 * max(1., std::abs(y[k]));
            yp[k] += hk;
            double rp[4];
            if (!residual(yp, rp))
            {
                yp[k] = y[k] - hk;
                if (!residual(yp, rp)) { res = nr; return false; }
                for (int i = 0; i < 4; i++) J[i][k] = (r[i] - rp[i]) / hk;
            }
            else
                for (int i = 0; i < 4; i++) J[i][k] = (rp[i] - r[i]) / hk;
        }
        double rhs[4] = {-r[0], -r[1], -r[2], -r[3]}, dy[4];
        if (!solve4(J, rhs, dy)) { res = nr; return false; }
        // limit the ln T step and backtrack on the residual norm
        const double lim = 1.0;
        double fac = 1.;
        if (std::abs(dy[0]) > lim) fac = lim / std::abs(dy[0]);
        bool accepted = false;
        for (int ls = 0; ls < 30; ls++)
        {
            double yt[4], rt[4];
            for (int k = 0; k < 4; k++) yt[k] = y[k] + fac * dy[k];
            if (residual(yt, rt) && norm(rt) < nr)
            {
                for (int k = 0; k < 4; k++) { y[k] = yt[k]; r[k] = rt[k]; }
                nr = norm(r);
                accepted = true;
                break;
            }
            fac *= 0.5;
        }
        if (!accepted) break;
    }

    p.T = std::exp(y[0]); p.mu = y[1]; p.v = y[2] / std::sqrt(g_zz); p.Pi = y[3];
    res = nr;
    return nr <= 1e-8;
}

// ---------------------------------------------------------------------------------------------------------------------

void FiniteTempMatter::prepare(BSSNSlice* slice)
{
    const int N = (int)slice->states2.size();
    pts.resize(N);
    f_sigma.resize(N); f_T.resize(N); f_mu.resize(N);
    coef_S.resize(N); flux_S.resize(N); atm.assign(N, 0); bath.assign(N, 0); atm_w.assign(N, 0.); flux_N.resize(N); flux_Sig.resize(N); flux_P.resize(N);
    if ((int)guess.size() != N) guess.assign(N, Point{});

    for (int j = 0; j < N; j++)
    {
        const CSF& c = slice->states2[j].csf;
        f_sigma[j] = std::sqrt(c.phi_re * c.phi_re + c.phi_im * c.phi_im);
    }

    c2p_failures = 0;
    theta_capped = 0;
    max_c2p_res = 0.;
    fluid_atm.assign(N, 0);
    double max_NN = 0.;
    for (int j = 0; j < N; j++)
        if (!(f_sigma[j] < fluid_sigma_min)) max_NN = max(max_NN, slice->states2[j].ftm.N_N); // evolved fluid only
    const double floor = max(n_floor, fluid_atm_frac * max_NN);
    for (int j = 0; j < N; j++)
    {
        State& s = slice->states2[j];
        const BSSNState& b = s.bssn;
        const double chi = max(st->min_chi, b.chi);
        const double g_zz = b.h_zz / chi, g_ZZ = chi / b.h_zz;
        const double sig = f_sigma[j];
        const double s2 = sig * sig;
        const double dre = slice->d_z(v_phi_re, j, 1), dim = slice->d_z(v_phi_im, j, 1);
        const double cross = s.csf.phi_re * dim - s.csf.phi_im * dre; // sigma^2 d_z psi
        const double dz_psi = (s2 > 1e-300) ? cross / s2 : 0.;
        const double dz_sigma = (sig > 1e-300) ? (s.csf.phi_re * dre + s.csf.phi_im * dim) / sig : 0.;

        const double alpha = b.alpha;
        atm[j] = !(sig >= sigma_atm);
        if (atm[j])
        {
            // CSF owns (phi, K_phi) here: slave K_sigma and N_S to it, normal component absent
            const double Kre = s.csf.K_phi_re, Kim = s.csf.K_phi_im;
            s.ftm.K_sigma = (sig > 1e-300) ? (s.csf.phi_re * Kre + s.csf.phi_im * Kim) / sig : 0.;
            s.ftm.N_S = 2. * (s.csf.phi_re * Kim - s.csf.phi_im * Kre);
            Point p{};
            p.sigma = sig; p.dz_sigma = dz_sigma; p.dz_psi = dz_psi; p.cross = cross;
            p.Pi = (s2 > 1e-300) ? -s.ftm.N_S / s2 : 0.;
            p.mu = -p.Pi; p.W = 1.; p.Ups = s2;
            p.Xi = -p.Pi * p.Pi + g_ZZ * dz_psi * dz_psi;
            pts[j] = p;
            fluid_atm[j] = 1;
            s.ftm.P_z = 0.;
            // only enter s T' + n mu' next to the fluid edge: continue them with the Tolman/Klein laws (T, mu ~ 1 / alpha)
            // so that the edge points stay in hydrostatic balance
            const double ar = (j > 0) ? slice->states2[j - 1].bssn.alpha / alpha : 1.;
            f_T[j] = (j > 0) ? f_T[j - 1] * ar : 0.;
            f_mu[j] = (j > 0) ? f_mu[j - 1] * ar : p.mu;
            coef_S[j] = alpha * g_ZZ;
            flux_S[j] = coef_S[j] * cross;
            flux_N[j] = flux_Sig[j] = flux_P[j] = 0.;
            continue;
        }

        Point p = guess[j];
        double res = 0.;
        bath[j] = (sig < fluid_sigma_min) && j < (int)T_alpha0.size() && T_alpha0[j] > 0.;
        fluid_atm[j] = !(s.ftm.N_N > floor) || !(s.ftm.Sigma > floor) || bath[j];
        // smooth taper of the normal-fluid evolution between floor and fluid_atm_taper * floor (1 = hard switch)
        if (fluid_atm[j]) atm_w[j] = 0.;
        else if (fluid_atm_taper > 1.)
        {
            const double x = min(1., max(0., (s.ftm.N_N / floor - 1.) / (fluid_atm_taper - 1.)));
            atm_w[j] = x * x * (3. - 2. * x);
        }
        else atm_w[j] = 1.;
        if (fluid_atm[j]) s.ftm.P_z = 0.;
        if (bath[j])
        {
            // 1D Newton for Pi at T = T_bath, mu = mu_bath, v = 0: -Upsilon(T, mu, Xi(Pi)) Pi = N_S
            const double Tb = T_alpha0[j] / b.alpha, mub = mu_alpha0[j] / b.alpha;
            const double gp2 = g_ZZ * dz_psi * dz_psi;
            auto fres = [&](double Pi) { return -(s2 + eos.thermal(Tb, mub, s2, -Pi * Pi + gp2).dUps) * Pi - s.ftm.N_S; };
            double Pi = (p.Pi < 0.) ? p.Pi : ((s2 > 1e-300) ? -s.ftm.N_S / s2 : -1.);
            bool ok = false;
            for (int it = 0; it < 40; it++)
            {
                const double f0 = fres(Pi);
                if (std::abs(f0) <= 1e-13 * max(1e-300, std::abs(s.ftm.N_S))) { ok = true; break; }
                const double h = 1e-7 * max(1., std::abs(Pi));
                const double df = (fres(Pi + h) - fres(Pi - h)) / (2. * h);
                if (!(std::abs(df) > 0.)) break;
                double step = -f0 / df;
                if (std::abs(step) > 0.2 * std::abs(Pi)) step = std::copysign(0.2 * std::abs(Pi), step);
                Pi += step;
            }
            if (!ok) c2p_failures++;
            p.T = Tb; p.mu = mub; p.Pi = Pi; p.v = 0.;
            res = 0.;
        }
        else if (!recover(s.ftm, floor, b, sig, dz_sigma, dz_psi, s.ftm.K_sigma, p, res))
            c2p_failures++;
        if (std::isfinite(res)) max_c2p_res = max(max_c2p_res, res);
        guess[j] = p;

        p.sigma = sig; p.dz_sigma = dz_sigma; p.dz_psi = dz_psi; p.cross = cross;
        p.Xi = -p.Pi * p.Pi + g_ZZ * dz_psi * dz_psi;
        const FiniteTempEoS::Thermo th = eos.thermal(p.T, p.mu, s2, p.Xi);
        const double u2 = min(1. - 1e-15, g_zz * p.v * p.v);
        p.W = 1. / std::sqrt(1. - u2);
        p.P_th = th.P_th; p.s = th.s; p.n = th.n; p.dF = th.dF; p.dUps = th.dUps;
        p.w = p.T * p.s + p.mu * p.n;
        p.Ups = s2 + th.dUps;
        if (bath[j]) { s.ftm.N_N = p.n; s.ftm.Sigma = p.s; }
        if (theta1 > 0. && !fluid_atm[j] && p.T > 0.)
        {
            p.J = p.W * (p.Pi + p.v * dz_psi) + p.mu;
            const double hmu = 1e-6 * max(1e-3, std::abs(p.mu));
            const double n_mu = (eos.thermal(p.T, p.mu + hmu, s2, p.Xi).n - eos.thermal(p.T, p.mu - hmu, s2, p.Xi).n) / (2. * hmu);
            const double rate = b.alpha * (1. / max(1e-300, p.W * n_mu) + 1. / max(1e-300, p.Ups));
            p.theta = theta_implicit ? 0. : theta1; // implicit: applied by relax_conversion after each step instead
            if (!theta_implicit && dt > 0. && theta1 * rate * dt > theta_cfl) { p.theta = theta_cfl / (rate * dt); theta_capped++; }
            if (!(n_mu > 0.) || !(p.Ups > 0.)) p.theta = 0.;
        }
        pts[j] = p;

        // derived SBSE momenta: K_phi = (phi / sigma) K_sigma - (i / 2) phi Pi_psi
        const double ks = (sig > 1e-300) ? s.ftm.K_sigma / sig : 0.;
        s.csf.K_phi_re = s.csf.phi_re * ks + 0.5 * s.csf.phi_im * p.Pi;
        s.csf.K_phi_im = s.csf.phi_im * ks - 0.5 * s.csf.phi_re * p.Pi;

        const double ar = (j > 0) ? slice->states2[j - 1].bssn.alpha / alpha : 1.;
        const bool cont = fluid_atm[j] && !bath[j] && j > 0 && !(p.T > 0.);
        f_T[j] = cont ? f_T[j - 1] * ar : p.T;
        f_mu[j] = cont ? f_mu[j - 1] * ar : p.mu;
        // alpha Upsilon gamma^zz d_z psi, with Upsilon d_z psi = (Upsilon / sigma^2) (sigma^2 d_z psi)
        coef_S[j] = (s2 > 1e-300) ? alpha * g_ZZ * (p.Ups / s2) : alpha * g_ZZ;
        flux_S[j] = coef_S[j] * cross;
        flux_N[j] = alpha * s.ftm.N_N * p.v;
        flux_Sig[j] = alpha * s.ftm.Sigma * p.v;
        flux_P[j] = alpha * s.ftm.P_z * p.v;
    }
}

void FiniteTempMatter::initialize_static(BSSNSlice& slice, const std::vector<FiniteTempPoint>& bg)
{
    const int N = (int)slice.states2.size();
    guess.assign(N, Point{});
    T_alpha0.assign(N, 0.);
    mu_alpha0.assign(N, 0.);
    for (int j = 0; j < N; j++)
    {
        State& s = slice.states2[j];
        const double s2 = s.csf.phi_re * s.csf.phi_re + s.csf.phi_im * s.csf.phi_im;
        const double T = (j < (int)bg.size()) ? bg[j].T : 0.;
        const double mu = (j < (int)bg.size()) ? bg[j].mu : 0.;
        const FiniteTempEoS::Thermo th = eos.thermal(T, mu, s2, -mu * mu);
        s.ftm.K_sigma = 0.;
        s.ftm.N_S = (s2 + th.dUps) * mu;
        s.ftm.N_N = th.n;
        s.ftm.Sigma = th.s;
        s.ftm.P_z = 0.;
        guess[j].T = T; guess[j].mu = mu; guess[j].v = 0.; guess[j].Pi = -mu;
        T_alpha0[j] = T * s.bssn.alpha;
        mu_alpha0[j] = mu * s.bssn.alpha;
    }
}

FTM FiniteTempMatter::rhs(BSSNSlice* slice, int j, double z, double d_z_alpha, double d_z_chi, double d_z_h_zz,
                          double d_z_h_ww, double d_z_beta, double chris_Zzz, double chris_Zww) const
{
    const State& s = slice->states2[j];
    const BSSNState& b = s.bssn;
    const Point& p = pts[j];
    const double D = st->D;
    const double d = D - 2.;
    const double min_z = st->min_z;
    const double alpha = b.alpha, beta = b.beta, K = b.K;
    const double chi = max(st->min_chi, b.chi);
    const double h_ZZ = 1. / b.h_zz, h_WW = 1. / b.h_ww;
    const double g_ZZ = chi * h_ZZ;

    // volume-element gradient used by Div[F] (exact for arbitrary h_zz, h_ww)
    const double g_ln = 0.5 * d_z_h_zz * h_ZZ + 0.5 * d * d_z_h_ww * h_WW - 0.5 * (D - 1.) * d_z_chi / chi;
    // divergence of an odd radial flux F (F/z -> F' at the origin) or of an even O(z^2) flux (F/z -> 0)
    auto div_odd = [&](const std::vector<double>& F) {
        const double Fz = deriv(slice, F, j, 1, true);
        return (z > min_z) ? (Fz + F[j] * (g_ln + d / z)) : ((d + 1.) * Fz + F[j] * g_ln);
    };
    auto div_even = [&](const std::vector<double>& F) {
        const double Fz = deriv(slice, F, j, 1, false);
        return (z > min_z) ? (Fz + F[j] * (g_ln + d / z)) : (Fz + F[j] * g_ln);
    };

    FTM r{};

    // condensate amplitude: d_t K_sigma = beta K_sigma' + alpha K K_sigma - (gamma^zz alpha' sigma' + alpha Lap sigma)/2
    //                                     + alpha (V' + Xi - delta F) sigma / 2
    const double sz = deriv(slice, f_sigma, j, 1, false);
    const double szz = deriv(slice, f_sigma, j, 2, false);
    const double sz_over_z = (z > min_z) ? sz / z : szz;
    const double lap = chi * (h_ZZ * (szz - chris_Zzz * sz) + d * h_WW * (sz_over_z - chris_Zww * sz))
                       - 0.5 * (D - 3.) * h_ZZ * d_z_chi * sz;
    const double sig = p.sigma;
    r.K_sigma = beta * slice->d_z(v_K_sigma, j, 1) + alpha * K * s.ftm.K_sigma
                - 0.5 * (g_ZZ * d_z_alpha * sz + alpha * lap)
                + 0.5 * alpha * (eos.dV(sig * sig) + p.Xi - p.dF) * sig;

    // superfluid charge
    // Div[alpha Upsilon gamma^zz psi'] written as Div[C X] with C = alpha gamma^zz Upsilon / sigma^2 (even) and
    // X = sigma^2 psi' = phi_re phi_im' - phi_im phi_re' (odd), using X' = phi_re phi_im'' - phi_im phi_re'' so that the
    // phase Laplacian uses compact second-derivative stencils (a derivative of a derivative leaves the 2h mode
    // undamped and goes unstable at the origin); at T = 0 this is exactly what the CSF Klein-Gordon equation discretizes.
    {
        const double C = coef_S[j], Cz = deriv(slice, coef_S, j, 1, false);
        const double X = p.cross;
        const double Xz = s.csf.phi_re * slice->d_z(v_phi_im, j, 2) - s.csf.phi_im * slice->d_z(v_phi_re, j, 2);
        const double div_S = (z > min_z) ? (Cz * X + C * Xz + C * X * (g_ln + d / z)) : ((d + 1.) * C * Xz + C * X * g_ln);
        r.N_S = beta * slice->d_z(v_N_S, j, 1) + alpha * K * s.ftm.N_S - div_S;
    }

    // normal number and entropy
    r.N_N = beta * slice->d_z(v_N_N, j, 1) + alpha * K * s.ftm.N_N - div_odd(flux_N);
    r.Sigma = beta * slice->d_z(v_Sigma, j, 1) + alpha * K * s.ftm.Sigma - div_odd(flux_Sig);

    // normal momentum (Valencia-type sources; pressure force s grad T + n grad mu)
    const double Tz = deriv(slice, f_T, j, 1, false);
    const double muz = deriv(slice, f_mu, j, 1, false);
    const double dz_gzz = d_z_h_zz / chi - b.h_zz * d_z_chi / (chi * chi);
    const double wW2 = p.w * p.W * p.W;
    r.P_z = beta * slice->d_z(v_P_z, j, 1) + alpha * K * s.ftm.P_z - div_even(flux_P) + s.ftm.P_z * d_z_beta
            + 0.5 * alpha * wW2 * p.v * p.v * dz_gzz - wW2 * d_z_alpha - alpha * (p.s * Tz + p.n * muz);
    if (p.theta > 0.)
    {
        // conversion (draft theta_1 terms); N_S + N_N is conserved, entropy is produced
        const double c = alpha * p.theta * p.J;
        r.N_S += c;
        r.N_N -= c;
        r.Sigma += c * p.J / p.T;
        r.P_z -= c * p.dz_psi;
    }
    if (fluid_atm[j]) { r.N_N = 0.; r.Sigma = 0.; r.P_z = 0.; }
    else if (atm_w[j] < 1.) { r.N_N *= atm_w[j]; r.Sigma *= atm_w[j]; r.P_z *= atm_w[j]; }
    return r;
}

void FiniteTempMatter::stress_energy(const BSSNState& b, int j, double& rho, double& j_z, double& S_zz, double& S_ww) const
{
    const Point& p = pts[j];
    const double chi = max(st->min_chi, b.chi);
    const double g_zz = b.h_zz / chi;
    const double wW2 = p.w * p.W * p.W;
    const double P_z = wW2 * g_zz * p.v;
    rho = wW2 - p.P_th + p.dUps * p.Pi * p.Pi;
    j_z = P_z - p.dUps * p.Pi * p.dz_psi;
    S_zz = P_z * g_zz * p.v + p.P_th * b.h_zz / chi + p.dUps * p.dz_psi * p.dz_psi;
    S_ww = p.P_th * b.h_ww / chi;
}

FTM FiniteTempMatter::atmosphere_rhs(const State& s, const State& r) const
{
    const CSF& c = s.csf;
    const CSF& d = r.csf;
    const double sg = std::sqrt(c.phi_re * c.phi_re + c.phi_im * c.phi_im);
    FTM out{0., 0., 0., 0., 0.};
    out.N_S = 2. * (d.phi_re * c.K_phi_im + c.phi_re * d.K_phi_im - d.phi_im * c.K_phi_re - c.phi_im * d.K_phi_re);
    if (sg > 1e-300)
    {
        const double sg_dot = (c.phi_re * d.phi_re + c.phi_im * d.phi_im) / sg;
        const double Ks = (c.phi_re * c.K_phi_re + c.phi_im * c.K_phi_im) / sg;
        out.K_sigma = (d.phi_re * c.K_phi_re + c.phi_re * d.K_phi_re + d.phi_im * c.K_phi_im + c.phi_im * d.K_phi_im) / sg
                      - Ks * sg_dot / sg;
    }
    return out;
}

// Implicit (backward-Euler) theta_1 conversion, applied by operator splitting after each full RK4 step:
// at every evolved-fluid point find the converted charge d with d = alpha theta_1 h J(U0 + d e), where
// e = (dN_S, dN_N, dSigma, dP_z) = (1, -1, J/T, -d_z psi). Unconditionally stable for any theta_1 h; for
// theta_1 h -> infinity it projects onto chemical equilibrium J = 0.
void FiniteTempMatter::relax_conversion(BSSNSlice* slice, double h)
{
    if (!(theta1 > 0.) || !(h > 0.)) return;
    prepare(slice);
    const int N = (int)slice->states2.size();
    double max_NN = 0.;
    for (int j = 0; j < N; j++)
        if (!(f_sigma[j] < fluid_sigma_min)) max_NN = max(max_NN, slice->states2[j].ftm.N_N);
    const double floor = max(n_floor, fluid_atm_frac * max_NN);
    for (int j = 0; j < N; j++)
    {
        if (atm[j] || fluid_atm[j]) continue;
        const Point p0 = pts[j];
        if (!(p0.T > 0.)) continue;
        State& st_j = slice->states2[j];
        const BSSNState& b = st_j.bssn;
        const FTM U0 = st_j.ftm;
        const double g_zz = b.h_zz / max(st->min_chi, b.chi);
        const double s2 = p0.sigma * p0.sigma;
        const double k = b.alpha * theta1 * h;
        const double J0 = p0.W * (p0.Pi + p0.v * p0.dz_psi) + p0.mu;
        const double JT0 = J0 / p0.T;
        FTM Ubest = U0;
        // evaluate g(d) = d - k J(U(d)); jt is the entropy factor J/T used for Sigma (trapezoidal in d)
        auto eval = [&](double d, double& g, double& jt, FTM& U) -> bool
        {
            double JT1 = JT0;
            Point q = p0;
            for (int pass = 0; pass < 2; pass++)
            {
                U = U0;
                U.N_S += d; U.N_N -= d; U.Sigma += d * 0.5 * (JT0 + JT1); U.P_z -= d * p0.dz_psi;
                if (!(U.N_N > floor) || !(U.Sigma > floor)) return false;
                double res = 0.;
                q = p0;
                if (!recover(U, floor, b, p0.sigma, p0.dz_sigma, p0.dz_psi, U0.K_sigma, q, res)) return false;
                if (!(q.T > 0.)) return false;
                const double W = 1. / std::sqrt(1. - min(1. - 1e-15, g_zz * q.v * q.v));
                const double J = W * (q.Pi + q.v * p0.dz_psi) + q.mu;
                JT1 = J / q.T;
                g = d - k * J;
            }
            jt = JT1;
            return true;
        };
        // linearised first guess: dJ/dd ~ -(1/(W n_mu) + 1/Upsilon)
        const double hmu = 1e-6 * max(1e-3, std::abs(p0.mu));
        const double n_mu = (eos.thermal(p0.T, p0.mu + hmu, s2, p0.Xi).n - eos.thermal(p0.T, p0.mu - hmu, s2, p0.Xi).n) / (2. * hmu);
        const double rate = 1. / max(1e-300, p0.W * n_mu) + 1. / max(1e-300, p0.Ups);
        if (!(n_mu > 0.) || !(p0.Ups > 0.)) continue;
        double da = 0., ga = 0. - k * J0;
        double db = k * J0 / (1. + k * rate), gb, jt;
        FTM U;
        bool ok = false;
        for (int it = 0; it < 30; it++)
        {
            int shrink = 0;
            while (!eval(db, gb, jt, U) && shrink < 20) { db = 0.5 * (da + db); shrink++; }
            if (shrink == 20) break;
            Ubest = U;
            if (std::abs(gb) <= 1e-12 * max(1e-12, std::abs(db)) + 1e-16) { ok = true; break; }
            const double denom = gb - ga;
            if (!(std::abs(denom) > 0.)) { ok = true; break; }
            const double dn = db - gb * (db - da) / denom;
            da = db; ga = gb; db = dn;
        }
        if (ok || Ubest.N_S != U0.N_S) st_j.ftm = Ubest;
    }
}
