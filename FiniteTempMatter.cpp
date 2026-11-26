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

    double B = 1. + (1. - c2) * mu * mu / Xs;
    if (std::abs(B) < 1e-12) B = (B >= 0.) ? 1e-12 : -1e-12;

    const double kph = M_PI * M_PI / 90.;
    const double T4 = T * T * T * T;
    const double Pph = kph * c * T4 / (B * B);
    // partial derivatives of P_ph
    const double Pph_T = 4. * Pph / T;
    const double Pph_mu = -2. * Pph / B * (2. * (1. - c2) * mu / Xs);
    const double Pph_c2 = Pph * (0.5 / c2 + 2. * mu * mu / (B * Xs));
    const double Pph_Xi = 2. * Pph * (1. - c2) * mu * mu / (B * Xs * Xs) + Pph_c2 * dc2_dXi;
    const double Pph_s2 = Pph_c2 * dc2_ds2;

    // free Boltzmann gas
    const double lnPg = 1.5 * std::log(m / (2. * M_PI)) + 2.5 * std::log(T) - (m - mu) / T;
    const double Pg_T_over_Pg = 2.5 / T + (m - mu) / (T * T);
    const double Pg_mu_over_Pg = 1. / T;

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
bool FiniteTempMatter::recover(const FTM& U, const BSSNState& b, double sigma, double dz_sigma, double dz_psi,
                               double K_sigma, Point& p, double& res) const
{
    (void)dz_sigma; (void)K_sigma; // the thermal EoS has no X dependence (delta Y = 0)
    const double chi = max(st->min_chi, b.chi);
    const double g_zz = b.h_zz / chi;
    const double g_ZZ = chi / b.h_zz;
    const double s2 = sigma * sigma;
    const double grad_psi2 = g_ZZ * dz_psi * dz_psi;

    const bool no_fluid = !(U.N_N > n_floor) || !(U.Sigma > n_floor);
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
    coef_S.resize(N); flux_S.resize(N); atm.assign(N, 0); flux_N.resize(N); flux_Sig.resize(N); flux_P.resize(N);
    if ((int)guess.size() != N) guess.assign(N, Point{});

    for (int j = 0; j < N; j++)
    {
        const CSF& c = slice->states2[j].csf;
        f_sigma[j] = std::sqrt(c.phi_re * c.phi_re + c.phi_im * c.phi_im);
    }

    c2p_failures = 0;
    max_c2p_res = 0.;
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
            f_T[j] = 0.;
            f_mu[j] = (j > 0) ? f_mu[j - 1] : p.mu; // only enters n mu' with n = 0 nearby; keep it smooth
            coef_S[j] = alpha * g_ZZ;
            flux_S[j] = coef_S[j] * cross;
            flux_N[j] = flux_Sig[j] = flux_P[j] = 0.;
            continue;
        }

        Point p = guess[j];
        double res = 0.;
        if (!recover(s.ftm, b, sig, dz_sigma, dz_psi, s.ftm.K_sigma, p, res))
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
        pts[j] = p;

        // derived SBSE momenta: K_phi = (phi / sigma) K_sigma - (i / 2) phi Pi_psi
        const double ks = (sig > 1e-300) ? s.ftm.K_sigma / sig : 0.;
        s.csf.K_phi_re = s.csf.phi_re * ks + 0.5 * s.csf.phi_im * p.Pi;
        s.csf.K_phi_im = s.csf.phi_im * ks - 0.5 * s.csf.phi_re * p.Pi;

        f_T[j] = p.T;
        f_mu[j] = p.mu;
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
