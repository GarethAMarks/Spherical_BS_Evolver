#ifndef BOGOLIUBOV_TABLE_H
#define BOGOLIUBOV_TABLE_H
// Tabulated Bogoliubov thermal pressure for the evolution (bogoliubov.h evaluates f by quadrature, ~16 us per call,
// too slow for primitive recovery at every point and RK stage).
//
// Table coordinates (all smooth across the physical domain):
//   x    = ln(1 + u), u = sqrt(delta (delta + 2 gam))   u is the quasiparticle gap over T (delta = phase-mode gap
//          parameter / T). At large gam, f varies with delta on the scale 1/gam, and the Landau limit has a
//          sqrt(delta) term, but both are smooth in u.
//   t    = ln(1 + gam)          gam = a / (2 mu_c T)
//   zeta = eta / eta_L(delta, gam)   counterflow as a fraction of the Landau limit, eta_L = min_k Omega^2 / (2 k^2)
//          = (u + delta)(u + delta + 2 gam) / (4 u), u = sqrt(delta (delta + 2 gam))  (gam/2 at delta = 0, delta at gam = 0)
// Stored: Lt = ln f + (1 - zeta) delta + 2 ln(1 - zeta) + 1.5 t, which removes the Boltzmann factor e^{eta - delta}, the
// phonon pole (1 - zeta)^{-2} and the phonon scaling gam^{-3/2}, so Lt varies slowly. Interpolation: tricubic
// B-spline through the nodes (C^2; a C^1 interpolant made the primitive-recovery Newton stall). Outside the box ln f is continued linearly (C^1) in each
// coordinate; in particular beyond zeta_max (close to and past the Landau limit) P stays finite and convex in eta.
#include "bogoliubov.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace bogoliubov
{
inline void landau_eta(double delta, double gam, double& etaL, double& etaL_d, double& etaL_g)
{
    // envelope theorem: derivatives of g(u) = (u + d)(u + d + 2 gam)/(4u) at the minimiser u*
    const double d = std::max(0., delta), g = std::max(0., gam);
    const double u = std::sqrt(d * (d + 2. * g));
    if (u < 1e-300)
    {
        etaL = 0.5 * g; etaL_d = (g > 0.) ? 1e300 : 1.; etaL_g = 0.5; return;
    }
    etaL = (u + d) * (u + d + 2. * g) / (4. * u);
    etaL_d = (2. * u + 2. * d + 2. * g) / (4. * u);
    etaL_g = (u + d) / (2. * u);
}

struct Table
{
    double s_max = 6.6, t_max = 7.6, z_min = -0.5, z_max = 0.9;
    int Ns = 111, Nt = 77, Nz = 71;
    double hs, ht, hz;
    std::vector<double> L; // [i][j][k] -> (i*Nt + j)*Nz + k
    std::vector<double> D; // d ln f / d delta at fixed (gam, eta), same layout
    bool built = false;

    void build()
    {
        hs = s_max / (Ns - 1); ht = t_max / (Nt - 1); hz = (z_max - z_min) / (Nz - 1);
        L.assign((size_t)Ns * Nt * Nz, 0.);
        D.assign((size_t)Ns * Nt * Nz, 0.);
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int i = 0; i < Ns; i++)
            for (int j = 0; j < Nt; j++)
            {
                const double t = j * ht, gam = std::expm1(t), u = std::expm1(i * hs);
                const double den = gam + std::sqrt(gam * gam + u * u);
                const double delta = (den > 0.) ? u * u / den : 0.; // = -gam + sqrt(gam^2 + u^2)
                double eL, eLd, eLg;
                landau_eta(delta, gam, eL, eLd, eLg);
                // step for d ln f / d delta: f varies with delta on the scale 1/(1 + 2 gam) near delta = 0
                const double hd = 1e-4 * (delta + 1. / (1. + 2. * gam));
                for (int k = 0; k < Nz; k++)
                {
                    const double z = z_min + k * hz, eta = z * eL;
                    const double lf = lnf_quad(delta, gam, eta);
                    const size_t id = ((size_t)i * Nt + j) * Nz + k;
                    L[id] = lf + (1. - z) * delta + 2. * std::log1p(-z) + 1.5 * t;
                    double dl;
                    if (delta - hd >= 0.)
                        dl = (lnf_quad(delta + hd, gam, eta) - lnf_quad(delta - hd, gam, eta)) / (2. * hd);
                    else
                        dl = (-3. * lf + 4. * lnf_quad(delta + hd, gam, eta) - lnf_quad(delta + 2. * hd, gam, eta)) / (2. * hd);
                    D[id] = dl;
                }
            }
        // nodes where the quadrature is unreliable (large gap at eta < 0: cancellation) are filled from the nearest
        // reliable node at smaller gap (the stored Lt has the Boltzmann factor removed, so it is nearly constant there)
        for (int j = 0; j < Nt; j++)
            for (int k = 0; k < Nz; k++)
                for (int i = 1; i < Ns; i++)
                {
                    const size_t id = ((size_t)i * Nt + j) * Nz + k, ip = ((size_t)(i - 1) * Nt + j) * Nz + k;
                    if (!std::isfinite(L[id])) L[id] = L[ip];
                    if (!std::isfinite(D[id])) D[id] = D[ip];
                }
        make_splines();
        built = true;
    }

    // ---- cubic B-spline interpolation (C^2, so the c2p Jacobian, built from second derivatives of P, is continuous) ----
    // Coefficients on the extended index range -1..N per axis, natural end conditions (zero second derivative at the
    // end nodes). Layout: ((i+1)*(Nt+2) + (j+1))*(Nz+2) + (k+1).
    std::vector<double> CL, CD;
    static void spline1d(std::vector<double>& y, std::vector<double>& c)
    {
        // y: N node values -> c: N+2 coefficients (c[0] is index -1)
        const int N = (int)y.size();
        c.assign(N + 2, 0.);
        c[1] = y[0]; c[N] = y[N - 1];
        const int M = N - 2; // interior unknowns c_1..c_{N-2} (array 2..N-1)
        if (M > 0)
        {
            std::vector<double> cp(M), dp(M);
            for (int q = 0; q < M; q++)
            {
                double rhs = 6. * y[q + 1];
                if (q == 0) rhs -= y[0];
                if (q == M - 1) rhs -= y[N - 1];
                const double b = 4., a = (q > 0) ? 1. : 0., cc = (q < M - 1) ? 1. : 0.;
                const double den = b - a * ((q > 0) ? cp[q - 1] : 0.);
                cp[q] = cc / den;
                dp[q] = (rhs - a * ((q > 0) ? dp[q - 1] : 0.)) / den;
            }
            for (int q = M - 1; q >= 0; q--) c[q + 2] = dp[q] - ((q < M - 1) ? cp[q] * c[q + 3] : 0.);
        }
        c[0] = 2. * c[1] - c[2];
        c[N + 1] = 2. * c[N] - c[N - 1];
    }
    void make_coeffs(const std::vector<double>& V, std::vector<double>& C) const
    {
        const int A = Ns + 2, B = Nt + 2, Z = Nz + 2;
        std::vector<double> T1((size_t)Ns * Nt * Z), T2((size_t)Ns * B * Z);
        C.assign((size_t)A * B * Z, 0.);
        std::vector<double> y, c;
        for (int i = 0; i < Ns; i++)
            for (int j = 0; j < Nt; j++)
            {
                y.assign(Nz, 0.);
                for (int k = 0; k < Nz; k++) y[k] = V[((size_t)i * Nt + j) * Nz + k];
                spline1d(y, c);
                for (int k = 0; k < Z; k++) T1[((size_t)i * Nt + j) * Z + k] = c[k];
            }
        for (int i = 0; i < Ns; i++)
            for (int k = 0; k < Z; k++)
            {
                y.assign(Nt, 0.);
                for (int j = 0; j < Nt; j++) y[j] = T1[((size_t)i * Nt + j) * Z + k];
                spline1d(y, c);
                for (int j = 0; j < B; j++) T2[((size_t)i * B + j) * Z + k] = c[j];
            }
        for (int j = 0; j < B; j++)
            for (int k = 0; k < Z; k++)
            {
                y.assign(Ns, 0.);
                for (int i = 0; i < Ns; i++) y[i] = T2[((size_t)i * B + j) * Z + k];
                spline1d(y, c);
                for (int i = 0; i < A; i++) C[((size_t)i * B + j) * Z + k] = c[i];
            }
    }
    void make_splines() { make_coeffs(L, CL); make_coeffs(D, CD); }

    static void bs(double x, double* w, double* dw)
    {
        const double x2 = x * x, x3 = x2 * x, om = 1. - x;
        w[0] = om * om * om / 6.;
        w[1] = (3. * x3 - 6. * x2 + 4.) / 6.;
        w[2] = (-3. * x3 + 3. * x2 + 3. * x + 1.) / 6.;
        w[3] = x3 / 6.;
        dw[0] = -0.5 * om * om;
        dw[1] = (9. * x2 - 12. * x) / 6.;
        dw[2] = (-9. * x2 + 6. * x + 3.) / 6.;
        dw[3] = 0.5 * x2;
    }
    // spline value and partials in (x, t, z) at a point inside the box (C = CL or CD)
    void eval(const std::vector<double>& C, double s, double t, double z, double& v, double& vs, double& vt, double& vz) const
    {
        const double xs = s / hs, xt = t / ht, xz = (z - z_min) / hz;
        const int i = std::max(0, std::min(Ns - 2, (int)std::floor(xs)));
        const int j = std::max(0, std::min(Nt - 2, (int)std::floor(xt)));
        const int k = std::max(0, std::min(Nz - 2, (int)std::floor(xz)));
        double ws[4], dws[4], wt[4], dwt[4], wz[4], dwz[4];
        bs(xs - i, ws, dws); bs(xt - j, wt, dwt); bs(xz - k, wz, dwz);
        const int B = Nt + 2, Z = Nz + 2;
        v = vs = vt = vz = 0.;
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
            {
                // coefficient index (i-1+a) -> array (i+a), etc.
                const double* row = &C[((size_t)(i + a) * B + (j + b)) * Z + k];
                double c = 0., cz = 0.;
                for (int e = 0; e < 4; e++) { c += wz[e] * row[e]; cz += dwz[e] * row[e]; }
                v += ws[a] * wt[b] * c; vs += dws[a] * wt[b] * c; vt += ws[a] * dwt[b] * c; vz += ws[a] * wt[b] * cz;
            }
        vs /= hs; vt /= ht; vz /= hz;
    }
    void interp(double s, double t, double z, double& v, double& vs, double& vt, double& vz) const { eval(CL, s, t, z, v, vs, vt, vz); }
    double interpD(double s, double t, double z) const { double v, a, b, c; eval(CD, s, t, z, v, a, b, c); return v; }

    // ln f(delta, gam, eta) and its partial derivatives. Linear continuation outside the box.
    void lnf(double delta, double gam, double eta, double& F, double& F_d, double& F_g, double& F_e) const
    {
        delta = std::max(0., delta); gam = std::max(0., gam);
        double eL, eLd, eLg;
        landau_eta(delta, gam, eL, eLd, eLg);
        const double u = std::sqrt(delta * (delta + 2. * gam));
        const double x = std::log1p(u), t = std::log1p(gam);
        const double xc = std::min(x, s_max), tc = std::min(t, t_max);
        // zeta; at delta = gam = 0 (eta_L = 0) only eta <= 0 is physical: use a tiny eta_L floor
        const double eLs = std::max(eL, 1e-12);
        const double z = eta / eLs;
        const double zc = std::min(std::max(z, z_min), z_max);
        double v, vx, vt, vz;
        interp(xc, tc, zc, v, vx, vt, vz);
        // ln f = Lt - (1 - z) delta(u, gam) - 2 ln(1 - z) - 1.5 t at the clamped point, and its partials in (x, t, z)
        const double gc = std::expm1(tc), uc = std::expm1(xc);
        const double sq = std::sqrt(gc * gc + uc * uc);
        const double dc = (gc + sq > 0.) ? uc * uc / (gc + sq) : 0.;
        const double dd_du = (sq > 0.) ? uc / sq : 1., dd_dg = (sq > 0.) ? gc / sq - 1. : 0.;
        double G = v - (1. - zc) * dc - 2. * std::log1p(-zc) - 1.5 * tc;
        double G_x = vx - (1. - zc) * dd_du * (1. + uc);
        double G_t = vt - 1.5 - (1. - zc) * dd_dg * (1. + gc);
        double G_z = vz + dc + 2. / (1. - zc);
        // linear continuation
        G += G_x * (x - xc) + G_t * (t - tc) + G_z * (z - zc);
        F = G;
        // chain rule: x = ln(1 + u(delta, gam)), t = ln(1 + gam), z = eta / eta_L(delta, gam).
        // du/d(delta) = (delta + gam)/u is singular at u = 0, and the singular parts of the x- and z-derivatives cancel; so
        // d ln f / d delta is tabulated separately (D) and the gam-derivative at fixed delta is assembled from the one at
        // fixed u: d/dgam|_delta = d/dgam|_u + (d/du) du/dgam|_delta, with (d/du) du/dgam = F_d (d delta/du)^-1 du/dgam = F_d delta/(delta + gam).
        const double dz_de = 1. / eLs;
        F_d = interpD(xc, tc, zc);
        const double eLg_u = eLg + eLd * ((std::sqrt(gam * gam + u * u) > 0.) ? (gam / std::sqrt(gam * gam + u * u) - 1.) : 0.);
        const double dz_dg_u = (eL > 1e-12) ? -z * eLg_u / eL : 0.;
        F_g = G_t / (1. + gam) + G_z * dz_dg_u + ((delta + gam > 0.) ? F_d * delta / (delta + gam) : 0.);
        F_e = G_z * dz_de;
    }
};

// The table is built on first use (about 1M quadratures; ~1 min serial) and cached in the file named by the environment
// variable BOGOLIUBOV_TABLE (default "bogoliubov_table.bin" in the working directory); a cache with different grid
// parameters is ignored and rewritten.
inline const Table& table()
{
    static Table T = []() {
        Table t;
        const char* env = std::getenv("BOGOLIUBOV_TABLE");
        const std::string path = env ? env : "bogoliubov_table.bin";
        const double hdr[8] = {4. /* format version */, t.s_max, t.t_max, t.z_min, t.z_max, (double)t.Ns, (double)t.Nt, (double)t.Nz};
        const size_t n = (size_t)t.Ns * t.Nt * t.Nz;
        if (std::FILE* fp = std::fopen(path.c_str(), "rb"))
        {
            double h[8];
            bool ok = std::fread(h, sizeof(double), 8, fp) == 8;
            for (int q = 0; ok && q < 8; q++) ok = (h[q] == hdr[q]);
            if (ok)
            {
                t.L.resize(n); t.D.resize(n);
                ok = std::fread(t.L.data(), sizeof(double), n, fp) == n && std::fread(t.D.data(), sizeof(double), n, fp) == n;
            }
            std::fclose(fp);
            if (ok)
            {
                t.hs = t.s_max / (t.Ns - 1); t.ht = t.t_max / (t.Nt - 1); t.hz = (t.z_max - t.z_min) / (t.Nz - 1);
                t.make_splines();
                t.built = true;
                return t;
            }
        }
        t.build();
        if (std::FILE* fp = std::fopen(path.c_str(), "wb"))
        {
            std::fwrite(hdr, sizeof(double), 8, fp);
            std::fwrite(t.L.data(), sizeof(double), n, fp);
            std::fwrite(t.D.data(), sizeof(double), n, fp);
            std::fclose(fp);
        }
        return t;
    }();
    return T;
}

// Thermal pressure and its partials from the table: P(T, mu, mu_c, a, gfac) with mass m and gap softening width w
// (see gap_delta and pressure in bogoliubov.h). Returns P and dP/dT, dP/dmu, dP/da, dP/dmu_c, dP/dgfac.
inline void pressure_table(double T, double mu, double mu_c, double a, double m, double gap_width, double gfac,
                           double& P, double& P_T, double& P_mu, double& P_a, double& P_c, double& P_g)
{
    P = P_T = P_mu = P_a = P_c = P_g = 0.;
    if (!(T > 0.) || !(mu_c > 0.)) return;
    a = std::max(0., a);
    double sig; // d delta / d delta_lin
    const double delta = gap_delta(T, mu_c, m, gap_width, gfac, &sig);
    const double d_lin = (m - mu_c) * gfac / T;
    const double dlin_T = -d_lin / T, dlin_c = -gfac / T, dlin_g = (m - mu_c) / T;
    const double M = kinetic_mass(mu_c, m, gfac), M_c = 1. - gfac, M_g = m - mu_c;
    const double gam = a / (2. * M * T), eta = (mu - mu_c) / T;
    double F, F_d, F_g, F_e;
    table().lnf(delta, gam, eta, F, F_d, F_g, F_e);
    const double lnpre = 1.5 * std::log(M * T / (2. * M_PI)) + std::log(T);
    P = std::exp(lnpre + F);
    // d ln P; gam depends on M
    const double dlP_M = 1.5 / M + F_g * (-gam / M);
    const double dlP_T = 2.5 / T + F_d * sig * dlin_T + F_g * (-gam / T) + F_e * (-eta / T);
    const double dlP_mu = F_e / T;
    const double dlP_a = F_g / (2. * M * T);
    const double dlP_c = dlP_M * M_c + F_d * sig * dlin_c + F_e * (-1. / T);
    const double dlP_g = dlP_M * M_g + F_d * sig * dlin_g;
    P_T = P * dlP_T; P_mu = P * dlP_mu; P_a = P * dlP_a; P_c = P * dlP_c; P_g = P * dlP_g;
}
} // namespace bogoliubov
#endif
