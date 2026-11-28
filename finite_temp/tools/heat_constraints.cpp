// Out-of-equilibrium "heated" initial data: take a relaxed finite-T profile, raise T_inf to a new value keeping the
// condensate (sigma(r), omega) fixed, and re-solve only the polar-areal metric equations
//     m' = 4 pi r^2 rho,  Xs = (1 - 2m/r)^{-1/2}               (Hamiltonian constraint)
//     (ln alpha)' = (Xs^2 - 1)/(2r) + 4 pi r Xs^2 (Y eta^2 + P)  (polar slicing; a gauge choice, alpha_N Xs_N = 1)
// with T_local = T_inf/alpha, mu = omega/alpha (Tolman/Klein) evaluated self-consistently. The normal fluid is at rest
// (momentum constraint satisfied), but the configuration is not a static solution of the field equations.
// Usage: heat_constraints FiniteTempRelaxer.par in_profile.dat T_new out_profile.dat [polar = 0] [refine = 1]
//   polar = 0 (default): keep the input lapse, solve only the Hamiltonian constraint; 1: also re-solve polar slicing
//   with alpha_N Xs_N = 1 by shooting on the central lapse
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#define private public
#include "../FiniteTempRelaxer.h"
#undef private

int main(int argc, char** argv)
{
    if (argc < 5) { std::cerr << "usage: heat_constraints par in_profile T_new out_profile\n"; return 1; }
    FiniteTempRelaxer R;
    R.read_parameters(argv[1], true);
    if (!R.load_profile(argv[2], true)) { std::cerr << "cannot load " << argv[2] << "\n"; return 1; }
    const double T0 = R.current_temperature_, Tn = std::atof(argv[3]);
    R.current_temperature_ = Tn;
    // optional grid refinement (factor argv[6]): sigma by cubic Hermite with its own slope sigma' = eta Xs, the lapse by
    // cubic Hermite with centred-difference slopes (both in areal r); Xs is re-solved below anyway
    const double refine = (argc > 6) ? std::atof(argv[6]) : 1.;
    if (refine > 1.)
    {
        const RelaxerState old = R.state_;
        const int No = (int)old.rbar.size();
        std::vector<double> ro(No), so(No), dso(No), ao(No), dao(No);
        for (int i = 0; i < No; i++) { ro[i] = R.areal_radius(old.rbar[i]); so[i] = old.sigma[i]; dso[i] = old.eta[i] * old.Xs[i]; ao[i] = old.alpha[i]; }
        for (int i = 0; i < No; i++)
        {
            const int a = std::max(0, i - 1), b = std::min(No - 1, i + 1);
            dao[i] = (i == 0) ? 0. : (ao[b] - ao[a]) / (ro[b] - ro[a]);
        }
        R.refine_grid(refine);
        RelaxerState& Sn = R.state_;
        for (size_t i = 0; i < Sn.rbar.size(); i++)
        {
            const double rr = R.areal_radius(Sn.rbar[i]);
            int k = (int)(std::upper_bound(ro.begin(), ro.end(), rr) - ro.begin()) - 1;
            k = std::max(0, std::min(No - 2, k));
            const double h = ro[k + 1] - ro[k], t = (rr - ro[k]) / h;
            const double h00 = (1 + 2 * t) * (1 - t) * (1 - t), h10 = t * (1 - t) * (1 - t), h01 = t * t * (3 - 2 * t), h11 = t * t * (t - 1);
            const double dh00 = 6 * t * t - 6 * t, dh10 = 3 * t * t - 4 * t + 1, dh01 = -6 * t * t + 6 * t, dh11 = 3 * t * t - 2 * t;
            Sn.sigma[i] = h00 * so[k] + h10 * h * dso[k] + h01 * so[k + 1] + h11 * h * dso[k + 1];
            const double ds = (dh00 * so[k] + dh10 * h * dso[k] + dh01 * so[k + 1] + dh11 * h * dso[k + 1]) / h;
            Sn.alpha[i] = h00 * ao[k] + h10 * h * dao[k] + h01 * ao[k + 1] + h11 * h * dao[k + 1];
            Sn.Xs[i] = 1.; Sn.eta[i] = ds; // eta = sigma' / Xs, fixed below once Xs is known
            Sn.omega[i] = old.omega.back();
        }
        R.state_.Xs.assign(R.state_.rbar.size(), 1.);
    }
    RelaxerState& S = R.state_;
    const int N = (int)S.rbar.size();
    std::vector<double> r(N), dsig(N);
    for (int i = 0; i < N; i++) { r[i] = R.areal_radius(S.rbar[i]); dsig[i] = S.eta[i] * S.Xs[i]; } // sigma' (fixed)
    const double om = S.omega.back();
    // polar = 1: outward integration for a given central lapse a_c (Heun steps on the relaxer grid): m and ln alpha together,
    // so rho(T_inf/alpha, omega/alpha) uses the alpha reached so far. Then a_c is chosen so that alpha_N Xs_N = 1.
    std::vector<double> m(N, 0.), lna(N, 0.), Xsv(N, 1.);
    auto src = [&](int i, double a, double Xs, double& rho, double& Srr) {
        a = std::max(1e-12, a);
        const double eta = dsig[i] / Xs;
        const FiniteTempRelaxer::PointThermo th = R.thermo_at(S.sigma[i], eta, a, om, r[i]);
        const double s2 = S.sigma[i] * S.sigma[i], X = eta * eta;
        const FiniteTempRelaxer::ThermoConjugates tc = R.conjugates_analytic(Tn / a, om / a, s2, X, -(om / a) * (om / a), r[i]);
        rho = th.rho; Srr = tc.Y * X + th.P;
    };
    auto rhs = [&](int i, double mm, double la, double& dm, double& dla) {
        const double Xs = (i == 0) ? 1. : 1. / std::sqrt(std::max(1e-12, 1. - 2. * mm / r[i]));
        double rho, Srr;
        src(i, std::exp(la), Xs, rho, Srr);
        dm = 4. * M_PI * r[i] * r[i] * rho;
        dla = (i == 0) ? 0. : (Xs * Xs - 1.) / (2. * r[i]) + 4. * M_PI * r[i] * Xs * Xs * Srr;
        return std::isfinite(dm) && std::isfinite(dla) && (i == 0 || 2. * mm < r[i]);
    };
    auto shoot = [&](double ac) -> double {
        m[0] = 0.; lna[0] = std::log(ac);
        for (int i = 1; i < N; i++)
        {
            const double h = r[i] - r[i - 1];
            double dm0, dl0, dm1, dl1;
            if (!rhs(i - 1, m[i - 1], lna[i - 1], dm0, dl0)) return NAN;
            const double mp = m[i - 1] + h * dm0, lp = lna[i - 1] + h * dl0;
            if (!rhs(i, mp, lp, dm1, dl1)) return NAN;
            m[i] = m[i - 1] + 0.5 * h * (dm0 + dm1); lna[i] = lna[i - 1] + 0.5 * h * (dl0 + dl1);
            if (2. * m[i] >= r[i]) return NAN;
        }
        const double XsN = 1. / std::sqrt(1. - 2. * m[N - 1] / r[N - 1]);
        return lna[N - 1] + std::log(XsN); // = ln(alpha_N Xs_N)
    };
    const bool polar = (argc > 5) && std::atoi(argv[5]) == 1;
    if (!polar)
    {
        // default: keep the equilibrium lapse (so T_local, mu_local are the old ones scaled by T_new/T_old and unchanged),
        // and solve only the Hamiltonian constraint for Xs (Heun)
        auto dm_at = [&](int i, double mm) {
            const double Xs = (i == 0) ? 1. : 1. / std::sqrt(std::max(1e-12, 1. - 2. * mm / r[i]));
            double rho, Srr; src(i, S.alpha[i], Xs, rho, Srr);
            return 4. * M_PI * r[i] * r[i] * rho;
        };
        m[0] = 0.;
        for (int i = 1; i < N; i++)
        {
            const double h = r[i] - r[i - 1];
            const double d0 = dm_at(i - 1, m[i - 1]);
            const double d1 = dm_at(i, m[i - 1] + h * d0);
            m[i] = m[i - 1] + 0.5 * h * (d0 + d1);
            if (!(2. * m[i] < r[i])) { std::cerr << "trapped surface in the constraint solution at r = " << r[i] << "\n"; return 2; }
        }
        for (int i = 0; i < N; i++) S.Xs[i] = (i == 0) ? 1. : 1. / std::sqrt(1. - 2. * m[i] / r[i]);
        for (int i = 0; i < N; i++) S.eta[i] = dsig[i] / S.Xs[i];
        R.write_profile(argv[4]);
        std::cerr << "T_inf " << T0 << " -> " << Tn << " (lapse kept), sigma0 " << S.sigma[0] << ", omega " << om << ", M(r_max) " << m[N - 1] << "\n";
        return 0;
    }
    // bracket: scan a_c on a grid, then bisect the sign change closest to the input central lapse
    const double a_in = S.alpha[0];
    double best_lo = -1., best_hi = -1.;
    double prev_a = 0., prev_g = NAN;
    for (int k = 1; k <= 400; k++)
    {
        const double a = k / 400.;
        const double g = shoot(a);
        if (std::isfinite(g) && std::isfinite(prev_g) && (g > 0.) != (prev_g > 0.))
        {
            if (best_lo < 0. || std::abs(0.5 * (a + prev_a) - a_in) < std::abs(0.5 * (best_lo + best_hi) - a_in)) { best_lo = prev_a; best_hi = a; }
        }
        std::cerr << "a_c " << a << " g " << g << " M " << m[N - 1] << "\n";
        prev_a = a; prev_g = g;
    }
    if (best_lo < 0.) { std::cerr << "no root: no polar-sliced constraint solution at this T\n"; return 2; }
    double lo = best_lo, hi = best_hi, glo = shoot(lo);
    for (int it = 0; it < 80; it++)
    {
        const double mid = 0.5 * (lo + hi), g = shoot(mid);
        if ((g > 0.) == (glo > 0.)) { lo = mid; glo = g; } else hi = mid;
    }
    shoot(0.5 * (lo + hi));
    for (int i = 0; i < N; i++)
    {
        S.alpha[i] = std::exp(lna[i]);
        S.Xs[i] = (i == 0) ? 1. : 1. / std::sqrt(1. - 2. * m[i] / r[i]);
    }
    for (int i = 0; i < N; i++) S.eta[i] = dsig[i] / S.Xs[i];
    R.write_profile(argv[4]);
    std::cerr << "T_inf " << T0 << " -> " << Tn << ", sigma0 " << S.sigma[0] << ", omega " << om << ", M(r_max) " << m[N - 1] << "\n";
    return 0;
}
