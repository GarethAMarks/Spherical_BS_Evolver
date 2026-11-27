#ifndef BOGOLIUBOV_H
#define BOGOLIUBOV_H
// One-loop thermal pressure of Bogoliubov quasiparticles of a locally uniform condensate (notes, section
// "A unified thermal EoS from the Bogoliubov spectrum"). Non-relativistic scaling form
//     P_th = (mu_c T / 2pi)^{3/2} T f(delta, gam, eta),
//     delta = Delta/T, Delta ~ (m - mu_c) exp(-(sigma/sigma_gap)^4), kept >= 0 (phase-mode gap: zero where there is a
//             condensate, m - mu_c for sigma -> 0; see gap_delta)
//     gam   = g/T,     g = a/(2 mu_c), a = sigma^2 V''(sigma^2)
//     eta   = (mu - mu_c)/T  (counterflow, mu = gamma mu_c; eta < 0 is the analytic continuation)
//     f = (2/sqrt(2pi)) int_0^inf kappa^2 G(Omega(kappa), 2 eta kappa^2) dkappa,
//     Omega = sqrt(e(e + 2 gam)), e = kappa^2/2 + delta,
//     G(Omega, b^2) = [Li2(e^{-Omega+b}) - Li2(e^{-Omega-b})]/(2b) = sum_k e^{-k Omega} sinh(k b)/(k^2 b),
// which is even in b, so G depends on b^2 = 2 eta kappa^2 and continues analytically to eta < 0.
// Limits: gam = 0 gives f = Li_{5/2}(e^{eta - delta}) (free Bose gas), gam >> 1, delta = eta = 0 gives
// f = (2pi)^{3/2} pi^2 / (90 gam^{3/2}) (phonons). The Landau limit is min_kappa [Omega(kappa) - b(kappa)] = 0.
// Quadrature: fixed Gauss-Legendre nodes on a smooth map kappa = L s/(1-s), so f is smooth in its arguments
// (finite differences of f are clean).
#include <cmath>
#include <complex>

namespace bogoliubov
{
// Real dilogarithm for x <= 1.
inline double li2(double x)
{
    if (x == 0.) return 0.;
    if (x == 1.) return M_PI * M_PI / 6.;
    if (x < -1.) { const double l = std::log(-x); return -M_PI * M_PI / 6. - 0.5 * l * l - li2(1. / x); }
    if (x > 0.5) return M_PI * M_PI / 6. - std::log(x) * std::log1p(-x) - li2(1. - x);
    if (x < -0.5) { return -li2(-x / (1. - x)) - 0.5 * std::log1p(-x) * std::log1p(-x); }
    // |x| <= 0.5: Bernoulli series in u = -ln(1-x)
    static const double B[] = {1., -0.25, 1. / 36., -1. / 3600., 1. / 211680., -1. / 10886400., 1. / 526901760.,
        -4.064761645144225524e-11, 8.921691020456452555e-13, -1.993929586072107569e-14, 4.518980029619918192e-16};
    const double u = -std::log1p(-x), u2 = u * u;
    double s = 0., up = u;
    s = B[0] * u + B[1] * u2;
    up = u;
    for (int k = 2; k < 11; k++) { up *= u2; s += B[k] * up; }
    return s;
}

// Complex dilogarithm for |z| < 1 (used for eta < 0).
inline std::complex<double> cli2(std::complex<double> z)
{
    typedef std::complex<double> C;
    static const double B[] = {1., -0.25, 1. / 36., -1. / 3600., 1. / 211680., -1. / 10886400., 1. / 526901760.,
        -4.064761645144225524e-11, 8.921691020456452555e-13, -1.993929586072107569e-14, 4.518980029619918192e-16};
    if (std::abs(z) < 1e-300) return C(0., 0.);
    if (z.real() > 0.5)
    {
        // Li2(z) = pi^2/6 - ln z ln(1-z) - Li2(1-z)
        return C(M_PI * M_PI / 6., 0.) - std::log(z) * std::log(C(1., 0.) - z) - cli2(C(1., 0.) - z);
    }
    if (std::abs(z) > 1.)
    {
        // Li2(z) = -pi^2/6 - ln^2(-z)/2 - Li2(1/z)
        const C l = std::log(-z);
        return C(-M_PI * M_PI / 6., 0.) - 0.5 * l * l - cli2(C(1., 0.) / z);
    }
    const C u = -std::log(C(1., 0.) - z), u2 = u * u;
    C s = B[0] * u + B[1] * u2, up = u;
    for (int k = 2; k < 11; k++) { up *= u2; s += B[k] * up; }
    return s;
}

// G(Omega, b2) as defined above, Omega > 0, b2 = b^2 (any sign), requires Omega > sqrt(b2) when b2 > 0.
inline double Gfun(double Om, double b2)
{
    const double q = std::exp(-Om);
    if (std::fabs(b2) < 1e-6 * (1. + Om * Om))
    {
        // Taylor in b^2: G = Li1(q) + b^2 Li_{-1}(q)/6 + b^4 Li_{-3}(q)/120
        const double om = -std::expm1(-Om);
        const double L1 = -std::log(om);
        const double Lm1 = q / (om * om);
        const double Lm3 = q * (1. + 4. * q + q * q) / (om * om * om * om);
        return L1 + b2 * Lm1 / 6. + b2 * b2 * Lm3 / 120.;
    }
    if (b2 > 0.)
    {
        const double b = std::sqrt(b2);
        return (li2(std::exp(-Om + b)) - li2(std::exp(-Om - b))) / (2. * b);
    }
    const double be = std::sqrt(-b2);
    const std::complex<double> zp = std::exp(std::complex<double>(-Om, be));
    // [Li2(q e^{i be}) - Li2(q e^{-i be})]/(2 i be) = Im Li2(q e^{i be}) / be
    return cli2(zp).imag() / be;
}

static const int kNodes = 160;

struct GLNodes
{
    double xw[2 * kNodes];
    GLNodes()
    {
        // Gauss-Legendre on [0,1]: xw[2i] = node, xw[2i+1] = weight
        for (int i = 0; i < kNodes; i++)
        {
            double x = std::cos(M_PI * (i + 0.75) / (kNodes + 0.5)), dp = 0.;
            for (int it = 0; it < 100; it++)
            {
                double p0 = 1., p1 = x;
                for (int k = 2; k <= kNodes; k++) { const double p2 = ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k; p0 = p1; p1 = p2; }
                dp = kNodes * (x * p1 - p0) / (x * x - 1.);
                const double dx = p1 / dp;
                x -= dx;
                if (std::fabs(dx) < 1e-15) break;
            }
            xw[2 * i] = 0.5 * (1. - x);
            xw[2 * i + 1] = 1. / ((1. - x * x) * dp * dp);
        }
    }
};

// Thread-safe (function-local static initialisation).
inline const double* gl_nodes()
{
    static const GLNodes nodes;
    return nodes.xw;
}

// f(delta, gam, eta). Returns +inf beyond the Landau limit.
inline double f(double delta, double gam, double eta)
{
    const double* xw = gl_nodes();
    const double L = 2.0 / std::sqrt(1. + gam) * std::pow(1. + 4. * eta * eta, 0.25); // smooth in eta (FD across eta = 0)
    double sum = 0.;
    for (int i = 0; i < kNodes; i++)
    {
        const double s = xw[2 * i], w = xw[2 * i + 1];
        const double kap = L * s / (1. - s), dk = L / ((1. - s) * (1. - s));
        const double e = 0.5 * kap * kap + delta;
        const double Om = std::sqrt(e * (e + 2. * gam));
        const double b2 = 2. * eta * kap * kap;
        if (b2 > 0. && Om * Om <= b2) return INFINITY;
        if (Om > 700.) continue;
        sum += w * dk * kap * kap * Gfun(Om, b2);
    }
    return 2. / std::sqrt(2. * M_PI) * sum;
}

// ln f by quadrature without underflow: the integrand is multiplied by e^S, S = max(0, delta - max(eta, 0)), and for
// Omega - S large only the k = 1 term of G is kept (in log form). Returns NaN where the result is not reliable (beyond
// the Landau limit, or f_scaled <= 0 from cancellation at eta < 0 with a large gap).
inline double lnf_quad(double delta, double gam, double eta)
{
    const double* xw = gl_nodes();
    const double L = 2.0 / std::sqrt(1. + gam) * std::pow(1. + 4. * eta * eta, 0.25); // smooth in eta (FD across eta = 0)
    const double S = std::max(0., delta - std::max(eta, 0.));
    double sum = 0.;
    for (int i = 0; i < kNodes; i++)
    {
        const double s = xw[2 * i], w = xw[2 * i + 1];
        const double kap = L * s / (1. - s), dk = L / ((1. - s) * (1. - s));
        const double e = 0.5 * kap * kap + delta;
        const double Om = std::sqrt(e * (e + 2. * gam));
        const double b2 = 2. * eta * kap * kap;
        if (b2 > 0. && Om * Om <= b2) return NAN;
        double Gs;
        if (Om > 40.)
        {
            // G ~ e^{-Om} sinh(b)/b (b2 > 0) or e^{-Om} sin(be)/be (b2 < 0), times e^S
            if (b2 > 1e-12)
            {
                const double b = std::sqrt(b2);
                const double lg = S - Om + b + std::log(-std::expm1(-2. * b) / (2. * b));
                Gs = (lg < -745.) ? 0. : std::exp(lg);
            }
            else if (b2 < -1e-12)
            {
                const double be = std::sqrt(-b2);
                const double lg = S - Om;
                Gs = (lg < -745.) ? 0. : std::exp(lg) * std::sin(be) / be;
            }
            else
            {
                const double lg = S - Om;
                Gs = (lg < -745.) ? 0. : std::exp(lg);
            }
        }
        else
            Gs = Gfun(Om, b2) * std::exp(S);
        sum += w * dk * kap * kap * Gs;
    }
    const double fs = 2. / std::sqrt(2. * M_PI) * sum;
    return (fs > 0. && std::isfinite(fs)) ? std::log(fs) - S : NAN;
}

// Phase-mode gap in units of T. Goldstone's theorem makes the phase mode gapless where there is a condensate, while for
// sigma -> 0 the excitations are free particles with gap m - mu_c. We interpolate with the condensate amplitude:
//     delta_lin = (m - mu_c) gfac / T,   gfac = exp(-(sigma / sigma_gap)^4)  (computed by the caller, bogoliubov::gap_factor)
//     delta     = w H(delta_lin / w)   (keeps delta >= 0; H(x) = x for x >= 1, 0 for x <= -1, C^2 in between)
// Why not the earlier (m - mu_c)_+ (or its softplus in mu_c): the gap then closes as a function of mu_c at mu_c = m, which
// makes P_th concave in Pi_psi there (P grows exponentially with mu_c below m and saturates above), and in our stars that
// concavity (~P/T^2) beats the condensate term sigma^2 at the mu_c = m shell, so the Hessian in (T, mu, Pi) is indefinite.
// With the gap linear in mu_c at fixed sigma that concavity is absent; the softplus only acts where gfac (m - mu_c) ~ T,
// i.e. not in equilibrium stars, where mu_c > m only for sigma >> sigma_gap. Delta = m - mu_c (rather than the NR
// (m^2 - mu_c^2)/(2 mu_c)) makes eta - delta = (mu - m)/T exactly at gfac = 1, so the sigma -> 0 limit is the free Bose gas
// independently of mu_c (no spurious, concave mu_c-dependence of O((m - mu_c)^2/T), which made Upsilon < 0 in the tail).
inline double gap_factor(double sigma_sq, double sigma_gap, double* dg_ds2 = nullptr)
{
    if (!(sigma_gap > 0.)) { if (dg_ds2) *dg_ds2 = 0.; return 1.; }
    const double x = std::max(0., sigma_sq) / (sigma_gap * sigma_gap);
    const double g = std::exp(-x * x);
    if (dg_ds2) *dg_ds2 = -g * 2. * x / (sigma_gap * sigma_gap);
    return g;
}
inline double gap_delta(double T, double mu_c, double m, double gap_width, double gfac, double* ddelta_dlin = nullptr)
{
    // delta = w H(delta_lin / w): H = 0 for x <= -1, H = x for x >= 1, and in between H' = S((x+1)/2) with the C^2 quintic
    // smoothstep S(t) = 6t^5 - 15t^4 + 10t^3, i.e. H = 2 (t^6 - 3t^5 + 5t^4/2). Unlike a softplus, delta = delta_lin exactly
    // once delta_lin >= w, so in the gas region P depends on mu only through mu - m (no leak into Upsilon at sigma -> 0).
    const double d_lin = (m - mu_c) * gfac / T;
    if (gap_width <= 0.) { if (ddelta_dlin) *ddelta_dlin = (d_lin > 0.) ? 1. : 0.; return std::max(0., d_lin); }
    const double x = d_lin / gap_width;
    if (x >= 1.) { if (ddelta_dlin) *ddelta_dlin = 1.; return d_lin; }
    if (x <= -1.) { if (ddelta_dlin) *ddelta_dlin = 0.; return 0.; }
    const double t = 0.5 * (x + 1.), t3 = t * t * t;
    if (ddelta_dlin) *ddelta_dlin = t3 * (10. + t * (-15. + 6. * t));
    return gap_width * 2. * t3 * t * (2.5 + t * (-3. + t));
}

// P_th(T, mu, mu_c, a, m) in the scaling form; a = sigma^2 V'', gfac = gap_factor(sigma^2, sigma_gap).
// The kinetic mass is mu_c for the Bogoliubov modes of a condensate and m for free particles; we use
// M = mu_c (1 - gfac) + m gfac in the prefactor (M T / 2pi)^{3/2} T and in gam = a / (2 M T), so that for gfac -> 1 the
// pressure is exactly the free Bose gas (m T/2pi)^{3/2} T Li_{5/2}(e^{(mu - m)/T}), with no dependence on mu_c (Pi_psi).
inline double kinetic_mass(double mu_c, double m, double gfac) { return mu_c * (1. - gfac) + m * gfac; }
inline double pressure(double T, double mu, double mu_c, double a, double m, double gap_width = 1.0, double gfac = 1.0)
{
    if (T <= 0. || mu_c <= 0.) return 0.;
    const double M = kinetic_mass(mu_c, m, gfac);
    const double pre = std::pow(M * T / (2. * M_PI), 1.5) * T;
    return pre * f(gap_delta(T, mu_c, m, gap_width, gfac), std::max(0., a) / (2. * M * T), (mu - mu_c) / T);
}
} // namespace bogoliubov
#endif
