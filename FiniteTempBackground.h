#ifndef FINITE_TEMP_BACKGROUND_H_
#define FINITE_TEMP_BACKGROUND_H_

// Read-in plumbing for finite-temperature boson-star backgrounds produced by the relaxation solver in
// finite_temp/ (ported from Finite_Temp_Relaxer). The relaxer works on the compactified coordinate
// rbar = r / (1 + r) in polar-areal gauge,
//     ds^2 = -alpha^2 dt^2 + Xs^2 dr^2 + r^2 dOmega^2,   eta = (d_r sigma) / Xs,
// which is exactly SBSE's BosonStar convention (A = sigma, X = Xs, phi = log(alpha), eta = eta, and the
// same potential V = mu^2 A^2 + lambda A^4 (+ solitonic factor), so lambda is passed through unchanged;
// note lambda_SBSE = lambda_relaxer = lambda_draft / 2).
//
// Sign conventions: the relaxer and the finite-T derivation (finiteT_BSSN_equations.tex) use the draft's
// phi = sigma e^{-i omega t}, so mu = omega / alpha > 0 and n is odd in mu. SBSE's standard initial data use
// phi ~ e^{+i omega t}. All thermodynamic quantities stored here are in the draft convention; with
// ft_draft_phase = 1 (default) read_BS_data also sets K_phi_im = +omega phi_re / (2 alpha) so the scalar field
// matches them.
//
// Only D = 4 is supported (the relaxer equations are four-dimensional).


// Static background thermodynamics at one radius (draft sign convention). POD so it can be broadcast as doubles.
struct FiniteTempPoint
{
    double T = 0.;       // local temperature, T_inf / alpha (Tolman)
    double mu = 0.;      // local chemical potential, omega / alpha (Klein)
    double P = 0.;       // total pressure P(T, mu, sigma^2, X, Xi) = L + P_thermal
    double s = 0.;       // entropy density
    double n = 0.;       // normal-component number density
    double Upsilon = 0.; // superfluid coefficient, sigma^2 + delta Upsilon
    double F = 0.;       // 2 dP/d(sigma^2)
    double w = 0.;       // normal enthalpy T s + mu n
    double rho = 0.;     // energy density n^a n^b T_ab = mu^2 Upsilon + w - P
};
constexpr int FT_POINT_NDOUBLES = sizeof(FiniteTempPoint) / sizeof(double);

// Initial values of the extra evolved matter variables of the finite-T formulation (tex, Sec. "Choice of evolved
// variables") for a static background: v = 0 (W = 1), Pi_psi = -mu, K_sigma = 0. Not evolved yet.
struct FiniteTempConserved
{
    double K_sigma = 0.; // -Pi_sigma / 2
    double N_S = 0.;     // -Upsilon Pi_psi  (superfluid charge density)
    double N_N = 0.;     // W n            (normal charge density)
    double Sigma = 0.;   // W s            (entropy density)
    double P_z = 0.;     // w W^2 v_z      (normal momentum density)
};
inline FiniteTempConserved ft_static_conserved(const FiniteTempPoint& p)
{
    FiniteTempConserved c;
    c.N_S = p.Upsilon * p.mu;
    c.N_N = p.n;
    c.Sigma = p.s;
    return c;
}

#endif /* FINITE_TEMP_BACKGROUND_H_ */
