#include "ComplexVectorField.h"
#include "EvolutionVariables.h"  // for BSSNState and Spacetime
#include <cmath>

// Elementwise arithmetic for CVF
CVF operator+(const CVF& a, const CVF& b) {
    return CVF{a.X_re + b.X_re, a.X_im + b.X_im, a.E_re + b.E_re, a.E_im + b.E_im,
               a.Xphi_re + b.Xphi_re, a.Xphi_im + b.Xphi_im, a.Z_re + b.Z_re, a.Z_im + b.Z_im};
}

CVF operator-(const CVF& a, const CVF& b) {
    return CVF{a.X_re - b.X_re, a.X_im - b.X_im, a.E_re - b.E_re, a.E_im - b.E_im,
               a.Xphi_re - b.Xphi_re, a.Xphi_im - b.Xphi_im, a.Z_re - b.Z_re, a.Z_im - b.Z_im};
}

CVF operator*(const CVF& a, const CVF& b) {
    return CVF{a.X_re * b.X_re, a.X_im * b.X_im, a.E_re * b.E_re, a.E_im * b.E_im,
               a.Xphi_re * b.Xphi_re, a.Xphi_im * b.Xphi_im, a.Z_re * b.Z_re, a.Z_im * b.Z_im};
}

CVF operator*(double c, const CVF& a) {
    return CVF{c * a.X_re, c * a.X_im, c * a.E_re, c * a.E_im,
               c * a.Xphi_re, c * a.Xphi_im, c * a.Z_re, c * a.Z_im};
}

// D_i V^i = d_z V^z + (D-2) V^z / z + V^z d_z ln sqrt(gamma), with
// d_z ln sqrt(gamma) = h^zz d_z h_zz / 2 + (D-2) h^ww d_z h_ww / 2 - (D-1) d_z chi / (2 chi).
// The (D-2) V^z / z term follows from d_a V^b = delta_a^b V^z / z (Cook et al. 2016, eq. A.4); at the origin
// V^z / z is replaced by its limit, supplied by the caller.
double ComplexVectorField::divergence(double V, double dV, double V_over_z_origin, const BSSNState& bssn,
                                      double d_z_chi, double d_z_h_zz, double d_z_h_ww, double z) const
{
    const double D = st->D;
    const double n = D - 2.0;
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double d_ln_sqrt_gamma = 0.5 * d_z_h_zz / bssn.h_zz + 0.5 * n * d_z_h_ww / bssn.h_ww - 0.5 * (D - 1.0) * d_z_chi / chi;
    const double V_over_z = (z <= st->min_z) ? V_over_z_origin : V / z;

    return dV + n * V_over_z + V * d_ln_sqrt_gamma;
}

// Energy density: rho = 1/2 E_i E^i + 1/2 mu^2 (X_phi^2 + X_i X^i), summed over real and imaginary parts
double ComplexVectorField::rho(const CVF& cvf, const BSSNState& bssn) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double g_zz = bssn.h_zz / chi;     // physical gamma_zz
    const double g_ZZ = chi / bssn.h_zz;     // physical gamma^zz
    const double mu2 = st->mu * st->mu;

    const double E2 = g_zz * (cvf.E_re * cvf.E_re + cvf.E_im * cvf.E_im);
    const double X2 = g_ZZ * (cvf.X_re * cvf.X_re + cvf.X_im * cvf.X_im);
    const double P2 = cvf.Xphi_re * cvf.Xphi_re + cvf.Xphi_im * cvf.Xphi_im;

    return 0.5 * E2 + 0.5 * mu2 * (P2 + X2);
}

// Momentum density: j_i = mu^2 X_phi X_i (the Poynting-like E x B part vanishes in spherical symmetry)
double ComplexVectorField::j_z(const CVF& cvf, const BSSNState& /*bssn*/) const {
    const double mu2 = st->mu * st->mu;
    return mu2 * (cvf.Xphi_re * cvf.X_re + cvf.Xphi_im * cvf.X_im);
}

// Stress: S_ij = -E_i E_j + 1/2 gamma_ij E^2 + mu^2 (X_i X_j - 1/2 gamma_ij (X_k X^k - X_phi^2))
double ComplexVectorField::S_zz(const CVF& cvf, const BSSNState& bssn) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double g_zz = bssn.h_zz / chi;
    const double g_ZZ = chi / bssn.h_zz;
    const double mu2 = st->mu * st->mu;

    const double E_up2 = cvf.E_re * cvf.E_re + cvf.E_im * cvf.E_im;   // (E^z)^2
    const double X_dn2 = cvf.X_re * cvf.X_re + cvf.X_im * cvf.X_im;   // (X_z)^2
    const double P2 = cvf.Xphi_re * cvf.Xphi_re + cvf.Xphi_im * cvf.Xphi_im;

    // -E_z E_z + 1/2 g_zz (g_zz E^z E^z) = -1/2 g_zz^2 (E^z)^2
    return -0.5 * g_zz * g_zz * E_up2 + mu2 * (X_dn2 - 0.5 * g_zz * (g_ZZ * X_dn2 - P2));
}

double ComplexVectorField::S_ww(const CVF& cvf, const BSSNState& bssn) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double g_zz = bssn.h_zz / chi;
    const double g_ZZ = chi / bssn.h_zz;
    const double g_ww = bssn.h_ww / chi;
    const double mu2 = st->mu * st->mu;

    const double E2 = g_zz * (cvf.E_re * cvf.E_re + cvf.E_im * cvf.E_im);
    const double X2 = g_ZZ * (cvf.X_re * cvf.X_re + cvf.X_im * cvf.X_im);
    const double P2 = cvf.Xphi_re * cvf.Xphi_re + cvf.Xphi_im * cvf.Xphi_im;

    return 0.5 * g_ww * E2 - 0.5 * mu2 * g_ww * (X2 - P2);
}

double ComplexVectorField::charge_density(const CVF& cvf) const {
    return cvf.E_im * cvf.X_re - cvf.E_re * cvf.X_im;
}

void ComplexVectorField::gauss_constraint(const CVF& cvf, const BSSNState& bssn, const CVF& d_z_cvf,
                                          double d_z_chi, double d_z_h_zz, double d_z_h_ww, double z,
                                          double& C_re, double& C_im) const
{
    const double mu2 = st->mu * st->mu;
    C_re = divergence(cvf.E_re, d_z_cvf.E_re, d_z_cvf.E_re, bssn, d_z_chi, d_z_h_zz, d_z_h_ww, z) + mu2 * cvf.Xphi_re;
    C_im = divergence(cvf.E_im, d_z_cvf.E_im, d_z_cvf.E_im, bssn, d_z_chi, d_z_h_zz, d_z_h_ww, z) + mu2 * cvf.Xphi_im;
}

// Full RHS for the Proca variables at a gridpoint
CVF ComplexVectorField::rhs_CVF(const CVF& cvf, const BSSNState& bssn, const CVF& d_z_cvf,
                                double d_z_alpha, double d_z_beta, double d_z_chi,
                                double d_z_h_zz, double d_z_h_ww, double z) const
{
    const double alpha = bssn.alpha;
    const double beta  = bssn.beta;
    const double K     = bssn.K;
    const double chi   = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double g_zz  = bssn.h_zz / chi;                 // physical gamma_zz
    const double g_ZZ  = chi / bssn.h_zz;                 // physical gamma^zz
    const double d_g_ZZ = d_z_chi / bssn.h_zz - chi * d_z_h_zz / (bssn.h_zz * bssn.h_zz);
    const double mu2   = st->mu * st->mu;
    const double kappa = st->proca_damping;

    // one real component of the field; the complex field is two independent copies
    auto rhs_component = [&](double X, double E, double P, double Z,
                             double dX, double dE, double dP, double dZ,
                             double& rX, double& rE, double& rP, double& rZ)
    {
        const double X_up = g_ZZ * X;                      // X^z
        const double dX_up = d_g_ZZ * X + g_ZZ * dX;       // d_z X^z; X^z / z -> d_z X^z at the origin (X odd)

        const double div_X = divergence(X_up, dX_up, dX_up, bssn, d_z_chi, d_z_h_zz, d_z_h_ww, z);
        const double div_E = divergence(E, dE, dE, bssn, d_z_chi, d_z_h_zz, d_z_h_ww, z);

        rX = beta * dX + X * d_z_beta - alpha * (g_zz * E + dP) - P * d_z_alpha;
        rE = beta * dE - E * d_z_beta + alpha * K * E + alpha * mu2 * X_up + alpha * g_ZZ * dZ;
        rP = beta * dP + alpha * K * P - alpha * div_X - X_up * d_z_alpha - alpha * Z;
        rZ = beta * dZ + alpha * (div_E + mu2 * P - kappa * Z);
    };

    CVF rhs{};
    rhs_component(cvf.X_re, cvf.E_re, cvf.Xphi_re, cvf.Z_re,
                  d_z_cvf.X_re, d_z_cvf.E_re, d_z_cvf.Xphi_re, d_z_cvf.Z_re,
                  rhs.X_re, rhs.E_re, rhs.Xphi_re, rhs.Z_re);
    rhs_component(cvf.X_im, cvf.E_im, cvf.Xphi_im, cvf.Z_im,
                  d_z_cvf.X_im, d_z_cvf.E_im, d_z_cvf.Xphi_im, d_z_cvf.Z_im,
                  rhs.X_im, rhs.E_im, rhs.Xphi_im, rhs.Z_im);
    return rhs;
}
