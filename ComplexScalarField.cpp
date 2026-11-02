#include "ComplexScalarField.h"
#include "EvolutionVariables.h"  // for BSSNState and Spacetime
#include <cmath>

// Elementwise arithmetic for CSF
CSF operator+(const CSF& a, const CSF& b) {
    return CSF{a.phi_re + b.phi_re, a.phi_im + b.phi_im,
               a.K_phi_re + b.K_phi_re, a.K_phi_im + b.K_phi_im};
}

CSF operator-(const CSF& a, const CSF& b) {
    return CSF{a.phi_re - b.phi_re, a.phi_im - b.phi_im,
               a.K_phi_re - b.K_phi_re, a.K_phi_im - b.K_phi_im};
}

CSF operator*(const CSF& a, const CSF& b) {
    return CSF{a.phi_re * b.phi_re, a.phi_im * b.phi_im,
               a.K_phi_re * b.K_phi_re, a.K_phi_im * b.K_phi_im};
}

CSF operator*(double c, const CSF& a) {
    return CSF{c * a.phi_re, c * a.phi_im, c * a.K_phi_re, c * a.K_phi_im};
}


// l-boson stars: sum over the multiplet of the angular gradients. With Phi_m = sqrt(4 pi/(2 ell+1)) f Y_{ell m},
// sum_m |Phi_m|^2 = |f|^2 and, in the cartoon w-direction at w = 0, sum_m |d_w Phi_m|^2 = ell(ell+D-3)|f|^2 / (n z^2)
// per transverse direction (n = D - 2). This helper returns ell(ell+D-3)|f|^2/z^2.
double ComplexScalarField::angular_grad2(const CSF& csf, const CSF& d_z_csf, double z) const {
    if (!st || st->ell == 0 || z < 0.0)
        return 0.0;
    const double LL = st->ell_LL();
    if (z <= st->min_z) // f ~ z^ell: f/z -> d_z f at the origin (non-zero only for ell = 1)
        return LL * (d_z_csf.phi_re * d_z_csf.phi_re + d_z_csf.phi_im * d_z_csf.phi_im);
    return LL * (csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im) / (z * z);
}

// Energy density rho = n^a n^b T_ab
double ComplexScalarField::rho(const CSF& csf, const BSSNState& bssn,
                               const CSF& d_z_csf, const CSF& /*d_zz_csf*/, double z) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_ZZ = 1.0 / bssn.h_zz;
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = d_z_csf.phi_re * d_z_csf.phi_re + d_z_csf.phi_im * d_z_csf.phi_im;
    double Vphi = st ? st->V(mod_phi) : 0.0;
    double r = 2.0 * K2 + 0.5 * h_ZZ * chi * grad2 + 0.5 * Vphi;
    if (st && st->ell > 0) // + (1/2) gamma^{AB} sum_m d_A Phi_m^* d_B Phi_m
        r += 0.5 * chi * angular_grad2(csf, d_z_csf, z) / bssn.h_ww;
    return r;
}

// Momentum density j_z = -n^a gamma^b{}_z T_ab
double ComplexScalarField::j_z(const CSF& csf, const BSSNState& /*bssn*/,
                               const CSF& d_z_csf, const CSF& /*d_zz_csf*/) const {
    return 2.0 * (csf.K_phi_re * d_z_csf.phi_re + csf.K_phi_im * d_z_csf.phi_im);
}

// Stress S_zz = gamma^a{}_z gamma^b{}_z T_ab
double ComplexScalarField::S_zz(const CSF& csf, const BSSNState& bssn,
                                const CSF& d_z_csf, const CSF& /*d_zz_csf*/, double z) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_zz = bssn.h_zz;
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = d_z_csf.phi_re * d_z_csf.phi_re + d_z_csf.phi_im * d_z_csf.phi_im;
    double Vphi = st ? st->V(mod_phi) : 0.0;
    double s = 0.5 * grad2 - 0.5 * h_zz * (Vphi - 4.0 * K2) / chi;
    if (st && st->ell > 0) // - (1/2) gamma_zz gamma^{AB} sum_m d_A Phi_m^* d_B Phi_m
        s -= 0.5 * h_zz * angular_grad2(csf, d_z_csf, z) / bssn.h_ww;
    return s;
}

// Stress S_ww (angular trace component)
double ComplexScalarField::S_ww(const CSF& csf, const BSSNState& bssn,
                                const CSF& d_z_csf, const CSF& /*d_zz_csf*/, double z) const {
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_ww = bssn.h_ww;
    const double h_ZZ = 1.0 / bssn.h_zz;
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = d_z_csf.phi_re * d_z_csf.phi_re + d_z_csf.phi_im * d_z_csf.phi_im;
    double Vphi = st ? st->V(mod_phi) : 0.0;
    double s = -0.5 * h_ww * (h_ZZ * grad2 + (Vphi - 4.0 * K2) / chi);
    if (st && st->ell > 0) // sum_m |d_w Phi_m|^2 - (1/2) gamma_ww gamma^{AB} (...) = ell(ell+D-3)|f|^2/z^2 (1/n - 1/2); zero in D = 4
        s += angular_grad2(csf, d_z_csf, z) * (1.0 / (st->D - 2.0) - 0.5);
    return s;
}

// Full RHS using available local data and field derivatives.
CSF ComplexScalarField::rhs_CSF(const CSF& csf, const BSSNState& bssn,
                                const CSF& d_z_csf, const CSF& d_zz_csf,
                                double d_z_alpha, double d_z_chi,
                                double chris_Zzz, double chris_Zww,
                                double z) const {
    CSF rhs{};

    const double alpha = bssn.alpha;
    const double beta  = bssn.beta;
    const double K     = bssn.K;
    const double chi   = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_ZZ  = 1.0 / bssn.h_zz;
    const double h_WW  = 1.0 / bssn.h_ww;
    const double Ddim  = st ? st->D : 4.0; // use st->D when available
    const double n     = Ddim - 2.0;

    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double dVphi   = st ? st->dV(mod_phi) : 0.0;

    rhs.phi_re   = beta * d_z_csf.phi_re - 2.0 * alpha * csf.K_phi_re;
    rhs.phi_im   = beta * d_z_csf.phi_im - 2.0 * alpha * csf.K_phi_im;

    // Covariant 2nd derivatives for a scalar: D_i D_j phi = d_ij phi - Gamma^k_ij d_k phi
    const double cD_zz_phi_re = d_zz_csf.phi_re - chris_Zzz * d_z_csf.phi_re;
    const double cD_zz_phi_im = d_zz_csf.phi_im - chris_Zzz * d_z_csf.phi_im;
    
    //helper: (d_phi/dz)/z with regularization at the origin
    const double d_phi_re_z_over_z = (z <= (st ? st->min_z : 0.0)) ? d_zz_csf.phi_re : (d_z_csf.phi_re / z);
    const double d_phi_im_z_over_z = (z <= (st ? st->min_z : 0.0)) ? d_zz_csf.phi_im : (d_z_csf.phi_im / z);

    const double cD_ww_phi_re = d_phi_re_z_over_z - chris_Zww * d_z_csf.phi_re;
    const double cD_ww_phi_im = d_phi_im_z_over_z - chris_Zww * d_z_csf.phi_im;

    // K_phi evolution including advection, trace coupling, potential, and geometric terms
    rhs.K_phi_re = beta * d_z_csf.K_phi_re
                                + alpha * K * csf.K_phi_re
                                + 0.5 * alpha * csf.phi_re * dVphi
                                - 0.5 * chi * (
                                            h_ZZ * d_z_alpha * d_z_csf.phi_re
                                        + alpha * ( h_ZZ * cD_zz_phi_re + n * h_WW * cD_ww_phi_re )
                                    )
                                + 0.25 * (Ddim - 3.0) * alpha * h_ZZ * d_z_chi * d_z_csf.phi_re;

    rhs.K_phi_im = beta * d_z_csf.K_phi_im
                                + alpha * K * csf.K_phi_im
                                + 0.5 * alpha * csf.phi_im * dVphi
                                - 0.5 * chi * (
                                            h_ZZ * d_z_alpha * d_z_csf.phi_im
                                        + alpha * ( h_ZZ * cD_zz_phi_im + n * h_WW * cD_ww_phi_im )
                                    )
                                + 0.25 * (Ddim - 3.0) * alpha * h_ZZ * d_z_chi * d_z_csf.phi_im;

    // l-boson stars: angular part of the Laplacian, Delta_S Y_lm = -ell(ell+D-3) Y_lm, gives a centrifugal term
    // -(1/2) alpha * (-chi h^ww ell(ell+D-3) f / z^2). f ~ z^ell, so f and K_phi vanish at the origin for ell > 0:
    // there the RHS is pinned to zero (see also slice_rhs, which re-pins after dissipation).
    if (st && st->ell > 0)
    {
        if (z <= st->min_z)
            return CSF{0.0, 0.0, 0.0, 0.0};
        const double cent = 0.5 * alpha * chi * h_WW * st->ell_LL() / (z * z);
        rhs.K_phi_re += cent * csf.phi_re;
        rhs.K_phi_im += cent * csf.phi_im;
    }

    return rhs;
}
