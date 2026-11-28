// Charged (U(1)-gauged) complex scalar field + Maxwell field in spherical symmetry, for the modified-cartoon BSSN/CCZ4
// evolution. Units G = 1, Gaussian Maxwell; the scalar field is normalised as in ComplexScalarField.
//
// Conventions. D_a phi = nabla_a phi + i q A_a phi. The code's static stars have phi = A(r) e^{+i omega t} (see
// BSSNSlice::read_BS_data), i.e. the complex conjugate of the e^{-i omega t} convention of the static charged solver
// (charged_bs/emkg.py, which uses D = nabla - i q A); the sign of the coupling is flipped accordingly so that both
// describe the same physical star, with the same A_t = Phi(r), and q > 0 gives Q > 0.
// Eulerian quantities (n^a the unit normal):
//   Phi_e = -n^a A_a,  A_z = spatial vector potential,  E^i = F^{ia} n_a,  Pi = n^a D_a phi = -2 K_phi.
// Gauge: Lorenz, nabla_a A^a = 0. Static stars are stationary in it (A_z = 0, Phi_e = -A_t / alpha static).
//
// Evolution equations (K = trace of extrinsic curvature, B = 0 in spherical symmetry):
//   d_t phi    = beta d_z phi - 2 alpha K_phi + i q alpha Phi_e phi
//   d_t K_phi  = beta d_z K_phi + i q alpha Phi_e K_phi + alpha K K_phi + (1/2) alpha V'(|phi|^2) phi
//                - (1/2) [ gamma^{ij} d_i alpha D_j phi + alpha gamma^{ij} (D_i D_j phi - Gamma^k_ij D_k phi) ]
//   d_t E^z    = beta d_z E^z - E^z d_z beta + alpha K E^z - 4 pi alpha J^z
//   d_t A_z    = beta d_z A_z + A_z d_z beta - alpha (E_z + d_z Phi_e) - Phi_e d_z alpha
//   d_t Phi_e  = beta d_z Phi_e + alpha K Phi_e - alpha D_i A^i - A^i d_i alpha
// with sources
//   rho_e = -n_a j^a = -2 q Im(phi^* K_phi),   J_z = -q Im(phi^* D_z phi),   j^a = -q Im(phi^* D^a phi),
// and constraint (Gauss)  D_i E^i = 4 pi rho_e. Divergences are taken in conservative form,
//   D_i V^i = (1/W) d_z (W V^z),  W = z^n sqrt(gamma) = z^n sqrt(h_zz h_ww^n) chi^{-(n+1)/2},
// by Spacetime::maxwell_divergence (the pointwise d_z V + n V / z form is unstable at the origin for the Lorenz system).
// The scalar Laplacian is the one of ComplexScalarField with every first derivative d_z phi replaced by D_z phi and
// d_zz phi by D_z D_z phi; in the cartoon angular term, d_w D_w phi -> (D_z phi) / z.
// Electromagnetic stress-energy (E^2 = gamma_zz (E^z)^2):
//   rho_em = E^2 / 8 pi,  j_em = 0,  S_zz = -gamma_zz E^2 / 8 pi,  S_ww = +gamma_ww E^2 / 8 pi.

#include "ChargedScalarField.h"
#include "EvolutionVariables.h"  // for BSSNState and Spacetime
#include <cmath>

EMF operator+(const EMF& a, const EMF& b) { return EMF{a.A_z + b.A_z, a.E_z + b.E_z, a.Phi_e + b.Phi_e}; }
EMF operator-(const EMF& a, const EMF& b) { return EMF{a.A_z - b.A_z, a.E_z - b.E_z, a.Phi_e - b.Phi_e}; }
EMF operator*(const EMF& a, const EMF& b) { return EMF{a.A_z * b.A_z, a.E_z * b.E_z, a.Phi_e * b.Phi_e}; }
EMF operator*(double c, const EMF& a)     { return EMF{c * a.A_z, c * a.E_z, c * a.Phi_e}; }

namespace {
inline double safe_chi(const BSSNState& b) { return (b.chi <= 0.0 ? 1e-300 : b.chi); }
}

CSF ChargedScalarField::cov_d_z(const CSF& csf, const EMF& em, const CSF& d_z_csf) const
{
    // D_z phi = d_z phi + i q A_z phi
    return CSF{d_z_csf.phi_re - q * em.A_z * csf.phi_im,
               d_z_csf.phi_im + q * em.A_z * csf.phi_re, 0.0, 0.0};
}

double ChargedScalarField::rho_scalar(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const
{
    const double chi = safe_chi(bssn);
    const CSF Dphi = cov_d_z(csf, em, d_z_csf);
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = Dphi.phi_re * Dphi.phi_re + Dphi.phi_im * Dphi.phi_im;
    const double Vphi = st ? st->V(mod_phi) : 0.0;
    return 2.0 * K2 + 0.5 * chi * grad2 / bssn.h_zz + 0.5 * Vphi;
}

double ChargedScalarField::rho_em(const EMF& em, const BSSNState& bssn) const
{
    const double E2 = bssn.h_zz * em.E_z * em.E_z / safe_chi(bssn);
    return E2 / (8.0 * M_PI);
}

double ChargedScalarField::rho(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const
{
    return rho_scalar(csf, em, bssn, d_z_csf) + rho_em(em, bssn);
}

double ChargedScalarField::j_z(const CSF& csf, const EMF& em, const BSSNState& /*bssn*/, const CSF& d_z_csf) const
{
    // j_z = 2 Re(K_phi^* D_z phi); the electromagnetic Poynting flux vanishes for B = 0
    const CSF Dphi = cov_d_z(csf, em, d_z_csf);
    return 2.0 * (csf.K_phi_re * Dphi.phi_re + csf.K_phi_im * Dphi.phi_im);
}

double ChargedScalarField::S_zz(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const
{
    const double chi = safe_chi(bssn);
    const CSF Dphi = cov_d_z(csf, em, d_z_csf);
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = Dphi.phi_re * Dphi.phi_re + Dphi.phi_im * Dphi.phi_im;
    const double Vphi = st ? st->V(mod_phi) : 0.0;
    const double gamma_zz = bssn.h_zz / chi;
    return 0.5 * grad2 - 0.5 * gamma_zz * (Vphi - 4.0 * K2) - gamma_zz * rho_em(em, bssn);
}

double ChargedScalarField::S_ww(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const
{
    const double chi = safe_chi(bssn);
    const CSF Dphi = cov_d_z(csf, em, d_z_csf);
    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double K2 = csf.K_phi_re * csf.K_phi_re + csf.K_phi_im * csf.K_phi_im;
    const double grad2 = Dphi.phi_re * Dphi.phi_re + Dphi.phi_im * Dphi.phi_im;
    const double Vphi = st ? st->V(mod_phi) : 0.0;
    const double gamma_ww = bssn.h_ww / chi;
    return -0.5 * bssn.h_ww * (grad2 / bssn.h_zz + (Vphi - 4.0 * K2) / chi) + gamma_ww * rho_em(em, bssn);
}

double ChargedScalarField::charge_density(const CSF& csf) const
{
    // rho_e = -2 q Im(phi^* K_phi)
    return -2.0 * q * (csf.phi_re * csf.K_phi_im - csf.phi_im * csf.K_phi_re);
}

double ChargedScalarField::current_z(const CSF& csf, const EMF& em, const CSF& d_z_csf) const
{
    // J_z = -q Im(phi^* D_z phi)
    const CSF Dphi = cov_d_z(csf, em, d_z_csf);
    return -q * (csf.phi_re * Dphi.phi_im - csf.phi_im * Dphi.phi_re);
}

double ChargedScalarField::gauss(const CSF& csf, double div_E) const
{
    return div_E - 4.0 * M_PI * charge_density(csf);
}

void ChargedScalarField::rhs(const CSF& csf, const EMF& em, const BSSNState& bssn,
                             const CSF& d_z_csf, const CSF& d_zz_csf, const EMF& d_z_em, double div_A,
                             const ChargedGeom& g, CSF& csf_rhs, EMF& em_rhs) const
{
    const double alpha = bssn.alpha;
    const double beta  = bssn.beta;
    const double K     = bssn.K;
    const double chi   = safe_chi(bssn);
    const double h_ZZ  = 1.0 / bssn.h_zz;
    const double h_WW  = 1.0 / bssn.h_ww;
    const double Ddim  = st ? st->D : 4.0;
    const double n     = Ddim - 2.0;
    const double min_z = st ? st->min_z : 0.0;
    const double z     = g.z;

    const double mod_phi = std::sqrt(csf.phi_re * csf.phi_re + csf.phi_im * csf.phi_im);
    const double dVphi   = st ? st->dV(mod_phi) : 0.0;

    const double A   = em.A_z;
    const double Phi = em.Phi_e;
    const double qA  = q * A;

    // gauge-covariant derivatives (complex numbers stored as {re, im})
    const double Dre = d_z_csf.phi_re - qA * csf.phi_im;               // D_z phi
    const double Dim = d_z_csf.phi_im + qA * csf.phi_re;
    // d_z (D_z phi) = d_zz phi + i q (d_z A) phi + i q A d_z phi
    const double dDre = d_zz_csf.phi_re - q * d_z_em.A_z * csf.phi_im - qA * d_z_csf.phi_im;
    const double dDim = d_zz_csf.phi_im + q * d_z_em.A_z * csf.phi_re + qA * d_z_csf.phi_re;
    // D_z D_z phi = d_z(D_z phi) + i q A D_z phi
    const double DDre = dDre - qA * Dim;
    const double DDim = dDim + qA * Dre;

    // covariant second derivatives (conformal Christoffels, as in ComplexScalarField)
    const double cD_zz_re = DDre - g.chris_Zzz * Dre;
    const double cD_zz_im = DDim - g.chris_Zzz * Dim;
    const double D_over_z_re = (z <= min_z) ? dDre : (Dre / z); // (D_z phi)/z, regularised at the origin
    const double D_over_z_im = (z <= min_z) ? dDim : (Dim / z);
    const double cD_ww_re = D_over_z_re - g.chris_Zww * Dre;
    const double cD_ww_im = D_over_z_im - g.chris_Zww * Dim;

    const double qaP = q * alpha * Phi;

    // d_t phi = beta d_z phi - 2 alpha K_phi + i q alpha Phi_e phi
    csf_rhs.phi_re = beta * d_z_csf.phi_re - 2.0 * alpha * csf.K_phi_re - qaP * csf.phi_im;
    csf_rhs.phi_im = beta * d_z_csf.phi_im - 2.0 * alpha * csf.K_phi_im + qaP * csf.phi_re;

    csf_rhs.K_phi_re = beta * d_z_csf.K_phi_re - qaP * csf.K_phi_im
                     + alpha * K * csf.K_phi_re
                     + 0.5 * alpha * csf.phi_re * dVphi
                     - 0.5 * chi * ( h_ZZ * g.d_z_alpha * Dre + alpha * ( h_ZZ * cD_zz_re + n * h_WW * cD_ww_re ) )
                     + 0.25 * (Ddim - 3.0) * alpha * h_ZZ * g.d_z_chi * Dre;

    csf_rhs.K_phi_im = beta * d_z_csf.K_phi_im + qaP * csf.K_phi_re
                     + alpha * K * csf.K_phi_im
                     + 0.5 * alpha * csf.phi_im * dVphi
                     - 0.5 * chi * ( h_ZZ * g.d_z_alpha * Dim + alpha * ( h_ZZ * cD_zz_im + n * h_WW * cD_ww_im ) )
                     + 0.25 * (Ddim - 3.0) * alpha * h_ZZ * g.d_z_chi * Dim;

    // Maxwell
    const double J_z = -q * (csf.phi_re * Dim - csf.phi_im * Dre);   // covariant current
    const double J_Z = chi * h_ZZ * J_z;                              // contravariant
    const double E_low = bssn.h_zz * em.E_z / chi;                    // covariant E_z

    em_rhs.E_z = beta * d_z_em.E_z - em.E_z * g.d_z_beta + alpha * K * em.E_z - 4.0 * M_PI * alpha * J_Z;

    em_rhs.A_z = beta * d_z_em.A_z + A * g.d_z_beta - alpha * (E_low + d_z_em.Phi_e) - Phi * g.d_z_alpha;

    // Lorenz gauge: d_t Phi_e = beta d_z Phi_e + alpha K Phi_e - alpha D_i A^i - A^i d_i alpha
    const double A_up = chi * h_ZZ * A;
    em_rhs.Phi_e = beta * d_z_em.Phi_e + alpha * K * Phi - alpha * div_A - A_up * g.d_z_alpha;
}
