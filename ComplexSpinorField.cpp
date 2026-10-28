#include "ComplexSpinorField.h"
#include "EvolutionVariables.h"  // for BSSNState and Spacetime
#include "DiracPotential.h"
#include <cmath>

CSpF operator+(const CSpF& a, const CSpF& b) { return CSpF{a.u_re + b.u_re, a.u_im + b.u_im, a.v_re + b.v_re, a.v_im + b.v_im}; }
CSpF operator-(const CSpF& a, const CSpF& b) { return CSpF{a.u_re - b.u_re, a.u_im - b.u_im, a.v_re - b.v_re, a.v_im - b.v_im}; }
CSpF operator*(const CSpF& a, const CSpF& b) { return CSpF{a.u_re * b.u_re, a.u_im * b.u_im, a.v_re * b.v_re, a.v_im * b.v_im}; }
CSpF operator*(double c, const CSpF& a) { return CSpF{c * a.u_re, c * a.u_im, c * a.v_re, c * a.v_im}; }

// scalar density S = sum psibar psi = (|u|^2 - |v|^2) / (2 pi), which sets mu_eff = mu + U'(S) and the potential terms
static inline double scalar_density(const CSpF& s)
{
    return (s.u_re * s.u_re + s.u_im * s.u_im - s.v_re * s.v_re - s.v_im * s.v_im) / (2. * M_PI);
}

// the self-interaction potential of the spacetime's Dirac field (DiracPotential.h)
DiracPotential ComplexSpinorField::potential() const
{
    return DiracPotential{st->mu, st->lambda, st->solitonic, st->sigma};
}

CSpF ComplexSpinorField::rhs_CSpF(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s,
                                  double d_z_alpha, double d_z_chi, double d_z_h_ww, double z) const
{
    const double alpha = bssn.alpha;
    const double beta = bssn.beta;
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double mu = potential().mu_eff(scalar_density(s)); // mu_eff

    const double cz = std::sqrt(chi / bssn.h_zz);   // 1 / sqrt(g_zz)
    const double cw = std::sqrt(chi / bssn.h_ww);   // 1 / sqrt(g_ww)
    const double L = 0.5 * (d_z_h_ww / bssn.h_ww - d_z_chi / chi + d_z_alpha / alpha);

    // v / z (v odd) and (cz - cw) u / z (u even, h_zz - h_ww = O(z^2)), regularised at the origin
    const bool origin = (z <= st->min_z);
    const double vz_re = origin ? d_z_s.v_re : s.v_re / z;
    const double vz_im = origin ? d_z_s.v_im : s.v_im / z;
    const double cuz = origin ? 0. : (cz - cw) / z;

    // P = cz (v' + L v) + (cz + cw) v / z,   Q = cz (u' + L u) + (cz - cw) u / z
    const double P_re = cz * (d_z_s.v_re + L * s.v_re) + (cz + cw) * vz_re;
    const double P_im = cz * (d_z_s.v_im + L * s.v_im) + (cz + cw) * vz_im;
    const double Q_re = cz * (d_z_s.u_re + L * s.u_re) + cuz * s.u_re;
    const double Q_im = cz * (d_z_s.u_im + L * s.u_im) + cuz * s.u_im;

    const double hK = 0.5 * alpha * bssn.K;
    const double am = alpha * mu;

    CSpF rhs{};
    // d_t u = beta u' + (alpha K/2) u - i alpha mu u + i alpha P
    rhs.u_re = beta * d_z_s.u_re + hK * s.u_re + am * s.u_im - alpha * P_im;
    rhs.u_im = beta * d_z_s.u_im + hK * s.u_im - am * s.u_re + alpha * P_re;
    // d_t v = beta v' + (alpha K/2) v + i alpha mu v - i alpha Q
    rhs.v_re = beta * d_z_s.v_re + hK * s.v_re - am * s.v_im + alpha * Q_im;
    rhs.v_im = beta * d_z_s.v_im + hK * s.v_im + am * s.v_re - alpha * Q_re;
    return rhs;
}

double ComplexSpinorField::number_density(const CSpF& s) const
{
    return (s.u_re * s.u_re + s.u_im * s.u_im + s.v_re * s.v_re + s.v_im * s.v_im) / (2. * M_PI);
}

// S_rr^ = Re(vbar u' - u vbar') / (2 pi sqrt(g_zz)) + P(S)
double ComplexSpinorField::S_rr_hat(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s) const
{
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double cz = std::sqrt(chi / bssn.h_zz);
    const double re = (s.v_re * d_z_s.u_re + s.v_im * d_z_s.u_im) - (s.u_re * d_z_s.v_re + s.u_im * d_z_s.v_im);
    return cz * re / (2. * M_PI) + potential().P(scalar_density(s));
}

// S_thth^ = -Re(u vbar) / (2 pi z sqrt(g_ww)) + P(S); v / z -> v' at the origin
double ComplexSpinorField::S_tt_hat(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const
{
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double cw = std::sqrt(chi / bssn.h_ww);
    const bool origin = (z <= st->min_z);
    const double vz_re = origin ? d_z_s.v_re : s.v_re / z;
    const double vz_im = origin ? d_z_s.v_im : s.v_im / z;
    return -cw * (s.u_re * vz_re + s.u_im * vz_im) / (2. * M_PI) + potential().P(scalar_density(s));
}

double ComplexSpinorField::rho(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const
{
    // rho = mu S + U(S) + kinetic stresses; S_rr_hat + 2 S_tt_hat include 3 P(S), which is subtracted
    const double S = scalar_density(s);
    const DiracPotential pot = potential();
    return st->mu * S + pot.U(S) - 3. * pot.P(S) + S_rr_hat(s, bssn, d_z_s) + 2. * S_tt_hat(s, bssn, d_z_s, z);
}

// j_z = sqrt(g_zz) j_r^ = Im(ubar u' + vbar v') / (2 pi)
double ComplexSpinorField::j_z(const CSpF& s, const BSSNState& /*bssn*/, const CSpF& d_z_s) const
{
    const double im = (s.u_re * d_z_s.u_im - s.u_im * d_z_s.u_re) + (s.v_re * d_z_s.v_im - s.v_im * d_z_s.v_re);
    return im / (2. * M_PI);
}

double ComplexSpinorField::S_zz(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s) const
{
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    return (bssn.h_zz / chi) * S_rr_hat(s, bssn, d_z_s);
}

double ComplexSpinorField::S_ww(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const
{
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    return (bssn.h_ww / chi) * S_tt_hat(s, bssn, d_z_s, z);
}
