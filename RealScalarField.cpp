#include "RealScalarField.h"
#include "EvolutionVariables.h"  // for BSSNState and Spacetime

// NOTE: This is a scaffolding implementation. Bodies are placeholders
// and will be filled with correct physics in subsequent iterations.

// Elementwise arithmetic for RSF
RSF operator+(const RSF& a, const RSF& b) {
    return RSF{a.psi + b.psi, a.K_psi + b.K_psi};
}

RSF operator-(const RSF& a, const RSF& b) {
    return RSF{a.psi - b.psi, a.K_psi - b.K_psi};
}

RSF operator*(const RSF& a, const RSF& b) {
    return RSF{a.psi * b.psi, a.K_psi * b.K_psi};
}

RSF operator*(double c, const RSF& a) {
    return RSF{c * a.psi, c * a.K_psi};
}

RSF RealScalarField::rhs_RSF(const RSF& rsf, const BSSNState& bssn,
                                                         const RSF& d_z_rsf, const RSF& d_zz_rsf,
                                                         double d_z_alpha, double d_z_chi,
                                                         double chris_Zzz, double chris_Zww,
                                                         double z) const {
        RSF rhs{};

        const double alpha = bssn.alpha;
        const double beta  = bssn.beta;
        const double K     = bssn.K;
        const double chi   = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
        const double h_ZZ  = 1.0 / bssn.h_zz;
        const double h_WW  = 1.0 / bssn.h_ww;
        const double Ddim  = st ? st->D : 4.0; // use spacetime D when available
        const double n     = Ddim - 2.0;

        // Advective + algebraic parts
        rhs.psi   = beta * d_z_rsf.psi - 2.0 * alpha * rsf.K_psi;

        // Covariant second derivatives for a scalar: D_i D_j psi = d_ij psi - Gamma^k_ij d_k psi
        const double cD_zz_psi = d_zz_rsf.psi - chris_Zzz * d_z_rsf.psi;

        // helper: (d_psi/dz)/z with regularization at the origin (use second derivative at z~0)
        const double d_psi_z_over_z = (z <= (st ? st->min_z : 0.0)) ? d_zz_rsf.psi : (d_z_rsf.psi / z);
        const double cD_ww_psi = d_psi_z_over_z - chris_Zww * d_z_rsf.psi;

        // K_psi evolution including advection, trace coupling, and geometric terms (no potential for massless field)
        rhs.K_psi = beta * d_z_rsf.K_psi
                                + alpha * K * rsf.K_psi
                                - 0.5 * chi * (
                                                h_ZZ * d_z_alpha * d_z_rsf.psi
                                            + alpha * ( h_ZZ * cD_zz_psi + n * h_WW * cD_ww_psi )
                                    )
                                + 0.25 * (Ddim - 3.0) * alpha * h_ZZ * d_z_chi * d_z_rsf.psi;

        return rhs;
}

double RealScalarField::rho(const RSF& rsf, const BSSNState& bssn,
                            const RSF& d_z_rsf, const RSF& /*d_zz_rsf*/) const {
    // Mirror ComplexScalarField::rho but for a single real component and no potential
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_ZZ = 1.0 / bssn.h_zz;
    const double K2 = rsf.K_psi * rsf.K_psi;
    const double grad2 = d_z_rsf.psi * d_z_rsf.psi;
    return 2.0 * K2 + 0.5 * h_ZZ * chi * grad2;
}

double RealScalarField::j_z(const RSF& rsf, const BSSNState& /*bssn*/,
                             const RSF& d_z_rsf, const RSF& /*d_zz_rsf*/) const {
    // Mirror ComplexScalarField::j_z coefficients
    return 2.0 * (rsf.K_psi * d_z_rsf.psi);
}

double RealScalarField::S_zz(const RSF& rsf, const BSSNState& bssn,
                              const RSF& d_z_rsf, const RSF& /*d_zz_rsf*/) const {
    // Mirror ComplexScalarField::S_zz with V=0 and single component
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_zz = bssn.h_zz;
    const double K2 = rsf.K_psi * rsf.K_psi;
    const double grad2 = d_z_rsf.psi * d_z_rsf.psi;
    return 0.5 * grad2 + 2.0 * h_zz * K2 / chi;
}

double RealScalarField::S_ww(const RSF& rsf, const BSSNState& bssn,
                              const RSF& d_z_rsf, const RSF& /*d_zz_rsf*/) const {
    // Mirror ComplexScalarField::S_ww with V=0 and single component
    const double chi = (bssn.chi <= 0.0 ? 1e-300 : bssn.chi);
    const double h_ww = bssn.h_ww;
    const double h_ZZ = 1.0 / bssn.h_zz;
    const double K2 = rsf.K_psi * rsf.K_psi;
    const double grad2 = d_z_rsf.psi * d_z_rsf.psi;
    return -0.5 * h_ww * (h_ZZ * grad2 - 4.0 * K2 / chi);
}
