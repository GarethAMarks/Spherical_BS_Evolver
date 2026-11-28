#ifndef CHARGED_SCALAR_FIELD_H
#define CHARGED_SCALAR_FIELD_H

#include <cstddef>
#include "ComplexScalarField.h" // CSF holds the (gauge-covariant) scalar variables

// Forward declarations to avoid heavy includes in this header
struct BSSNState;   // defined in EvolutionVariables.h
class Spacetime;    // defined in EvolutionVariables.h/.cpp

// Maxwell field in spherical symmetry (no magnetic field), on a slice. Gaussian units, G = 1: the action is
//   S = int sqrt(-g) [ R / 16 pi - (1/2) ( g^{ab} (D_a phi)^* D_b phi + V(|phi|) ) - F_{ab} F^{ab} / 16 pi ],
// so the scalar field keeps the normalisation of ComplexScalarField. Conventions (see ChargedScalarField.cpp):
//   D_a = nabla_a + i q A_a  (the code's stars rotate as phi ~ e^{+i omega t}, so q > 0 gives Q > 0),
//   E^i = F^{ia} n_a  (Eulerian electric field), Phi_e = -n^a A_a (Eulerian scalar potential).
struct EMF {
    double A_z;    // radial component of the spatial vector potential, covariant (odd across z = 0)
    double E_z;    // radial electric field, contravariant E^z (odd across z = 0)
    double Phi_e;  // scalar potential -n^a A_a (even); evolved with the Lorenz gauge condition
};

// Elementwise arithmetic for EMF (mirroring CSF)
EMF operator+(const EMF& a, const EMF& b);
EMF operator-(const EMF& a, const EMF& b);
EMF operator*(const EMF& a, const EMF& b);       // termwise multiply
EMF operator*(double c, const EMF& a);           // scalar multiply (left)
inline EMF operator/(const EMF& a, double c) {   // scalar divide
    return (1.0 / c) * a;
}

// Geometric data at a gridpoint needed by the charged-matter equations (z-derivatives of the BSSN variables and the
// conformal Christoffel symbols used by ComplexScalarField)
struct ChargedGeom {
    double d_z_alpha, d_z_chi, d_z_h_zz, d_z_h_ww, d_z_beta;
    double chris_Zzz, chris_Zww;
    double z;
};

// Matter model: complex scalar field minimally coupled to a U(1) gauge field, plus the Maxwell field itself.
// Analogue of ComplexScalarField; with q = 0 and E = A = Phi_e = 0 its scalar equations reduce to ComplexScalarField's.
class ChargedScalarField {
public:
    explicit ChargedScalarField(Spacetime* st_) : st(st_) {}

    ChargedScalarField() = delete;
    ChargedScalarField(const ChargedScalarField&) = default;
    ChargedScalarField(ChargedScalarField&&) = default;
    ChargedScalarField& operator=(const ChargedScalarField&) = default;
    ChargedScalarField& operator=(ChargedScalarField&&) = default;
    ~ChargedScalarField() = default;

    double q = 0.; // gauge coupling

    // Right-hand sides of the scalar (csf_rhs) and Maxwell (em_rhs) variables at a gridpoint.
    // d_z_csf needs all four first derivatives, d_zz_csf only the phi second derivatives;
    // d_z_em needs all three first derivatives. div_A = D_i A^i is passed in: it must be evaluated in the conservative
    // form (1/W) d_z (W A^z), W = z^n sqrt(gamma) (Spacetime::maxwell_divergence); the pointwise form
    // d_z A^z + n A^z / z + ... makes the first-order Lorenz-gauge system (A_z, Phi_e) unstable at the origin.
    void rhs(const CSF& csf, const EMF& em, const BSSNState& bssn,
             const CSF& d_z_csf, const CSF& d_zz_csf, const EMF& d_z_em, double div_A,
             const ChargedGeom& g, CSF& csf_rhs, EMF& em_rhs) const;

    // Gauge-covariant radial derivative D_z phi = d_z phi + i q A_z phi (returned as {Re, Im} in phi_re / phi_im)
    CSF cov_d_z(const CSF& csf, const EMF& em, const CSF& d_z_csf) const;

    // Stress-energy projections (scalar + electromagnetic), as in ComplexScalarField. Only the phi first derivatives
    // of d_z_csf are used.
    double rho (const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const;
    double j_z (const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const;
    double S_zz(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const;
    double S_ww(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const;

    // scalar-field and electromagnetic parts of rho separately (rho = rho_scalar + rho_em)
    double rho_scalar(const CSF& csf, const EMF& em, const BSSNState& bssn, const CSF& d_z_csf) const;
    double rho_em(const EMF& em, const BSSNState& bssn) const;

    // charge density rho_e = -n_a j^a and the covariant radial current J_z = gamma_z^a j_a
    double charge_density(const CSF& csf) const;
    double current_z(const CSF& csf, const EMF& em, const CSF& d_z_csf) const;

    // Gauss constraint D_i E^i - 4 pi rho_e, given div_E = D_i E^i (Spacetime::maxwell_divergence)
    double gauss(const CSF& csf, double div_E) const;

    Spacetime* spacetime() const { return st; }

private:
    Spacetime* st; // non-owning pointer to the spacetime context
};

#endif // CHARGED_SCALAR_FIELD_H
