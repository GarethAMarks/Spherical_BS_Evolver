#ifndef COMPLEX_VECTOR_FIELD_H
#define COMPLEX_VECTOR_FIELD_H

#include <cstddef>

// Forward declarations to avoid heavy includes in this header
struct BSSNState;   // defined in EvolutionVariables.h
class Spacetime;    // defined in EvolutionVariables.h/.cpp

/*
 Complex Proca field, L = -1/4 F_ab Fbar^ab - 1/2 mu^2 A_a Abar^a (same normalisation as ProcaStar), evolved in the
 3+1 form of Zilhao, Witek & Cardoso (2015) with a Gauss-constraint damping variable Z. Each real component obeys
     d_t X_i   = Lie_beta X_i - alpha (E_i + D_i X_phi) - X_phi D_i alpha
     d_t E^i   = Lie_beta E^i + alpha K E^i + alpha mu^2 X^i + alpha D^i Z       (curl term vanishes in spherical symmetry)
     d_t X_phi = Lie_beta X_phi + alpha K X_phi - alpha D_i X^i - X^i D_i alpha - alpha Z
     d_t Z     = Lie_beta Z + alpha (D_i E^i + mu^2 X_phi - kappa Z)
 with X_i = gamma_i^a A_a, X_phi = -n^a A_a, E_i = gamma_i^a F_ab n^b, Gauss constraint C = D_i E^i + mu^2 X_phi = 0.
 In spherical symmetry only X_z (lower index) and E^z (upper index) are nonzero on the w = 0 line; both are odd in z,
 while X_phi and Z are even. The modified-cartoon reductions of the divergences and Lie derivatives used here are
 checked against brute-force Cartesian calculations in derivations/check_cartoon_reduction.py, and the 3+1 system
 itself against the covariant field equations in derivations/check_proca_3p1.py.
 Being first order in space, the system needs Kreiss-Oliger dissipation (damping_factor ~ 1) to be stable at the origin.
*/

// Container for complex Proca field variables on a slice
struct CVF {
    double X_re;     // Re(X_z), radial covariant component of the spatial projection of A
    double X_im;     // Im(X_z)
    double E_re;     // Re(E^z), radial contravariant component of the electric field
    double E_im;     // Im(E^z)
    double Xphi_re;  // Re(X_phi) = -n^a Re(A_a)
    double Xphi_im;  // Im(X_phi)
    double Z_re;     // Re(Z), Gauss-constraint damping variable
    double Z_im;     // Im(Z)
};

// Elementwise arithmetic for CVF
CVF operator+(const CVF& a, const CVF& b);
CVF operator-(const CVF& a, const CVF& b);
CVF operator*(const CVF& a, const CVF& b);       // termwise multiply
CVF operator*(double c, const CVF& a);           // scalar multiply (left)
inline CVF operator/(const CVF& a, double c) {   // scalar divide
    return (1.0 / c) * a;
}

// Matter model: Complex Proca field coupled to a given Spacetime
class ComplexVectorField {
public:
    explicit ComplexVectorField(Spacetime* st_) : st(st_) {}

    ComplexVectorField() = delete;
    ComplexVectorField(const ComplexVectorField&) = default;
    ComplexVectorField(ComplexVectorField&&) = default;
    ComplexVectorField& operator=(const ComplexVectorField&) = default;
    ComplexVectorField& operator=(ComplexVectorField&&) = default;
    ~ComplexVectorField() = default;

    // Right-hand-sides for the CVF evolution variables at a gridpoint, given the field state, its z-derivatives
    // (d_z_cvf) and the needed geometric derivatives.
    CVF rhs_CVF(const CVF& cvf, const BSSNState& bssn, const CVF& d_z_cvf,
                double d_z_alpha, double d_z_beta, double d_z_chi,
                double d_z_h_zz, double d_z_h_ww, double z) const;

    // Stress-energy projections (no derivatives are needed in spherical symmetry since the magnetic part vanishes)
    double rho(const CVF& cvf, const BSSNState& bssn) const;   // energy density n^a n^b T_ab
    double j_z(const CVF& cvf, const BSSNState& bssn) const;   // momentum density -gamma_z^a n^b T_ab
    double S_zz(const CVF& cvf, const BSSNState& bssn) const;  // gamma_z^a gamma_z^b T_ab
    double S_ww(const CVF& cvf, const BSSNState& bssn) const;  // angular component

    // Noether charge density with respect to the coordinate volume element sqrt(gamma) d^{D-1}x:
    // -n_a j^a = Im(E)^z Re(X)_z - Re(E)^z Im(X)_z
    double charge_density(const CVF& cvf) const;

    // Gauss constraint C = D_i E^i + mu^2 X_phi (real and imaginary parts), for diagnostics
    void gauss_constraint(const CVF& cvf, const BSSNState& bssn, const CVF& d_z_cvf,
                          double d_z_chi, double d_z_h_zz, double d_z_h_ww, double z,
                          double& C_re, double& C_im) const;

    Spacetime* spacetime() const { return st; }

private:
    Spacetime* st; // non-owning pointer to the spacetime context

    // D_i V^i for a radial vector with upstairs component V^z, given d_z V^z and (for regularisation at the
    // origin) the value to use for V^z / z there
    double divergence(double V, double dV, double V_over_z_origin, const BSSNState& bssn,
                      double d_z_chi, double d_z_h_zz, double d_z_h_ww, double z) const;
};

#endif // COMPLEX_VECTOR_FIELD_H
