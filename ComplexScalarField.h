#ifndef COMPLEX_SCALAR_FIELD_H
#define COMPLEX_SCALAR_FIELD_H

#include <cstddef>

// Forward declarations to avoid heavy includes in this header
struct BSSNState;   // defined in EvolutionVariables.h
class Spacetime;    // defined in EvolutionVariables.h/.cpp

// Container for complex scalar field variables on a slice
struct CSF {
    double phi_re;    // Re(phi)
    double phi_im;    // Im(phi)
    double K_phi_re;  // Re(conjugate momentum)
    double K_phi_im;  // Im(conjugate momentum)
};

// Elementwise arithmetic for CSF
CSF operator+(const CSF& a, const CSF& b);
CSF operator-(const CSF& a, const CSF& b);
CSF operator*(const CSF& a, const CSF& b);       // termwise multiply
CSF operator*(double c, const CSF& a);           // scalar multiply (left)
inline CSF operator/(const CSF& a, double c) {   // scalar divide
    return (1.0 / c) * a;
}

// Matter model: Complex Scalar Field coupled to a given Spacetime
class ComplexScalarField {
public:
    // Construct with a pointer to the owning/associated Spacetime
    explicit ComplexScalarField(Spacetime* st_) : st(st_) {}

    // Default rule of five (no special ownership)
    ComplexScalarField() = delete;
    ComplexScalarField(const ComplexScalarField&) = default;
    ComplexScalarField(ComplexScalarField&&) = default;
    ComplexScalarField& operator=(const ComplexScalarField&) = default;
    ComplexScalarField& operator=(ComplexScalarField&&) = default;
    ~ComplexScalarField() = default;

    // Right-hand-sides for the CSF evolution variables at a gridpoint
    // Inputs are the current field state (csf) and the corresponding BSSN state (bssn)
    CSF rhs_CSF(const CSF& csf, const BSSNState& bssn,
            const CSF& d_z_csf, const CSF& d_zz_csf,
            double d_z_alpha, double d_z_chi,
            double chris_Zzz, double chris_Zww,
            double z) const;

    // Stress-energy projections T_{ab} against 3+1 split
    // These return the local contributions at a gridpoint, given csf and metric data in bssn
        double rho(const CSF& csf, const BSSNState& bssn,
                   const CSF& d_z_csf, const CSF& d_zz_csf) const;   // energy density
        double j_z(const CSF& csf, const BSSNState& bssn,
                   const CSF& d_z_csf, const CSF& d_zz_csf) const;   // momentum density along z (radial)
        double S_zz(const CSF& csf, const BSSNState& bssn,
                    const CSF& d_z_csf, const CSF& d_zz_csf) const;  // stress in zz
        double S_ww(const CSF& csf, const BSSNState& bssn,
                    const CSF& d_z_csf, const CSF& d_zz_csf) const;  // stress in angular directions (trace over sphere)

    // Access to associated spacetime
    Spacetime* spacetime() const { return st; }

private:
    Spacetime* st; // non-owning pointer to the spacetime context
};

#endif // COMPLEX_SCALAR_FIELD_H
