#ifndef REAL_SCALAR_FIELD_H
#define REAL_SCALAR_FIELD_H

#include <cstddef>

// Forward declarations to avoid heavy includes in this header
struct BSSNState;   // defined in EvolutionVariables.h
class Spacetime;    // defined in EvolutionVariables.h/.cpp

// Container for real scalar field variables on a slice
struct RSF {
    double psi;     // field value (real)
    double K_psi;   // conjugate momentum (real)
};

// Elementwise arithmetic for RSF (mirroring CSF operators)
RSF operator+(const RSF& a, const RSF& b);
RSF operator-(const RSF& a, const RSF& b);
RSF operator*(const RSF& a, const RSF& b);       // termwise multiply
RSF operator*(double c, const RSF& a);           // scalar multiply (left)
inline RSF operator/(const RSF& a, double c) {   // scalar divide
    return (1.0 / c) * a;
}

// Matter model: Massless Real Scalar Field coupled to a given Spacetime
class RealScalarField {
public:
    // Construct with a pointer to the owning/associated Spacetime
    explicit RealScalarField(Spacetime* st_) : st(st_) {}

    // Default rule of five (no special ownership)
    RealScalarField() = delete;
    RealScalarField(const RealScalarField&) = default;
    RealScalarField(RealScalarField&&) = default;
    RealScalarField& operator=(const RealScalarField&) = default;
    RealScalarField& operator=(RealScalarField&&) = default;
    ~RealScalarField() = default;

    // Right-hand-sides for the RSF evolution variables at a gridpoint
    RSF rhs_RSF(const RSF& rsf, const BSSNState& bssn,
                const RSF& d_z_rsf, const RSF& d_zz_rsf,
                double d_z_alpha, double d_z_chi,
                double chris_Zzz, double chris_Zww,
                double z) const;

    // Stress-energy projections T_{ab} against 3+1 split
    // These return the local contributions at a gridpoint, given rsf and metric data in bssn
    double rho(const RSF& rsf, const BSSNState& bssn,
               const RSF& d_z_rsf, const RSF& d_zz_rsf) const;   // energy density
    double j_z(const RSF& rsf, const BSSNState& bssn,
               const RSF& d_z_rsf, const RSF& d_zz_rsf) const;   // momentum density along z (radial)
    double S_zz(const RSF& rsf, const BSSNState& bssn,
                const RSF& d_z_rsf, const RSF& d_zz_rsf) const;  // stress in zz
    double S_ww(const RSF& rsf, const BSSNState& bssn,
                const RSF& d_z_rsf, const RSF& d_zz_rsf) const;  // stress in angular directions (trace over sphere)

    // Access to associated spacetime
    Spacetime* spacetime() const { return st; }

private:
    Spacetime* st; // non-owning pointer to the spacetime context
};

#endif // REAL_SCALAR_FIELD_H
