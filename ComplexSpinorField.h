#ifndef COMPLEX_SPINOR_FIELD_H
#define COMPLEX_SPINOR_FIELD_H

// Forward declarations to avoid heavy includes in this header
struct BSSNState;   // defined in EvolutionVariables.h
class Spacetime;    // defined in EvolutionVariables.h/.cpp

/*
 Two Dirac spinors of mass mu (the matter of DiracStar), (gamma^a nabla_a + mu) psi = 0 with {gamma^a, gamma^b} = 2 eta^{ab}
 in the (-,+,+,+) signature, occupying the j = 1/2, m_j = +-1/2 s-wave states so that the stress tensor is spherically
 symmetric (D = 4 only):
     psi_m = ( u Omega_{-1,m},  i v Omega_{+1,m} ),
 with u, v complex orthonormal-frame amplitudes (for a stationary star u = e^{-i omega t} g(r), v = e^{-i omega t} f(r), with
 f, g the physical spinor functions written by DiracStar). On the cartoon z-axis u is even and v is odd in z.

 The tetrad is fixed by the time gauge e_0 = n and the symmetric spatial triad (aligned with the diagonal h_zz, h_ww), so it
 is an algebraic function of the BSSN variables. The evolution equations are (' = d_z)
     d_t u = beta u' + (alpha K/2 - i alpha mu_eff) u + i alpha [ sqrt(chi/h_zz) (v' + L v) + (sqrt(chi/h_zz) + sqrt(chi/h_ww)) v / z ]
     d_t v = beta v' + (alpha K/2 + i alpha mu_eff) v - i alpha [ sqrt(chi/h_zz) (u' + L u) + (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z ]
     L = (h_ww'/h_ww - chi'/chi + alpha'/alpha) / 2,
 with mu_eff = mu + lambda S for the quartic self-interaction U(S) = (lambda/2) S^2 (lambda = 0: free field), and, with
 overbars for complex conjugates and hats for orthonormal-frame components,
     n        = (|u|^2 + |v|^2) / (2 pi)                               (number density; Q = int n sqrt(gamma) d^3x)
     S        = (|u|^2 - |v|^2) / (2 pi)                               (scalar density sum psibar psi)
     S_rr^    = Re( vbar u' - u vbar' ) / (2 pi sqrt(g_zz)) + (lambda/2) S^2
     S_thth^  = -Re( u vbar ) / (2 pi z sqrt(g_ww)) + (lambda/2) S^2
     rho      = mu S + (lambda/2) S^2 + [kinetic parts of S_rr^ + 2 S_thth^]
     j_r^     = Im( ubar u' + vbar v' ) / (2 pi sqrt(g_zz)).
 Derived and checked in derivations/check_dirac_3p1.py (free field), derivations/check_dirac_selfint_3p1.py (the
 self-interaction terms) and, for the stationary equations, derivations/check_dirac_star.py.
 At z = 0: v / z -> v', and (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z -> 0 since h_zz - h_ww = O(z^2).
 As for the Proca system, the equations are first order in space; Kreiss-Oliger dissipation (damping_factor ~ 1) is recommended.
*/

struct CSpF {
    double u_re;   // Re(u), upper (s-wave) amplitude, even in z
    double u_im;   // Im(u)
    double v_re;   // Re(v), lower (p-wave) amplitude, odd in z
    double v_im;   // Im(v)
};

CSpF operator+(const CSpF& a, const CSpF& b);
CSpF operator-(const CSpF& a, const CSpF& b);
CSpF operator*(const CSpF& a, const CSpF& b);   // termwise multiply
CSpF operator*(double c, const CSpF& a);
inline CSpF operator/(const CSpF& a, double c) { return (1.0 / c) * a; }

class ComplexSpinorField {
public:
    explicit ComplexSpinorField(Spacetime* st_) : st(st_) {}

    ComplexSpinorField() = delete;
    ComplexSpinorField(const ComplexSpinorField&) = default;
    ComplexSpinorField(ComplexSpinorField&&) = default;
    ComplexSpinorField& operator=(const ComplexSpinorField&) = default;
    ComplexSpinorField& operator=(ComplexSpinorField&&) = default;
    ~ComplexSpinorField() = default;

    // right-hand sides at a gridpoint, given the field, its z-derivatives and the needed geometric derivatives
    CSpF rhs_CSpF(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s,
                  double d_z_alpha, double d_z_chi, double d_z_h_ww, double z) const;

    // 3+1 projections of the stress tensor (summed over both spinors); need z-derivatives of u, v
    double rho(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const;   // n^a n^b T_ab
    double j_z(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s) const;             // -gamma_z^a n^b T_ab (lower z index)
    double S_zz(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s) const;            // gamma_z^a gamma_z^b T_ab
    double S_ww(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const;  // angular component (lower w index)

    // number density n = sum psi^dag psi (orthonormal frame); the particle number is int n sqrt(gamma) d^3x
    double number_density(const CSpF& s) const;

    Spacetime* spacetime() const { return st; }

private:
    Spacetime* st; // non-owning pointer to the spacetime context

    // orthonormal-frame stresses S_rr^, S_thth^ (shared by rho and the S_ij)
    double S_rr_hat(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s) const;
    double S_tt_hat(const CSpF& s, const BSSNState& bssn, const CSpF& d_z_s, double z) const;
};

#endif // COMPLEX_SPINOR_FIELD_H
