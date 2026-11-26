#ifndef FINITE_TEMP_MATTER_H
#define FINITE_TEMP_MATTER_H

// Finite-temperature superfluid matter (the draft's two-component model) for SBSE, following
// finiteT_BSSN_equations.tex with theta_1 = theta_2 = 0:
//   condensate phi = sigma e^{i psi} (phi_re, phi_im still live in CSF), plus the evolved densities
//     K_sigma = -Pi_sigma / 2,  N_S = -Upsilon Pi_psi,  N_N = W n,  Sigma = W s,  P_z = w W^2 v_z,
//   primitives {T, mu, v, Pi_psi} recovered algebraically at every RHS evaluation. SBSE's K_phi_re/K_phi_im are then
//   derived quantities, reconstructed from K_sigma and Pi_psi so that ComplexScalarField's phi equation and stress-energy
//   can be reused unchanged.
// Draft sign convention throughout: phi ~ sigma e^{-i mu_inf t}, mu > 0 (see FiniteTempBackground.h).
// Serial only (no MPI ghost exchange of the new variables).

#include <vector>
#include "FiniteTempBackground.h"

struct BSSNState;
class Spacetime;
class BSSNSlice;
struct State;

// Evolved finite-T variables at a gridpoint (in addition to phi_re, phi_im in CSF)
struct FTM {
    double K_sigma; // -Pi_sigma / 2 (even)
    double N_S;     // superfluid charge density -Upsilon Pi_psi (even)
    double N_N;     // normal charge density W n (even)
    double Sigma;   // entropy density W s (even)
    double P_z;     // normal momentum density w W^2 v_z (odd)
};

FTM operator+(const FTM& a, const FTM& b);
FTM operator-(const FTM& a, const FTM& b);
FTM operator*(const FTM& a, const FTM& b);
FTM operator*(double c, const FTM& a);
inline FTM operator/(const FTM& a, double c) { return (1.0 / c) * a; }

// Equation of state of the relaxer (finite_temp/FiniteTempRelaxer.cpp), with analytic derivatives:
//   P = -(X + sigma^2 Xi + V(sigma^2)) / 2 + P_th,  P_th = (P_ph^-q + P_gas^-q)^(-1/q),
//   P_ph = pi^2 c T^4 / (90 B^2), c^2 = (a + eps) / (a - 2 Xi), a = sigma^2 V'', B = 1 + (1 - c^2) mu^2 / Xi,
//   P_gas = (m / 2 pi)^{3/2} T^{5/2} e^{(mu - m) / T}.
// V is SBSE's potential (V = mu^2 A^2 + lambda A^4, optionally solitonic), written as a function of sigma^2.
struct FiniteTempEoS
{
    double m = 1.;          // scalar mass (SBSE mu)
    double lambda = 0.;
    bool solitonic = false;
    double sigma_scale = 1.;
    double cph_epsilon = 1e-8;
    double q = 1.;          // pressure_interp_power

    struct Thermo
    {
        double P_th = 0.;    // thermal pressure
        double s = 0.;       // dP/dT
        double n = 0.;       // dP/dmu
        double dF = 0.;      // delta F = 2 dP_th/d(sigma^2)
        double dUps = 0.;    // delta Upsilon = -2 dP_th/dXi
    };
    // thermal part only (the scalar part is the Lagrangian, handled by SBSE's ComplexScalarField)
    Thermo thermal(double T, double mu, double sigma_sq, double Xi) const;

    double V(double s2) const;   // V(sigma^2)
    double dV(double s2) const;  // dV/d(sigma^2)
    double ddV(double s2) const; // d^2V/d(sigma^2)^2
    double dddV(double s2) const;
};

class FiniteTempMatter
{
public:
    explicit FiniteTempMatter(Spacetime* st_) : st(st_) {}
    FiniteTempMatter() = delete;

    FiniteTempEoS eos;
    double n_floor = 1e-40;
    // atmosphere: where sigma < sigma_atm the polar split phi = sigma e^{i psi} is ill-conditioned (Pi = -N_S / sigma^2), so
    // there the CSF evolves (phi, K_phi) itself, K_sigma and N_S are slaved to it, and the normal fluid is frozen and ignored
    double sigma_atm = 1e-5;
    std::vector<char> atm;   // below this N_N (and Sigma) the normal component is treated as absent (T = 0, v = 0)
    int c2p_max_iter = 40;
    double c2p_tol = 1e-12;

    // Per-slice primitive and auxiliary caches (filled by prepare(); valid for the slice last prepared)
    struct Point
    {
        double T = 0., mu = 0., v = 0., Pi = 0.;   // primitives
        double sigma = 0., dz_sigma = 0., dz_psi = 0.;
        double W = 1., w = 0., P_th = 0., s = 0., n = 0., dF = 0., dUps = 0., Ups = 0.;
        double Xi = 0.;
        double cross = 0.;   // sigma^2 d_z psi = phi_re d_z phi_im - phi_im d_z phi_re
    };
    std::vector<Point> pts;
    std::vector<double> f_sigma, f_T, f_mu, coef_S, flux_S, flux_N, flux_Sig, flux_P; // arrays that get differentiated
    int c2p_failures = 0;     // in the last prepare()
    double max_c2p_res = 0.;  // worst relative residual in the last prepare()

    // Primitive recovery on every point of the slice, reconstruction of csf.K_phi_re/K_phi_im, and the flux arrays.
    void prepare(BSSNSlice* slice);

    // Initial FTM values for a static background (K_sigma = 0, v = 0, Pi_psi = -mu, Xi = -mu^2 at t = 0).
    void initialize_static(BSSNSlice& slice, const std::vector<FiniteTempPoint>& bg);

    // RHS of the FTM variables at point j (prepare() must have been called on slice)
    FTM rhs(BSSNSlice* slice, int j, double z, double d_z_alpha, double d_z_chi, double d_z_h_zz, double d_z_h_ww,
            double d_z_beta, double chris_Zzz, double chris_Zww) const;

    // FTM RHS in the atmosphere: K_sigma, N_S rates implied by the CSF rates (chain rule), fluid frozen
    FTM atmosphere_rhs(const State& s, const State& rhs) const;

    // Thermal stress-energy contributions at point j (added to the CSF ones)
    void stress_energy(const BSSNState& bssn, int j, double& rho, double& j_z, double& S_zz, double& S_ww) const;

    // derivative of a cached array with the slice's stencil machinery (origin mirror with given parity, refinement)
    double deriv(const BSSNSlice* slice, const std::vector<double>& f, int j, int order, bool odd) const;

private:
    Spacetime* st;
    std::vector<Point> guess; // previous primitives, used as the Newton starting guess
    bool recover(const FTM& U, const BSSNState& b, double sigma, double dz_sigma, double dz_psi, double K_sigma, Point& p, double& res) const;
};

#endif
