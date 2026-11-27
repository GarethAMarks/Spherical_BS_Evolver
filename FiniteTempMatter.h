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
    double K_sigma; // -Y Pi_sigma / 2 (even); = K_sigma for delta Y = 0
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
    double B_floor = 0.;   // >0: convex floor on the phonon B factor at B_floor*c^2 (evolution only)
    double q = 1.;          // pressure_interp_power
    int gas_statistics = 0; // 0 Boltzmann, 1 Bose-Einstein free gas (bose_gas.h)
    int crossover_mode = 0; // 0: q-norm (draft eq. 31); 1: convex sigma blend; 2: Bogoliubov (see FiniteTempRelaxer.h)
    double bogoliubov_gap_width = 1.; // crossover_mode 2: softening of the phase-mode gap at mu_c = m, in units of T
    double bogoliubov_sigma_gap = 0.005; // crossover_mode 2: gap switched off where sigma >> bogoliubov_sigma_gap
    double sigma_phonon = 0.023, sigma_gas = 0.017, crossover_exponent = 20.;
    // crossover_mode 3 (EoS R): w(sigma) P_qp + (1 - w) P_BE, w = 1 - exp(-(sigma/sigma_decouple)^decouple_exponent);
    // P_qp = (mu_c T / 2 pi)^{3/2} T f(delta, gam, eta), delta = X / (2 m sigma^2 T) (gradient gap, softened for X <= 0 by
    // gap_soft), gam = sigma^2 V'' / (2 mu_c T), eta = (mu - mu_c) / T; f from bogoliubov_table.h.
    double sigma_decouple = 0.003, decouple_exponent = 8., gap_soft = 0.;

    struct Thermo
    {
        double P_th = 0.;    // thermal pressure
        double s = 0.;       // dP/dT
        double n = 0.;       // dP/dmu
        double dF = 0.;      // delta F = 2 dP_th/d(sigma^2)
        double dUps = 0.;    // delta Upsilon = -2 dP_th/dXi
        double dY = 0.;      // delta Y = -2 dP_th/dX (crossover_mode 3 only)
    };
    // thermal part only (the scalar part is the Lagrangian, handled by SBSE's ComplexScalarField)
    Thermo thermal(double T, double mu, double sigma_sq, double Xi, double X = 0.) const;
    // smoothed gradient gap delta_eff(delta) >= 0 and its derivative
    double soft_gap(double delta, double& d_ddelta) const;
    void phonon_B(double mu, double sigma_sq, double Xi, double& c2, double& B) const;

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
    double n_floor = 1e-40;   // below this N_N (and Sigma) the normal component is treated as absent (T = 0, v = 0)
    double fluid_atm_taper = 1.; // >1: normal-fluid RHS tapered smoothly from 0 at the floor to 1 at taper*floor
    std::vector<double> atm_w;
    bool fluid_atm_fixed = false; // floor relative to the initial max N_N instead of the current one
    double max_NN0 = 0.;
    double av_quad = 0., av_lin = 0., av_cs = 0.1; // artificial viscosity coefficients (0 = off)
    std::vector<double> av_Q, av_theta;
    bool hrsc = false;        // finite-volume LLF + MC transport for N_N, Sigma, P_z
    double hrsc_cmax = 1.;    // LLF wave-speed bound added to |V| (times alpha)
    std::vector<double> vel_V; // alpha v - beta
    double fluid_atm_frac = 1e-6; // fluid atmosphere: N_N < fluid_atm_frac * max N_N -> normal fluid frozen, v = 0
    // heat bath: where sigma < fluid_sigma_min (outer gas region of the convex crossover) the normal component is not evolved;
    // it is held at the Tolman/Klein values T = (T alpha)_initial / alpha, mu = (mu alpha)_initial / alpha, v = 0, and
    // N_N, Sigma are slaved to it
    double fluid_sigma_min = 0.;
    std::vector<double> T_alpha0, mu_alpha0;
    std::vector<char> bath;
    std::vector<char> fluid_atm;
    // conversion between superfluid and normal charge (draft theta_1 term; theta_2 = 0): sources
    // +-alpha theta_1 J in N_S / N_N, alpha theta_1 J^2 / T in Sigma, -alpha theta_1 J d_z psi in P_z, J = W (Pi + v d_z psi) + mu.
    // The relaxation rate of J is ~ alpha theta_1 (1 / (W dn/dmu) + 1 / Upsilon); theta_1 is capped locally at
    // theta_cfl / (dt * rate) so that the explicit RK4 step stays stable (dt = time step of the evolution).
    double theta1 = 0.;
    double theta_cfl = 0.5;
    bool theta_implicit = false;
    double theta_pole_p = 0.; // >0: implicit conversion rate theta_1 (c^2/B)^p // theta_1 by implicit operator splitting (relax_conversion) instead of RHS sources
    double dt = 0.;
    int theta_capped = 0;     // points where theta_1 was capped in the last prepare()
    // atmosphere: where sigma < sigma_atm the polar split phi = sigma e^{i psi} is ill-conditioned (Pi = -N_S / sigma^2), so
    // there the CSF evolves (phi, K_phi) itself, K_sigma and N_S are slaved to it, and the normal fluid is frozen and ignored
    double sigma_atm = 1e-5;
    std::vector<char> atm;
    int c2p_max_iter = 40;
    double c2p_tol = 1e-12;

    // Per-slice primitive and auxiliary caches (filled by prepare(); valid for the slice last prepared)
    struct Point
    {
        double T = 0., mu = 0., v = 0., Pi = 0.;   // primitives
        double sigma = 0., dz_sigma = 0., dz_psi = 0.;
        double W = 1., w = 0., P_th = 0., s = 0., n = 0., dF = 0., dUps = 0., Ups = 0.;
        double Xi = 0.;
        double Ksig = 0.;  // primitive K_sigma = (evolved Y K_sigma) / Y
        double X = 0.;     // X = -4 K_sigma^2 + gamma^zz (d_z sigma)^2
        double dY = 0.;    // delta Y = Y - 1
        double cross = 0.;
        double J = 0., theta = 0.; // deviation from Josephson and the (capped) theta_1 used here   // sigma^2 d_z psi = phi_re d_z phi_im - phi_im d_z phi_re
    };
    std::vector<Point> pts;
    std::vector<double> f_sigma, f_T, f_mu, f_Y, coef_S, flux_S, flux_N, flux_Sig, flux_P; // arrays that get differentiated
    int c2p_failures = 0;     // in the last prepare()
    double max_c2p_res = 0.;  // worst relative residual in the last prepare()

    // Primitive recovery on every point of the slice, reconstruction of csf.K_phi_re/K_phi_im, and the flux arrays.
    void prepare(BSSNSlice* slice);
    // implicit theta_1 conversion over a step h (operator split, after each full RK4 step)
    void relax_conversion(BSSNSlice* slice, double h);

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

    // K_sigma from the evolved Y K_sigma (Kt) at fixed (T, mu, sigma^2, Xi): fixed point K = Kt / Y(X(K)), with
    // X = -4 K^2 + gzz_sz2; returns K and the thermodynamics at the final X
    double solve_Ksig(double T, double mu, double s2, double Xi, double Kt, double gzz_sz2, FiniteTempEoS::Thermo& th, double& X) const;

private:
    Spacetime* st;
    std::vector<Point> guess; // previous primitives, used as the Newton starting guess
    bool recover(const FTM& U, double floor, const BSSNState& b, double sigma, double dz_sigma, double dz_psi, double K_sigma, Point& p, double& res) const;
};

#endif
