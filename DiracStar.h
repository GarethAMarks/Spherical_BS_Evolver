#ifndef DIRACSTAR_HPP_
#define DIRACSTAR_HPP_

#include<vector>
#include<math.h>
#include "DiracPotential.h"
#include <iostream>
#include<string>
#include <fstream>

/*
 Construction of spherical Dirac stars (Einstein-Dirac solitons), mirroring the BosonStar and ProcaStar classes.

 Conventions (identical to BosonStar/ProcaStar where they overlap), signature (-,+,+,+):
   metric:  ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,   alpha = exp(phi), r areal, D = 4 only
   matter:  two Dirac spinors of mass mu, L = -psibar (gamma^a nabla_a + mu) psi each ({gamma^a, gamma^b} = 2 eta^{ab}),
            occupying the j = 1/2 (kappa = -1, s-wave) states with m_j = +1/2 and -1/2, so that T_ab is spherically
            symmetric (the configuration of Finster, Smoller & Yau 1999 and Herdeiro, Pombo & Radu 2017), with an optional
            self-interaction -U(S) in L, S = sum_m psibar_m psi_m (DiracPotential.h): the quartic (Soler-type)
            U = (lambda/2) S^2 and/or, with solitonic = 1, the degenerate-vacuum potential mu S + U = mu S (1 - S/S0)^2,
            S0 = sigma^2 / (4 pi). Each spinor then obeys (gamma^a nabla_a + mu_eff) psi = 0 with mu_eff = mu + U'(S).
            For the quartic, lambda > 0 is repulsive (with the sign conventions of Dzhunushaliev & Folomeev 2019,
            lambda = -lambda_DF and lambda_bar_DF = -lambda / (4 pi G) for mu = 1).
   ansatz:  psi_m = e^{-i omega t} alpha^{-1/2} r^{-1} ( G(r) Omega_{-1,m}, i F(r) Omega_{+1,m} )   (Dirac representation,
            Omega_{kappa,m} the spinor spherical harmonics, Omega_{+1,m} = -(sigma . n) Omega_{-1,m}).

 With this ansatz the Dirac equation separates into
       G' = X [  G/r + (mu_eff + omega/alpha) F ],      F' = X [ -F/r + (mu_eff - omega/alpha) G ],
 with mu_eff = mu + U'(S) and S = (G^2 - F^2) / (2 pi alpha r^2), and the stress tensor summed over both spinors is
 (orthonormal frame; T_ab = T^kin_ab + g_ab P(S), P = S U' - U)
       n   = (F^2 + G^2) / (2 pi alpha r^2)                      (number density; Q = int n sqrt(gamma) d^3x)
       rho = (omega / alpha) n - P(S)
       p_r = [ (omega/alpha)(F^2 + G^2) - mu (G^2 - F^2) + 2 F G / r ] / (2 pi alpha r^2) - U(S)
       p_t = -F G / (2 pi alpha r^3) + P(S),
 so that, as for the boson and Proca stars,
       X'   = X (1 - X^2) / (2 r) + 4 pi G r X^3 rho,      phi' = (X^2 - 1) / (2 r) + 4 pi G r X^2 p_r.
 These are derived (and the separation, symmetry and conservation of T_ab checked) in derivations/check_dirac_star.py.
 At strong repulsive coupling (lambda mu^2 >> G^-1) the star approaches a fluid with rho = mu n + lambda n^2 / 2,
 p = lambda n^2 / 2 (mu_eff = omega/alpha locally, F -> 0), whose TOV sequence has M_max = 0.1158 sqrt(lambda) / G^(3/2)
 (Dzhunushaliev & Folomeev: 0.41 sqrt(|lambda_bar|) M_Pl^2 / mu) and M/R = 0.21 at the maximum mass.

 Regularity at the origin: G = g0 sqrt(alpha0) r + O(r^3), F = (mu_eff0 - omega/alpha0) g0 sqrt(alpha0) r^2 / 3 + O(r^4), with
 mu_eff0 = mu_eff(g0^2 / (2 pi)), where
 g0 = |upper spinor component| * sqrt(4 pi) at r = 0 is the physical central amplitude (invariant under the rescaling
 alpha -> k alpha, omega -> k omega, (F, G) -> sqrt(k) (F, G) used to impose alpha(infinity) = 1). A_central = g0.

 Shooting on omega counts the nodes of G: the ground state has nodeless F and G, eigen = n gives n nodes of G.

 The stored state variables are F, G, X and phi. D = 4 only.
 For small sigma the shooting fails (thin walls); relaxation/dirac_relax.py (collocation with the frequency as a grid
 variable, see README) constructs those models.
*/

struct DiracFieldState
{
    double F {}; // lower radial spinor function (includes the alpha^{1/2} r factor), vanishes as r^2 at the origin
    double G {}; // upper radial spinor function (includes the alpha^{1/2} r factor), vanishes as r at the origin
    double X {}; // sqrt(g_rr)
    double phi {}; // log(lapse)
    DiracFieldState() = default;
    DiracFieldState(double F_, double G_, double X_, double phi_) : F(F_), G(G_), X(X_), phi(phi_) {}
};

DiracFieldState operator+(const DiracFieldState& s1, const DiracFieldState& s2);
DiracFieldState operator-(const DiracFieldState& s1, const DiracFieldState& s2);
DiracFieldState operator*(double c, const DiracFieldState& s);

class DiracStar
{
    private:
        double mu; //fermion mass
        double lambda; //quartic self-coupling, U(S) = (lambda/2) S^2 with S = sum psibar psi (lambda > 0 repulsive)
        bool solitonic; //degenerate-vacuum potential mu S (1 - S/S0)^2, S0 = sigma^2 / (4 pi)
        DiracPotential pot; //U(S), U'(S), set from mu, lambda, solitonic, sigma by read_parameters
        double G; //Newton's constant
        double sigma;
        int eigen; //desired number of nodes in G, 0 for ground state
        double alpha_central; //central lapse (initial guess, rescaled after shooting)

        double frequency_guess; //initial guess for eigenfrequency; must exceed the (unrescaled) true value
        double r_match_fac;

        std::vector<DiracFieldState> state;
        std::vector<double> radius_array; // radius

        std::vector<double> phi_iso_array; // log(lapse) in isotropic coords
        std::vector<double> psi_iso_array; // conformal factor in isotropic coords
        std::vector<double> f_iso_array; // physical lower spinor function F / (sqrt(alpha) r) at isotropic radii (a scalar under r -> R)
        std::vector<double> g_iso_array; // physical upper spinor function G / (sqrt(alpha) r) at isotropic radii
        std::vector<double> r_iso_array; // isotropic radii given in terms of areal indices

        DiracFieldState state_RHS(const double radius, long double frequency, DiracFieldState s, bool asymptotic_region);
        void matter_terms(double r, long double frequency, const DiracFieldState& s, double& rho, double& p_r, double& n);
        double scalar_density(double r, const DiracFieldState& s); //S = sum psibar psi = (G^2 - F^2) / (2 pi alpha r^2)

        double D;
        double R; //max radius of computational domain
        int n_gridpoints; // number of spatial gridpoints
        double courant_factor; //evolution parameters, read for passing to Spacetime
        double stop_time;
        int blowup_point; //gridpoint at which the solution first blows up
        double freq_epsilon;//tolerance for the frequency finder

        int count_zero_crossings();
        long double find_frequency(bool quiet = 0);
        int find_last_minimum();
        bool fill_asymptotic(bool quiet = 0);
        double m(int j);
        void rescale_lapse (double phi_shift);
        double f_RHS(double r, double f);
        double r_areal (int j_iso);
        void check_supported();

    public:
        double A_central; //central amplitude g0 (see above)
        long double omega_pre_rescale; //frequency before rescaling to enforce lapse BC
        long double omega; //solved frequency
        double M; //solved mass
        double binding_energy; //M - mu*N
        double noether_charge; //total particle number of both spinors
        double compactness;
        double r_99; //solved radius containing 99% of mass

        //Gaussian perturbation of the spinor amplitudes, (u, v) -> (1 + perturb_amp exp(-(r - perturb_center)^2 / perturb_spread^2)) (u, v),
        //applied before evolution with the metric re-solved (same parameters as for the boson and Proca stars)
        bool perturb = 0;
        double perturb_amp = 0.;
        double perturb_spread = 1.;
        double perturb_center = 0.;

        double A0;
        double dA;
        int n_stars;

        bool isotropic; //whether we will use isotropic coordinates (otherwise polar-areal)
        bool cycle_only = 0;

        DiracStar() = default;
        void rk4_solve (const long double freq);
        void read_parameters(bool quiet = 0, std::string filename = "BSParams.par");
        void write_field(std::string filename="DSdata.dat");
        void convergence_test(long double freq = 0.);
        bool solve(bool quiet = 0);
        double get_noether_charge();
        double get_c4(); //compactness from Collodel et. al's C4 definition, as in BosonStar
        double get_cmax(); //get maximum of aspect mass / r
        double g_central() const {return A_central;}

        void fill_isotropic_arrays();
        void write_isotropic();

        void cycle_models(int n_stars, double A_0, double delta_A);

        void add_perturbation(double a, double k, double center, bool quiet = 0);

        void broadcast();

    friend class BSSNSlice;  //declare friend classes so we can construct them from a dirac_star object accessing private members
    friend class Spacetime;
};

#endif /* DIRACSTAR_HPP_ */
