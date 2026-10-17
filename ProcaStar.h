#ifndef PROCASTAR_HPP_
#define PROCASTAR_HPP_

#include<vector>
#include<math.h>
#include <iostream>
#include<string>
#include <fstream>

/*
 Construction of spherical Proca stars (PSs), mirroring the BosonStar class.

 Conventions (identical to BosonStar where they overlap):
   metric:  ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega_{D-2}^2,   alpha = exp(phi), r areal
   action:  S = int sqrt(-g) [ R / (16 pi G) - 1/4 F_ab Fbar^ab - 1/2 V(Psi) ],   Psi = A_a Abar^a
   field:   A = e^{-i omega t} ( f(r) dt + i g(r) dr )   =>   Psi = -f^2 / alpha^2 + g^2 / X^2

 V is the same family of potentials used for the BS, with |phi|^2 replaced by Psi:
   mini/quartic: V = mu^2 Psi + lambda Psi^2;   solitonic: V = mu^2 Psi (1 - 2 Psi / sigma^2)^2 + lambda Psi^2
 and dV, ddV are derivatives with respect to Psi (as dV, ddV are with respect to |phi|^2 in BosonStar).

 Evolved ODE variables are f, X, phi and E := f' - omega g (so that F_tr = -E e^{-i omega t}).
 The radial component of the field equation nabla_a F^ab = V'(Psi) A^b is algebraic,
       V'(Psi) g = -omega E / alpha^2,
 so g is not an ODE variable; it is reconstructed from (f, X, phi, E) by solve_g(). For the free (mini)
 potential this is explicit; with self-interactions it is solved by Newton iteration and can fail to have
 a solution (signalled by a nan), which is treated exactly like a blowup of the shooting solution.

 Shooting on omega counts nodes of E (equivalently g), not f: since g ~ -f(0) r near the origin while g and f share a
 sign asymptotically, the ground state has nodeless g and E but f with one node. E then plays the role of the BS amplitude A.

 A_central is defined as f(0) / alpha(0) = sqrt(-Psi(0)), which is invariant under the lapse/frequency
 rescaling performed after shooting, playing the same role as the central scalar amplitude for a BS.
 Under that rescaling (alpha -> k alpha, omega -> k omega) f and E also scale by k while g is unchanged.
*/

//holds A_t, conformal factor, log(lapse), and E = f' - omega g
struct ProcaFieldState
{
    double f {}; // time component of the vector potential (A_t, up to the e^{-i omega t} phase)
    double X {}; // conformal factor (squared)
    double phi {}; // log(lapse)
    double E {}; // f' - omega g; (minus) the electric field
    ProcaFieldState() = default;
    ProcaFieldState(double f_, double X_, double phi_, double E_) : f(f_), X(X_), phi(phi_), E(E_) {}
};

//overloads for adding/subtracting//scalar mult by ProcaFieldState struct
ProcaFieldState operator+(const ProcaFieldState& s1, const ProcaFieldState& s2);
ProcaFieldState operator-(const ProcaFieldState& s1, const ProcaFieldState& s2);
ProcaFieldState operator*(double c, const ProcaFieldState& s);

class ProcaStar
{
    private:
        double mu; //vector field mass
        double lambda; //quartic interaction term coefficient
        bool solitonic; //1 for solitonic-type potential, 0 for mini/quartic
        double G;
        double sigma;
        int eigen; //desired number of nodes in E (equivalently g), 0 for ground state (whose f has one node)
        double alpha_central; //central lapse (initial guess, rescaled after shooting)

        double frequency_guess; //initial guess for eigenfrequency; must exceed the (unrescaled) true value
        double r_match_fac;

        std::vector<ProcaFieldState> state;
        std::vector<double> radius_array; // radius
        std::vector<double> g_array; // radial component g of the vector potential, reconstructed after solving

        std::vector<double> phi_iso_array; // log(lapse) in isotropic coords
        std::vector<double> psi_iso_array; // conformal factor in isotropic coords
        std::vector<double> f_iso_array; // A_t in isotropic coords
        std::vector<double> g_iso_array; // A_R (radial component w.r.t. isotropic radius) in isotropic coords
        std::vector<double> r_iso_array; // isotropic radii given in terms of areal indices

        double V (const double Psi);
        double dV (const double Psi);
        double ddV (const double Psi);
        bool self_interacting() const {return (solitonic || lambda != 0.);}

        double solve_g(const long double frequency, const ProcaFieldState& s, double& W, double& Psi);
        ProcaFieldState state_RHS(const double radius, long double frequency, ProcaFieldState s, bool asymptotic_region);
        ProcaFieldState state_expansion(const double radius, long double frequency);

        double D;
        double R; //max radius of computational domain
        int n_gridpoints; // number of spatial gridpoints
        int blowup_point; //gridpoint at which the solution first blows up
        double freq_epsilon;//tolerance for the frequency finder
        int max_newton_iterations = 50; //for the algebraic solve for g with self-interactions
        double newton_tol = 1e-14;
        bool g_solve_failed = 0; //set if solve_g failed during the current rk4_solve
        double g_breakdown_r = -1.; //radius at which solve_g first failed in the current rk4_solve (negative if it did not)

        int count_zero_crossings();
        long double find_frequency(bool quiet = 0);
        int find_last_minimum();
        bool fill_asymptotic(bool quiet = 0);
        void fill_g_array();
        double m(int j);
        void rescale_lapse (double phi_shift);
        double f_RHS(double r, double f);
        double r_areal (int j_iso);

    public:
        double A_central; //central gauge-invariant amplitude f(0)/alpha(0) = sqrt(-Psi(0))
        long double omega_pre_rescale; //frequency before rescaling to enforce lapse BC
        long double omega; //solved frequency
        double M; //solved mass
        double binding_energy; //M - mu*N
        double noether_charge;
        double compactness;
        double r_99; //solved radius containing 99% of mass

        double A0;
        double dA;
        int n_stars;

        bool isotropic; //whether we will use isotropic coordinates (otherwise polar-areal)
        bool cycle_only = 0;

        ProcaStar() = default;
        void rk4_solve (const long double freq);
        void read_parameters(bool quiet = 0, std::string filename = "BSParams.par");
        void write_field(std::string filename="PSdata.dat");
        void convergence_test(long double freq = 0.);
        void double_resolution();
        bool solve(bool quiet = 0);
        double get_noether_charge();
        double get_c4(); //compactness from Collodel et. al's C4 definition, as in BosonStar
        double get_cmax(); //get maximum of aspect mass / r
        double get_g(int j) const {return g_array[j];}
        double f_central() const {return state[0].f;} //f(0) in the gauge alpha(infinity) = 1, after solving

        void fill_isotropic_arrays();
        void write_isotropic();

        void cycle_models(int n_stars, double A_0, double delta_A);

        void broadcast();

    friend class BSSNSlice;  //declare friend classes so we can construct them from a proca_star object accessing private members
    friend class Spacetime;
};

#endif /* PROCASTAR_HPP_ */
