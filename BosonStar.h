#ifndef BOSONSTAR_HPP_
#define BOSONSTAR_HPP_

#include<vector>
#include<math.h>
#include <iostream>
#include<string>
#include <fstream>
#include <quadmath.h>


//holds field amplitude, conformal factor, their derivative, and the lapse
struct FieldState
{
    double A {}; // scalar field modulus
    double X {}; // conformal factor (squared)
    double phi {}; // log(lapse)
    double eta {}; // conformal factor gradient
};

//overloads for adding/subtracting//scalar mult by FieldState struct
FieldState operator+(const FieldState& s1, const FieldState& s2);
FieldState operator-(const FieldState& s1, const FieldState& s2);
FieldState operator*(double c, const FieldState& s);

class BosonStar
{
    private:
        double mu; //scalar mass
        double lambda; //interaction term coefficient
        bool solitonic; //1 for solitonic BS, 0 for mini
        double G;
        double sigma;
        int eigen; //desired eigenvalue, 0 for ground state
        double alpha_central; //central lapse
        bool isotropic; //whether we will use isotropic coordinates (otherwise polar-areal)

        double frequency_guess; //initial guess for eigenfrequency
        double r_match_fac;


        std::vector<FieldState> state;
        std::vector<double> radius_array; // radius

        std::vector<double> phi_iso_array; // log(lapse) in isotropic coords
        std::vector<double> psi_iso_array; // conformal factor in isotropic coords
        std::vector<double> A_iso_array; // field amplitude in isotropic coords

         std::vector<double> r_iso_array; // isotropic radii given in terms of areal indices


        double V (const double p);
        double dV (const double p);

        FieldState state_RHS(const double radius, const long double frequency, FieldState  s, bool asymptotic_region, bool given_A = 0);


        double R; //max radius of computational domain
        int n_gridpoints; // number of spatial gridpoints
        double courant_factor;
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

    public:
        double A_central; //central field modulus
        long double omega_pre_rescale; //frequency before rescaling to enforce lapse BC, suitable
        long double omega; //solved frequency
        double M; //solved mass
        double binding_energy; //M - mu*N
        double noether_charge;
        double compactness;
        double r_99; //solved radius containing 99% of mass

        BosonStar() = default;
        //BosonStar(const BosonStar& boson_star); //copy constructor(don't seem to need for now)
        void rk4_solve (const long double freq);
        void read_parameters(bool quiet = 0);
        void read_thinshell();
        void write_field();
        void convergence_test(long double freq = 0.);
        void double_resolution();
        bool solve(bool quiet = 0);
        bool solve_finding_A(long double freq, double A_guess, double A_range, bool quiet = 0);
        double get_noether_charge();

        void fill_isotropic_arrays();
        void write_isotropic();

        void fill_given_A( const double freq);

        void cycle_models(int n_stars, double A_0, double delta_A);





    //declare BSSNSlice a friend class so we can construct one from a boson_star object
    friend class BSSNSlice;
    friend class Spacetime;


};





#endif /* BOSONSTAR_HPP_ */
