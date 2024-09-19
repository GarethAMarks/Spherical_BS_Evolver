#ifndef BOSONSTAR_HPP_
#define BOSONSTAR_HPP_

#include<vector>
#include<math.h>
#include <iostream>
#include<string>
#include <fstream>


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
        bool solitonic; //1 for solitonic BS, 0 for mini
        double G;
        double sigma;
        int eigen; //desired eigenvalue, 0 for ground state
        double alpha_central; //central lapse
        bool isotropic; //whether we will use isotropic coordinates (otherwise polar-areal)

        double frequency_guess; //initial guess for eigenfrequency


        std::vector<FieldState> state;
        std::vector<double> radius_array; // radius

        std::vector<double> phi_iso_array; // log(lapse) in isotropic coords
        std::vector<double> psi_iso_array; // conformal factor in isotropic coords
        std::vector<double> A_iso_array; // field amplitude in isotropic coords

         std::vector<double> r_iso_array; // isotropic radii given in terms of areal indices


        double V (const double p);
        double dV (const double p);

        FieldState state_RHS(const double radius, const double frequency, FieldState  s, bool asymptotic_region);

        double R; //max radius of computational domain
        int n_gridpoints; // number of spatial gridpoints
        double courant_factor;
        double stop_time;
        int blowup_point; //gridpoint at which the solution first blows up
        double freq_epsilon;//tolerance for the frequency finder

        int count_zero_crossings();
        double find_frequency(bool quiet = 0);
        int find_last_minimum();
        bool fill_asymptotic(bool quiet = 0);
        double m(int j);
        void rescale_lapse (double phi_shift);
        double f_RHS(double r, double f);
        double r_areal (int j_iso);

    public:
        double A_central; //central field modulus
        double omega_pre_rescale; //frequency before rescaling to enforce lapse BC, suitable
        double omega; //solved frequency
        double M; //solved mass
        double binding_energy; //M - mu*N
        double noether_charge;
        double r_99; //solved radius containing 99% of mass

        BosonStar() = default;
        //BosonStar(const BosonStar& boson_star); //copy constructor(don't seem to need for now)
        void rk4_solve (const double freq);
        void read_parameters(bool quiet = 0);
        void write_field();
        void convergence_test(double freq = 0.);
        void double_resolution();
        bool solve(bool quiet = 0);

        void cycle_models(int n_stars, double A_0, double delta_A);

        double get_noether_charge();

        void fill_isotropic_arrays();
        void write_isotropic();


    //declare BSSNSlice a friend class so we can construct one from a boson_star object
    friend class BSSNSlice;
    friend class Spacetime;


};





#endif /* BOSONSTAR_HPP_ */
