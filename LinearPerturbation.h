#ifndef LINEARPERTURBATION_HPP_
#define LINEARPERTURBATION_HPP_

#include<vector>
#include<math.h>
#include <iostream>
#include <fstream>
#include <quadmath.h>

//F, L are shorthand for delta psi_1 and delta lambda, the (fractional) perturbations in the scalar field real part and areal rr metric component. Lp, Fp are radial derivatives
struct PertState
{
    double F {};
    double Fp {};
    double L {};
    double Lp {};
};

using namespace std;

//class to hold constructed linear perturbations of a background boson star, a pointer to which is stored as a member variable and can be used in constructor
class LinearPerturbation
{
    public:
        BosonStar* bg; //pointer to the background boson star

        LinearPerturbation(BosonStar* bg, double given_chi_sq0, long double given_gamma0, double given_chi_range, double given_cutoff_radius) //constructor takes a pointer to BosonStar and autofills inherited params
        :bg{bg}, chi_sq0{given_chi_sq0}, gamma0{given_gamma0}, n_gridpoints{bg->n_gridpoints}, R{bg->R}, sigma{bg->sigma}, chi_range{given_chi_range},
        solitonic{bg->solitonic}, omega{bg->omega}, A_central{bg->A_central}, dr{R/(n_gridpoints - 1)}, chi_epsilon{0.00000000001}, cutoff_radius{given_cutoff_radius}
        {}

        void rk4_solve(double chi_sq, long double gamma); //runs RK4 solver on F, L using given chi, gamma
        long double get_best_gamma(double chi_sq, bool quiet = 0);
        long double get_best_gamma_neg(double chi_sq, bool quiet = 0);
        double get_chi_sq();
        double get_chi_sq_newton();
        double get_noether_perturbation(); //returns the perturbation to the Noether charge associated with the computed perturbation
        void write_pert(string filename = "pert.dat");
        void pert_cycle(double A0, double dA, int n_stars);
        void write_chi_results();
        void read_parameters(bool quiet = 0);
        PertState test_rhs (double r, FieldState f, PertState s, double chi_sq, long double gamma);

    private:
        std::vector<PertState> pert; //array to hold the perturbations


        double chi_sq0; // initial guess for characteristic frequency of radial oscillations
        long double gamma0; //initial guess for undetermined constant affecting central initial values

        //inherited BS params
        int n_gridpoints;
        double R;
        double sigma;
        double chi_range;
        bool solitonic;
        long double omega;
        double A_central;
        double dr;

        double chi_epsilon;
        double cutoff_radius; //radius outside of which to look for zero crossings


        double solved_chi_sq;
        long double solved_gamma;

        int blowup_point;

        double noether_perturbation;

        PertState pert_rhs (double r, PertState s, double chi_sq, long double gamma); //returns rhs of the radial ODEs that F, L satisfy
        double last_extremum(); //returns absolute value of last extremum in L
        int count_zero_crossings();



};



#endif /* LINEARPERTURBATION_HPP_ */
