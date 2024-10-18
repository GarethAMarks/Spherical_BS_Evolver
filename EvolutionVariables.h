#ifndef EVOLUTIONVARIABLES_HPP_
#define EVOLUTIONVARIABLES_HPP_
#include <vector>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include "BosonStar.h"

//enum class for the BSSN evolved variables. v_ prefix is to prevent ambiguity with numerical state values
enum bssn_var
{
    v_chi,
    v_h_zz,
    v_h_ww,
    v_A_zz,
    v_A_ww,
    v_K,
    v_c_chris_Z,
    v_phi_re,
    v_phi_im,
    v_K_phi_re,
    v_K_phi_im,
    v_alpha,
    v_beta
};

//a set of BSSN evolution variables at a point in spacetime
struct BSSNState
{
    double chi; //BSSN conformal factor
    double h_zz; //radial component of the conformally rescaled (unit det) metric
    double h_ww; // w denotes auxiliary metric components; everything is finally evaluated on w = 0
    double A_zz; //conformally rescaled traceless extrinsic curvature
    double A_ww;
    double K; //trace of extrinsic curvature
    double c_chris_Z; //z component of contracted conformal christoffel symbols

    double phi_re; //real + imaginary part of scalar field
    double phi_im;

    double K_phi_re; //real + imaginary part of scalar field momentum
    double K_phi_im;

    double alpha; //lapse
    double beta; //shift (only z component nonzero)
};

//the BSSN variables, along with some auxiliary quantities, on a single time slice.
class BSSNSlice
{
    //public vars are BSSN evolution/gauge variables; privates are auxiliary quantities used in evolution/diagnostic equations
    private:

    bool has_BH = 0;
    double make_tangherlini(double m, double min_chi);
    double refinement_level;



    public:
        std::vector<BSSNState> states; //array of states on a BSSN time slice

        std::vector<int> refinement_points; //points at which refinement should be halved; to be copied from spacetime version at start of each step

        void read_BS_data(BosonStar& boson_star,int BS_resolution_factor = 1., bool isotropic = 1);
        void read_checkpoint(int time, int n_gridpoints);
        void write_slice( std::string file_name = "SliceData.dat");

        int get_refinement_level(int j, std::vector<int>& refinement_points);

        bool smooth_lapse();

        //BSSNState end_state1; //states off the boundary, computed by assuming outward spherical wave propagation and used in the evolution of the endpoints via Sommerfeld BCs
        //BSSNState end_state2;

        double R;
        double d_z(bssn_var var, int index, int order); //takes z derivative of a BSSN variable on the current slice centred on given index
        double d_zz(bssn_var var, int index);  //second z-derivative, basically syntactic sugar for d_z(var, index, 2)


    friend class Spacetime; //Spacetime should be able to access private BSSNSlice members
    friend BSSNSlice operator+(const BSSNSlice& slice1, const BSSNSlice& slice2);
    friend BSSNSlice operator*(double c, const BSSNSlice& slice);
};

//holds state information for entire spacetime, and auxiliary quantities on current slice.
class Spacetime
{
    public:
        std::vector<BSSNSlice> slices;//vector of entire time slices

        void initialize(BosonStar& boson_star);
        void write_diagnostics();
        void evolve();

        void halve_resolution();
        void fourier_transform_A0();

//auxiliary variables held on a particular time slice-- CONVENTION: upper/lowercase r,w denote upstairs/downstairs indices where relevent
    private:

        BSSNSlice* current_slice_ptr; //pointer to current slice, updated on each time step

        double D;

        std::vector<double> h_WW;//inverse metric components
        std::vector<double> h_ZZ;


        std::vector<double> chris_Zww;
        std::vector<double> chris_Zzz;//christoffel symbol components we need in evolution eqns

        //maybe unnecessary!
        /*std::vector<double> chris_zww;
        std::vector<double> chris_wwz;
        std::vector<double> chris_Wwz;*/

        std::vector<double> aux_constraint;//auxiliary constraint times divergence of shift

        std::vector<double> D_zz_alpha; //physical 2nd covariant 3-derivatives of lapste
        std::vector<double> D_ww_alpha;
        std::vector<double> D_zz_alpha_TF;//trace-free versions
        std::vector<double> D_ww_alpha_TF;

        std::vector<double> R_zz; //3D Riemann tensor components and trace-free versions
        std::vector<double> R_ww;
        std::vector<double> R_zz_TF;
        std::vector<double> R_ww_TF;

        std::vector<double>  rho; //energy density
        std::vector<double>  j_z; //momentum density
        std::vector<double>  S_zz; //stress tensor
        std::vector<double>  S_ww;

        std::vector<double>  S_zz_TF; //trace and trace-free parts
        std::vector<double>  S_ww_TF;
        std::vector<double>  S;

        //Hamiltionian and Momentum constraints, and determinant of h
        std::vector<double>  Ham;
        std::vector<double>  Mom_Z;
        std::vector<double>  det_h;

        bool store_A0;
        std::vector<double>  A0_values; // if store_A0 true will hold an array of central field amplitudes


        std::vector<int> refinement_points; //points at which refinement should be halved
        std::vector<bool> active_points; //determines whether each point needs to be updated; 1 if active 0 if not

        std::vector<int> refinement_levels; //lists the refinement levels of each point.

        int last_active_j;//outermost active gridpoint, for bdry purposes

        //some commonly-used derivatives that we store to avoid re-evaluating across update procedure.
        double d_z_chi ;
        double d_z_h_zz;
        double d_z_h_ww;
        double d_z_phi_re;
        double d_z_phi_im;
        double d_z_alpha;
        double d_z_beta;

        double rho0_init; //initial central field energy density

        //L2 norms of hamiltonian/momentum constraints
        double Ham_L2;
        double Mom_L2;

        int max_stored_slices; //will only hold this many previous slices in memory
        int BS_resolution_factor;
        double min_chi;

        double R; //max radius of computational domain
        int n_gridpoints; // number of spatial gridpoints
        double courant_factor;
        double stop_time;
        double dt;
        double dr;

        int start_time;
        int checkpoint_time;

        //inherited BS parameters
        double sigma;
        double mu;
        bool solitonic;
        double omega;
        bool isotropic; //whether we will use isotropic coordinates (otherwise polar-areal)
        double M;

        double sigma_BSSN; //multiplier of constraint term in the gamma evolution equation
        double eta; //gamma driver condition parameter
        double damping_factor; //multiplier for kreiss-oliger damping
        int write_interval;
        int write_CN_interval;
        double min_z; //minimum z-value below which z-> regularizations will be used
        bool evolve_shift; //if false, shift is forced to be zero
        bool run_quietly;
        bool read_thinshell;
        double cutoff_frac;

        bool make_tangherlini;

        void read_parameters(bool quiet = 0); //read auxiliary parameters


        //TODO: add checkpointing

        void auxiliary_quantities_at_point(BSSNSlice* slice_ptr, int j);
        void compute_auxiliary_quantities(BSSNSlice* slice_ptr, bool derivatives_computed = 0);
        BSSNSlice slice_rhs(BSSNSlice* slice_ptr);
        void make_A_traceless(BSSNSlice* slice_ptr);
        void compute_diagnostics(BSSNSlice* slice_ptr);
        double slice_mass(BSSNSlice* slice_ptr);
        void update_outer_boundary(double time_step);

        void fill_active_points();
        void fill_refinement_levels();

        double V(const double A);
        double dV(const double A);

        double chi_asymp(double r);
        double alpha_asymp(double r);

        double d_chi_asymp(double r);
        double d_alpha_asymp(double r);


        double d_z(bssn_var var, int index, int order); //just syntactic sugar to call d_z, d_zz on current slice
        double d_zz(bssn_var var, int index);



        double test_ctr;

};




#endif /*EVOLUTIONVARIABLES*/
