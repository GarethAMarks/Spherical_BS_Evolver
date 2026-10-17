#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include <iostream>
#include <math.h>
#include <sstream>
#include "BosonStar.h"
#include "ProcaStar.h"
#include "EvolutionVariables.h"
#include "LinearPerturbation.h"
#include "mathutils.h"
#include "mpi_helpers.h"
#include <iomanip>


using namespace std;

void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh);
void gauss_initialize(BosonStar& boson_star);
void compute_linear_perturbation(BosonStar& boson_star, double A0, double dA, double n_stars);
bool read_proca_flag();
int run_proca_star();

int main(int argc, char** argv)
{
    bss_mpi::init(&argc, &argv);
    const int mpi_rank = bss_mpi::rank();

    //if proca = 1 in BSParams.par, construct Proca stars instead of boson stars
    if (read_proca_flag())
    {
        int status = run_proca_star();
        bss_mpi::finalize();
        return status;
    }

    //construct initial BS solution
    BosonStar boson_star{};
    boson_star.read_parameters();

    if (mpi_rank == 0)
    {
        boson_star.solve();
        if (boson_star.isotropic)
        {
            boson_star.fill_isotropic_arrays();
            boson_star.write_isotropic();
        }
        boson_star.write_field();
    }
    boson_star.broadcast();

    if (boson_star.gaussian_start)
        gauss_initialize(boson_star);

    //main run: either cycles models/ radial oscillation freqs/ dynamical evolution depending on param choices
    if (boson_star.cycle_only)
        boson_star.cycle_models(boson_star.n_stars, boson_star.A0, boson_star.dA);
    else if (boson_star.pert_only)
        compute_linear_perturbation(boson_star, boson_star.A0, boson_star.dA, boson_star.n_stars);
    else
    {
        Spacetime st{};
        st.initialize(boson_star);
        // MPI: both writes are collective (all ranks gather; rank 0 writes).
        // Using write_current_slice instead of slices[0].write_slice() so ranks
        // don't all try to open the same file simultaneously.
        st.write_current_slice();
        st.write_diagnostics();
        
        if (st.critical_study)
        {
            if (st.critical_gaussian_start)
            {
                boson_star.add_real_field = false;
                double& tuning_param = boson_star.perturb_amp;
                st.tune_to_critical(tuning_param, st.hi_guess, st.lo_guess, &boson_star);
            }
            else
            {
                double& tuning_param = boson_star.real_amp; //example: tune the perturbation amplitude
                st.tune_to_critical(tuning_param, st.hi_guess, st.lo_guess, &boson_star);
            }
        }
        else
            st.evolve();
    }
    if (mpi_rank == 0)
    {
        cout <<  "c4 = " << boson_star.get_c4() << endl;
        cout << "Ending..."  << endl;
    }
    bss_mpi::finalize();
    return 0;
}

//reads the proca flag from BSParams.par (defaults to 0 if absent)
bool read_proca_flag()
{
    bool proca = 0;
    ifstream params{ "BSParams.par" };
    string current_line{};

    while (getline(params, current_line))
        fill_parameter(current_line, "proca = ", proca, 1);

    return proca;
}

//constructs a Proca star (or a family of them if cycle_only = 1). Evolution and linear perturbations are not yet implemented for Proca stars.
int run_proca_star()
{
    const int mpi_rank = bss_mpi::rank();

    ProcaStar proca_star{};
    proca_star.read_parameters();

    if (mpi_rank == 0)
    {
        proca_star.solve();
        if (proca_star.isotropic)
        {
            proca_star.fill_isotropic_arrays();
            proca_star.write_isotropic();
        }
        proca_star.write_field();
    }
    proca_star.broadcast();

    if (proca_star.cycle_only)
    {
        if (mpi_rank == 0)
            proca_star.cycle_models(proca_star.n_stars, proca_star.A0, proca_star.dA);
    }
    else if (mpi_rank == 0)
        cout << "NOTE: evolution and linear perturbations are not yet implemented for Proca stars; only the initial model was constructed (PSdata.dat)." << endl;

    if (mpi_rank == 0)
    {
        cout <<  "c4 = " << proca_star.get_c4() << endl;
        cout << "Ending..."  << endl;
    }
    return 0;
}

void gauss_initialize(BosonStar& boson_star)
{
    const int mpi_rank = bss_mpi::rank();
    std::ofstream nm_file;
    if (mpi_rank == 0) nm_file.open("mass_charge.dat");

    int num_gaussians = 1; //change this to produce a gaussian model cycle to explore 1-parameter families
    double step = 0.0001;

    double& step_var = boson_star.perturb_amp;

    for (int k = 0; k < num_gaussians; k++)
    {
        boson_star.clear_BS();
        boson_star.omega = boson_star.enforced_freq;
        boson_star.add_perturbation(boson_star.perturb_amp, boson_star.perturb_spread, 0.);

        //if (!boson_star.fill_given_A(boson_star.omega, 0) )
            boson_star.default_metric_vars();

        
        boson_star.fill_isotropic_arrays();
        boson_star.write_isotropic();
        boson_star.write_field();

        Spacetime st_gauss{};
        st_gauss.initialize(boson_star);
        st_gauss.write_current_slice();
        st_gauss.write_diagnostics();

        // slice_mass / slice_charge are collective (Allreduce); nm_file write is rank-0-only.
        const double m_local = st_gauss.slice_mass(&st_gauss.slices[0]);
        const double q_local = st_gauss.slice_charge(&st_gauss.slices[0]);
        if (mpi_rank == 0)
            nm_file << m_local << "    " << q_local <<  "    " <<  step_var << endl;

        step_var += step;
        }
}

//computes oscillation frequencies for sequence of models starting at A0 separated by dA
void compute_linear_perturbation(BosonStar& boson_star, double A0, double dA, double n_stars)
{
    LinearPerturbation lp{&boson_star, -2e-5, 0.025, 0.0001, 100.}; //s=0.08, A = 0.06: 0.00003, 0.18, 0.00003

    //lp.rk4_solve(-0.00001, 0.0235);
    //lp.get_best_gamma(-4.5e-5);
    //lp.get_chi_sq_newton();

    lp.read_parameters(0);
    //lp.get_chi_sq();
    //lp.get_best_gamma(4.5e-5);
    lp.pert_cycle(A0, dA, n_stars); //lp.pert_cycle(0.098, 0.0005, 200);
    cout << "Noether charge perturbation is " << lp.get_noether_perturbation() << endl;
    lp.write_pert();
}

/*BSSNSlice sl{}; BSSNSlice sm{}; BSSNSlice sh{};
    sl.read_BS_data(boson_star);

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sm.read_BS_data(boson_star);

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sh.read_BS_data(boson_star);

    slice_convergence_test(sl, sm, sh);*/

  /*Spacetime stl{}; Spacetime stm{}; Spacetime sth{};
    stl.initialize(boson_star);
    stl.evolve();
    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    stm.initialize(boson_star);
    stm.evolve();

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sth.initialize(boson_star);
    sth.evolve();

    slice_convergence_test(stl.slices[stl.slices.size() - 1], stm.slices[stm.slices.size() - 1], sth.slices[sth.slices.size() - 1]);*/

    //LinearPerturbation lp{&boson_star, 0.0, 0.0, 0.0001, 150.}; //s=0.08, A = 0.06: 0.00003, 0.18, 0.00003

    //FieldState f = (FieldState){0.1, 0.2, 0.3, 0.4};
    //PertState p = (PertState){0.1, 0.2, 0.3, 0.4};
    //lp.test_rhs(1.0, f, p, 0.0145, 6.3);

    //lp.rk4_solve(0.0145, 6.3); //0.00003519, 0.3097 // 0.0, 0.30965
    //lp.get_best_gamma(0.000084); //0.000083
    //lp.get_chi_sq();
    //lp.pert_cycle(0.0001, 0.0001, 150); //lp.pert_cycle(0.098, 0.0005, 200);
    //cout << "Noether charge perturbation is " << lp.get_noether_perturbation() << endl;
    //lp.write_pert();
    //lp.write_chi_results();


    //st.slices[st.slices.size() - 1].write_slice();
