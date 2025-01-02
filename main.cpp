#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include <iostream>
#include <math.h>
#include <sstream>
#include "BosonStar.h"
#include "EvolutionVariables.h"
#include "mathutils.h"
#include <iomanip>


using namespace std;

void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh);
void gauss_initialize(BosonStar& boson_star);

enum evolution_system
{
    BSSN,
    CCZ4
};


int main()
{
    //enum evolution_system system = BSSN;

    BosonStar boson_star{};
    boson_star.read_parameters();


    //cout << boson_star.solve() << endl;
    //boson_star.rk4_solve(50.0);
    //boson_star.write_field();
    //boson_star.fill_isotropic_arrays();a
    //boson_star.write_isotropic();

    //boson_star.read_thinshell();
    boson_star.solve();
    boson_star.fill_isotropic_arrays();
    boson_star.write_isotropic();
    boson_star.write_field();

    //boson_star.cycle_models(5000, 0.001, 0.001);

    if (boson_star.gaussian_start)
        gauss_initialize(boson_star);

    //boson_star.convergence_test();

    Spacetime st{};
    st.initialize(boson_star);
    st.slices[0].write_slice();
    st.write_diagnostics();

    st.evolve();

    //st.fourier_transform_A0();

    //st.slices[st.slices.size() - 1].write_slice();

    cout << "Ending..."  << endl;

    return 0;
}

void gauss_initialize(BosonStar& boson_star)
{

    //TODO NEXT: add proper logic for when to add perturbations + maybe start with random mini BS model chi + alpha or similar

    std::ofstream nm_file{"mass_charge.dat"};

    int num_gaussians = 1; //change this to produce a gaussian model cycle to explore 1-parameter families
    double step = 0.0001;

    double& step_var = boson_star.perturb_amp;

    for (int k = 0; k < num_gaussians; k++)
    {
        boson_star.clear_BS();
        boson_star.omega = boson_star.enforced_freq;
        boson_star.add_perturbation(boson_star.perturb_amp, boson_star.perturb_spread, 0.);

        if (!boson_star.fill_given_A(boson_star.omega, 0) )
        {
            //boson_star.solitonic = 0;
            //boson_star.solve();
           //boson_star

           boson_star.default_metric_vars();
        }
        boson_star.fill_isotropic_arrays();
        boson_star.write_isotropic();
        boson_star.write_field();

        Spacetime st_gauss{};
        st_gauss.initialize(boson_star);
        st_gauss.slices[0].write_slice();
        st_gauss.write_diagnostics();

        nm_file << st_gauss.slice_mass(&st_gauss.slices[0]) << "    " << st_gauss.slice_charge(&st_gauss.slices[0]) <<  "    " <<  step_var << endl;

        step_var += step;
        }

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

