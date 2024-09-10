#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include <iostream>
#include <math.h>
#include <sstream>
//#include "DimensionMacros.h"
#include "BosonStar.h"
#include "EvolutionVariables.h"
#include "mathutils.h" //no idea why this is necessary. gcc is screwy
#include <iomanip>


using namespace std;

void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh);

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

    //boson_star.rk4_solve(2.85); //2.799473
    //boson_star.write_field();
    //boson_star.count_zero_crossings();
    //boson_star.convergence_test(2.799473);


    boson_star.solve();
    //boson_star.rk4_solve(50.0);
    boson_star.write_field();
    boson_star.fill_isotropic_arrays();
    boson_star.write_isotropic();

    //boson_star.cycle_models(500, 0.01, 0.0005);
    //boson_star.write_field();

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

    //boson_star.convergence_test();

    Spacetime st{};
    st.initialize(boson_star);

    //st.slices[0].make_tangherlini(1.);

    st.slices[0].write_slice();
    st.write_diagnostics();

    st.evolve();
    st.fourier_transform_A0();


    //st.slices[st.slices.size() - 1].write_slice();

    /*BSSNSlice slice{};
    slice.read_BS_data(boson_star);
    slice.write_slice();*/



    cout << "Ending..."  << endl;

    return 0;
}


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

