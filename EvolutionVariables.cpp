#ifndef EVOLUTIONVARIABLES_CPP_
#define EVOLUTIONVARIABLES_CPP_

#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include "EvolutionVariables.h"
#include "mathutils.h"
#include "spline.h"
#include <iomanip>
#include<algorithm>
#include<filesystem>
#include <complex.h>
#include <array>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::rotate;
// Spacetime constructor: initialize reusable CSF model with this context
Spacetime::Spacetime() : csf_model(this), rsf_model(this) {}

//overloads for addition/ scalar multiplication of BSSNState sets and the slice arrays containing them
// Forward declarations for BSSNState operators used by State operators below
BSSNState operator+(const BSSNState& s1, const BSSNState& s2);
BSSNState operator-(const BSSNState& s1, const BSSNState& s2);
BSSNState operator*(const BSSNState& s1, const BSSNState& s2);
BSSNState operator*(double c, const BSSNState& s);

// Composite State elementwise operations (applies to both csf and bssn)
State operator+(const State& a, const State& b) {
    return State{a.csf + b.csf, a.bssn + b.bssn, a.rsf + b.rsf};
}

State operator-(const State& a, const State& b) {
    return State{a.csf - b.csf, a.bssn - b.bssn, a.rsf - b.rsf};
}

State operator*(const State& a, const State& b) {
    return State{a.csf * b.csf, a.bssn * b.bssn, a.rsf * b.rsf};
}

State operator*(double c, const State& a) {
    return State{c * a.csf, c * a.bssn, c * a.rsf};
}

inline State operator/(const State& a, double c) {
    return (1.0 / c) * a;
}

BSSNState operator+(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi + s2.chi, s1.h_zz + s2.h_zz, s1.h_ww + s2.h_ww, s1.A_zz + s2.A_zz, s1.A_ww + s2.A_ww, s1.K + s2.K, s1.c_chris_Z + s2.c_chris_Z,
                      s1.alpha + s2.alpha, s1.beta + s2.beta};
}

BSSNState operator-(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi - s2.chi, s1.h_zz - s2.h_zz, s1.h_ww - s2.h_ww, s1.A_zz - s2.A_zz, s1.A_ww - s2.A_ww, s1.K - s2.K, s1.c_chris_Z - s2.c_chris_Z,
                      s1.alpha - s2.alpha, s1.beta - s2.beta};
}

//termwise multiplication for convenient use with characteristic speeds etc.
BSSNState operator*(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi * s2.chi, s1.h_zz * s2.h_zz, s1.h_ww * s2.h_ww, s1.A_zz * s2.A_zz, s1.A_ww * s2.A_ww, s1.K * s2.K, s1.c_chris_Z * s2.c_chris_Z,
                      s1.alpha * s2.alpha, s1.beta * s2.beta};
}

BSSNState operator*(double c, const BSSNState& s)
{
    return (BSSNState){c * s.chi, c * s.h_zz, c * s.h_ww , c * s.A_zz, c * s.A_ww, c * s.K, c * s.c_chris_Z,
                      c * s.alpha, c * s.beta};
}

inline BSSNState operator/(const BSSNState& s, double c)
{
    return (1. / c) * s;
}

BSSNSlice operator+(const BSSNSlice& slice1, const BSSNSlice& slice2)
{
    BSSNSlice return_slice;
    unsigned int length = slice1.states2.size();

    if (length != slice2.states2.size() || slice1.R != slice2.R)
    {
        cerr << "ERROR: attempted to add data from slices of different sizes!" << endl;
        exit(1);
    }

    return_slice.states2.resize(length);

    if (slice1.use_CCZ4)
        return_slice.theta.resize(length);


    return_slice.R = slice1.R;
    return_slice.has_BH = slice1.has_BH;
    return_slice.use_CCZ4 = slice1.use_CCZ4;
    return_slice.refinement_points = slice1.refinement_points;

    for (unsigned int j = 0; j < length; j++)
    {
        return_slice.states2[j] = slice1.states2[j] + slice2.states2[j];

        if (slice1.use_CCZ4)
            return_slice.theta[j] = slice1.theta[j] + slice2.theta[j];
    }

    return return_slice;
}

BSSNSlice operator*(double c, const BSSNSlice& slice)
{
    BSSNSlice return_slice;

    int length = slice.states2.size();
    return_slice.states2.resize(length);

    if (slice.use_CCZ4)
        return_slice.theta.resize(length);

    return_slice.R = slice.R;
    return_slice.has_BH = slice.has_BH;
    return_slice.refinement_points = slice.refinement_points;
    return_slice.use_CCZ4 = slice.use_CCZ4;

    for (int j = 0; j < length; j++)
    {
        return_slice.states2[j] = c * slice.states2[j];

        if (slice.use_CCZ4)
            return_slice.theta[j] = c * slice.theta[j];
    }

    return return_slice;
}

//radial partial derivative of a given bssn var at index. Order is an optional argument that allows higher z-derivatives to be taken, default is 1.
//should add chacks on refinement levels
double BSSNSlice::d_z(bssn_var var, int index, int order = 1 )
{
    // Use composite state vector size for grid length to ensure compatibility with RK temp slices
    int n_gridpoints = states2.size();
    // Use cell-centred grid spacing to match read_BS_data (innermost point at z = 0.5*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    int ref_level = get_refinement_level(index, refinement_points); //refinement level; starts at 1 for no refinement and halves every time.

    // Use bit-shift instead of pow for 2^(ref_level-1)
    int res_fac = 1 << (ref_level - 1);

    //check index is valid and error out if not
   if (index < 0 || index >= n_gridpoints )
    {
        cerr << "ERROR: invalid index requested in radial derivative" << endl;
        exit(1);
    }

    //set of indices using which derivative will be evaluated
    vector<int> J{index - 2 * res_fac, index - res_fac, index, index + res_fac, index + 2 * res_fac};

    //enforce BC at z = 0 by using symmetry
    if (index <= 1)
        J[0] = -J[0];
    if (index == 0)
        J[1] = -J[1];


    //At outer edge just use central value for rightmost two; this is cut off in initialize()/ overwritten with BCs anyway
    if (index >= n_gridpoints - res_fac)
        J[3] = J[2];
    if (index >= n_gridpoints - 2 * res_fac)
        J[4] = J[3];

    //set of values used to compute derivative
    vector<double> F{0., 0., 0., 0., 0.};

    //Fill out five values; note that at outer boundary this gives unreliable results
    if (index < n_gridpoints)
    {
        for (int j  = 0; j < 5; j++)
        {
            switch (var)
            {
                case v_chi:
                    F[j] = states2[J[j]].bssn.chi;
                    break;
                case v_h_zz:
                    F[j] = states2[J[j]].bssn.h_zz;
                    break;
                case v_h_ww:
                    F[j] = states2[J[j]].bssn.h_ww;
                    break;
                case v_A_zz:
                    F[j] = states2[J[j]].bssn.A_zz;
                    break;
                case v_A_ww:
                    F[j] = states2[J[j]].bssn.A_ww;
                    break;
                case v_K:
                    F[j] = states2[J[j]].bssn.K;
                    break;
                case v_c_chris_Z:
                    F[j] = states2[J[j]].bssn.c_chris_Z;
                    break;
                case v_phi_re:
                    F[j] = states2[J[j]].csf.phi_re;
                    break;
                case v_phi_im:
                    F[j] = states2[J[j]].csf.phi_im;
                    break;
                case v_K_phi_re:
                    F[j] = states2[J[j]].csf.K_phi_re;
                    break;
                case v_K_phi_im:
                    F[j] = states2[J[j]].csf.K_phi_im;
                    break;
                case v_alpha:
                    F[j] = states2[J[j]].bssn.alpha;
                    break;
                case v_beta:
                    F[j] = states2[J[j]].bssn.beta;
                    break;
                case v_theta:
                    F[j] = theta[J[j]];
                    break;
                case v_psi:
                    F[j] = states2[J[j]].rsf.psi;
                    break;
                case v_K_psi:
                    F[j] = states2[J[j]].rsf.K_psi;
                    break;
                default:
                    cerr << "ERROR: invalid variable requested for differentiation" << endl;
                    exit(0);
            }
        }

        bool parity_is_odd = ((var == v_c_chris_Z || var == v_beta)/*|| (has_BH && (var == v_chi || var == v_alpha))*/);
       //if (var == v_alpha) cout << parity_is_odd<< endl;

        //account for odd parity of contracted christoffel symbols and beta across z = 0
            if (index <= 1 && parity_is_odd )
                F[0] = -F[0];
            if (index == 0 && parity_is_odd)
                F[1] = -F[1];
    }
    return fivePointDeriv(res_fac * dr, order, F[0],F[1],F[2],F[3],F[4]);
}

//second z-derivative, basically syntactic sugar for d_z(var, index, 2)
double BSSNSlice::d_zz(bssn_var var, int index )
{
    return d_z(var, index, 2);
}

//converts the current slice to a tangherlini BH of mass m; must have set states vector size already. Also returns mass so it can be set in spacetime object
double BSSNSlice::make_tangherlini (double m, double min_chi, double D)
{
    int n_gridpoints = states2.size();

    // Use cell-centred grid spacing to match slice convention (z = (j+0.5)*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    double psi_power = -4. / (D - 3.); //12. / ((D - 3.) * (1. - D));

    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = ( (j + 0.5) * dr);

        //conformal factor in isotropic convention
        double psi = 1. + pow(m, D - 3.) / (4 * pow(r, D - 3.));

        //conformal factor and conformally rescaled metric components are all powers of X in this gauge
    states2[j].bssn.chi = max(min_chi, pow(psi, psi_power));
    states2[j].bssn.h_zz = 1.;
    states2[j].bssn.h_ww = 1.;

    states2[j].bssn.A_zz = 0.;
    states2[j].bssn.A_ww = 0.;
    states2[j].bssn.K = 0.;
    states2[j].bssn.c_chris_Z = 0.;
    states2[j].csf.phi_re = 0.;
    states2[j].csf.phi_im = 0.;
    states2[j].csf.K_phi_re = 0.;
    states2[j].csf.K_phi_im = 0.;

    states2[j].bssn.alpha = pow(states2[j].bssn.chi, 0.5);
    states2[j].bssn.beta = 0.;

    }
    return m;
}

// Helper to apply high-order damping at a single gridpoint j
void Spacetime::apply_damping_at_point(int j,
                                       BSSNSlice* slice_ptr,
                                       std::vector<State>& rhs_states2,
                                       std::vector<double>& rhs_theta)
{
    // refinement factor local to j
    int res_fac = 1 << (slice_ptr->get_refinement_level(j, refinement_points) - 1);

    // indices for 7-point centered stencil with symmetry fix near the origin
    std::array<int, 7> J = {j - 3 * res_fac, j - 2 * res_fac, j - 1 * res_fac, j, j + 1 * res_fac, j + 2 * res_fac, j + 3 * res_fac};
    if (j < 3) {
        for (int k = 0; k < 7; ++k) J[k] = std::max(-J[k], J[k]);
    }

    const double h = dr * res_fac;
    const double d_mult = damping_factor * ((h * h * h * h * h) / 64.0);

    // references for State stencil
    State& sJ0 = slice_ptr->states2[J[0]];
    State& sJ1 = slice_ptr->states2[J[1]];
    State& sJ2 = slice_ptr->states2[J[2]];
    State& sJ3 = slice_ptr->states2[J[3]];
    State& sJ4 = slice_ptr->states2[J[4]];
    State& sJ5 = slice_ptr->states2[J[5]];
    State& sJ6 = slice_ptr->states2[J[6]];

    State damping_corr = d_mult * sevenPointDeriv(h, 6, sJ0, sJ1, sJ2, sJ3, sJ4, sJ5, sJ6);
    rhs_states2[j] = rhs_states2[j] + damping_corr;

    if (use_CCZ4 && rhs_theta.size() == static_cast<size_t>(n_gridpoints))
    {
        double& tJ0 = slice_ptr->theta[J[0]];
        double& tJ1 = slice_ptr->theta[J[1]];
        double& tJ2 = slice_ptr->theta[J[2]];
        double& tJ3 = slice_ptr->theta[J[3]];
        double& tJ4 = slice_ptr->theta[J[4]];
        double& tJ5 = slice_ptr->theta[J[5]];
        double& tJ6 = slice_ptr->theta[J[6]];
        double theta_corr = d_mult * sevenPointDeriv<double>(h, 6, tJ0, tJ1, tJ2, tJ3, tJ4, tJ5, tJ6);
        rhs_theta[j] += theta_corr;
    }
}

// Apply Sommerfeld-like outer boundary conditions to rhs
// Updates rhs.states2 at [last_active_j] and [last_active_j - res_fac]
void Spacetime::sommerfeld_BC(BSSNSlice& rhs, int res_fac_outer, BSSNSlice* slice_ptr)
{
    int res_fac = res_fac_outer; //number of points skipped per slot in stencil

    // asymptotic states and their derivatives where relevant
    State asymp_state{}; // matter part zero by default
    asymp_state.bssn = (BSSNState){chi_asymp(R - dr * res_fac), 1., 1., 0., 0., 0., 0., alpha_asymp(R - dr * res_fac), 0.};
    State asymp_deriv{}; // r-derivative of asymptotic expansion (matter zero)
    asymp_deriv.bssn = (BSSNState){d_chi_asymp(R - dr * res_fac), 0., 0., 0., 0., 0., 0., d_alpha_asymp(R - dr * res_fac), 0.};

    //version without spatially varying asymptotics
    if (!spatially_varying_BC)
    {
        asymp_state.bssn = (BSSNState){1., 1., 1., 0., 0., 0., 0., 1., 0.};
        asymp_deriv.bssn = (BSSNState){0., 0., 0., 0., 0., 0., 0., 0., 0.};
    }

    //maybe adjust to account for purported 1/r^2 decay in K!
    State N{};
    N.bssn = (BSSNState){1., 1., 1., 1., 1., 1., 1., 1., 1.};
    N.csf = (CSF){1., 1., 1., 1.};
    
    if (add_real_field)
        N.rsf = (RSF){1., 1.};

    N = pow(R, (4. - D))  * N;

    //outermost 5 active gridpoints
    State& s1 = slice_ptr->states2[last_active_j - 4 * res_fac];
    State& s2 = slice_ptr->states2[last_active_j - 3 * res_fac];
    State& s3 = slice_ptr->states2[last_active_j - 2 * res_fac];
    State& s4 = slice_ptr->states2[last_active_j - 1 * res_fac];
    State& s5 = slice_ptr->states2[last_active_j];

    //limiting characteristic speeds at infinity
    // characteristic speeds
    State char_speeds{};
    double cs = sqrt(2. / s4.bssn.alpha);
    char_speeds.bssn = (BSSNState){cs, 1., 1., 1., 1., cs, 1., cs, 1.};
    char_speeds.csf = (CSF){1., 1., 1., 1.};

    if (add_real_field)
        char_speeds.rsf = (RSF){1., 1.};

    rhs.states2[last_active_j - res_fac] = (-1.) * char_speeds * ( p4_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                         + (0.5 * D - 1.) * N * (s4 - asymp_state) /  (dr * (last_active_j - res_fac)) - asymp_deriv);

    //update variable asymptotic states to outermost edge
    if (spatially_varying_BC)
    {asymp_state.bssn.chi = chi_asymp(R);  asymp_state.bssn.alpha = alpha_asymp(R);
    asymp_deriv.bssn.chi = d_chi_asymp(R);  asymp_deriv.bssn.alpha = d_alpha_asymp(R);}

    rhs.states2[last_active_j] = (-1) * char_speeds * ( p5_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                 + (0.5 * D - 1.) * N * (s5 - asymp_state) / (dr * last_active_j) - asymp_deriv);

    if (use_CCZ4)
    {
        double& t1 = slice_ptr->theta[last_active_j - 4 * res_fac];
        double& t2 = slice_ptr->theta[last_active_j - 3 * res_fac];
        double& t3 = slice_ptr->theta[last_active_j - 2 * res_fac];
        double& t4 = slice_ptr->theta[last_active_j - 1 * res_fac];
        double& t5 = slice_ptr->theta[last_active_j];

        rhs.theta[last_active_j - res_fac] = -(p4_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                             + (0.5 * D - 1.) * t4 / (dr * (last_active_j - res_fac)));

        rhs.theta[last_active_j] = -(p5_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                   + (0.5 * D - 1.) * t5 / (dr * (last_active_j)));
    }
}

//enforce tracelessness of A
void Spacetime::make_A_traceless(BSSNSlice* slice_ptr)
{
    double n = D - 2.;

    for (int j = 0; j < n_gridpoints; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

        const double h_zz = slice_ptr->states2[j].bssn.h_zz;
        const double h_ww = slice_ptr->states2[j].bssn.h_ww;
        const double A_zz = slice_ptr->states2[j].bssn.A_zz;
        const double A_ww = slice_ptr->states2[j].bssn.A_ww;

            double A = A_zz / h_zz + n * A_ww / h_ww;
        slice_ptr->states2[j].bssn.A_zz = A_zz - h_zz * A / (D - 1.);
        slice_ptr->states2[j].bssn.A_ww = A_ww - h_ww * A / (D - 1.);
    }

}


//reads data from a boson_star object we have already solved for into initial BSSN slice.
void BSSNSlice::read_BS_data (BosonStar& boson_star, int BS_resolution_factor, bool isotropic, bool cell_centered)
{
    //prepare to fill states array

    if ((boson_star.n_gridpoints + BS_resolution_factor - 1) % BS_resolution_factor != 0 )
        cout << "WARNING: Incompatible resolution factor used in read_BS_data" << endl;

    int n_gridpoints = (boson_star.n_gridpoints + BS_resolution_factor - 1) / BS_resolution_factor;
    states2.resize(n_gridpoints); // initialize composite State vector

    if (use_CCZ4)
       theta.resize(n_gridpoints);

    R = boson_star.R;
    // Grid spacing; when cell_centered==true, cell centers are at z=(j+0.5)*dr
    double dr = R / ((double)n_gridpoints - 1.0);

    double D = boson_star.D;

    //
    // Precompute BosonStar spacing and sizes for interpolation
    int N_bs = boson_star.n_gridpoints;
    double dr_bs = boson_star.R / (double)(max(1, N_bs - 1));

    // Build x-grid for spline (uniform nodes j*dr_bs)
    std::vector<double> X_nodes(N_bs);
    for (int k = 0; k < N_bs; ++k) X_nodes[k] = k * dr_bs;

    // Build contiguous arrays and their splines
    tk::spline s_X, s_A, s_phi, s_psi_iso, s_A_iso, s_phi_iso, s_pert, s_pert_iso;
    if (isotropic) {
        if (!boson_star.psi_iso_array.empty()) {
            s_psi_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_psi_iso.set_points(X_nodes, boson_star.psi_iso_array, tk::spline::cspline);
        }
        if (!boson_star.A_iso_array.empty())   {
            s_A_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_A_iso.set_points(X_nodes, boson_star.A_iso_array, tk::spline::cspline);
        }
        if (!boson_star.phi_iso_array.empty()) {
            s_phi_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_phi_iso.set_points(X_nodes, boson_star.phi_iso_array, tk::spline::cspline);
        }
        if (boson_star.perturb && !boson_star.pert_iso_array.empty()) {
            s_pert_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_pert_iso.set_points(X_nodes, boson_star.pert_iso_array, tk::spline::cspline);
        }
    } else {
        std::vector<double> bs_X(N_bs), bs_A(N_bs), bs_phi(N_bs);
        for (int k = 0; k < N_bs; ++k) {
            bs_X[k] = boson_star.state[k].X;
            bs_A[k] = boson_star.state[k].A;
            bs_phi[k] = boson_star.state[k].phi;
        }
        s_X.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_A.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_phi.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_X.set_points(X_nodes, bs_X, tk::spline::cspline);
        s_A.set_points(X_nodes, bs_A, tk::spline::cspline);
        s_phi.set_points(X_nodes, bs_phi, tk::spline::cspline);
        if (boson_star.perturb && !boson_star.pert_array.empty()) {
            s_pert.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_pert.set_points(X_nodes, boson_star.pert_array, tk::spline::cspline);
        }
    }

    for (int j = 0; j < n_gridpoints; j++)
    {
        // cell- or vertex-centered coordinate
        double z = (j + 0.5 * (double)cell_centered) * dr;

        // Interpolate X (areal conformal factor) or psi (isotropic conformal factor)
        double X_interp = 1.;
        double psi_interp = 1.;
        if (isotropic)
        {
            // interpolate isotropic psi (conformal factor)
            psi_interp = s_psi_iso(z);
            states2[j].bssn.chi = pow(psi_interp, -4.);
            states2[j].bssn.h_zz = 1.; //isotropic gauge is conformally flat
            states2[j].bssn.h_ww = 1.;
        }
        else
        {
            X_interp = s_X(z);
            states2[j].bssn.chi = pow(X_interp, -2. / (D - 1.));
            states2[j].bssn.h_zz = pow(X_interp, (2. * D - 4.) / (D - 1.));
            states2[j].bssn.h_ww = pow(X_interp, -2. / (D - 1.));
        }

        //for time-independent static, spherically symmetric boson stars K_ij = 0-- may need to do more work when we consider perturbations.
        states2[j].bssn.A_zz = 0.;
        states2[j].bssn.A_ww = 0.;
        states2[j].bssn.K = 0.;

        // start at t = 0 so the scalar field is real regardless of phase
        double A_interp = 0.;
        double phi_interp = 0.;
        if (isotropic)
        {
            A_interp = s_A_iso(z);
            phi_interp = s_phi_iso(z);
            states2[j].csf.phi_re = A_interp;
            states2[j].csf.phi_im = 0.;
            states2[j].bssn.alpha = exp(phi_interp);
            states2[j].bssn.beta = 0.;
        }
        else
        {
            A_interp = s_A(z);
            phi_interp = s_phi(z);
            states2[j].csf.phi_re = A_interp;
            states2[j].csf.phi_im = 0.;
            states2[j].bssn.alpha = exp(phi_interp);
            states2[j].bssn.beta = 0.;
        }

        //starting BS real means its momentum is imaginary (with 0 starting shift)
    states2[j].csf.K_phi_re = 0.;

        double pert_correction = 0.;
        if (boson_star.perturb)
        {
            double pert_interp = 0.;
            if (isotropic)
                pert_interp = s_pert_iso(z);
            else
                pert_interp = s_pert(z);

            pert_correction = 1. * (boson_star.omega / states2[j].bssn.alpha) * pert_interp; //conserves Noether charge at 1st order in perturbation
        }

    states2[j].csf.K_phi_im = - boson_star.omega * states2[j].csf.phi_re / (2. * states2[j].bssn.alpha) + pert_correction;

        if (use_CCZ4)
            theta[j] = 0.;
    }

    //separate loop to fill in values for the contracted Christoffel symbols as they require derivatives of h (hence info off j value)
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j) * dr;

        // inverse metric components, for convenience
        double h_ZZ = 1 / states2[j].bssn.h_zz;
        double h_WW = 1 / states2[j].bssn.h_ww;

        states2[j].bssn.c_chris_Z = 0.5 * h_ZZ * h_ZZ * d_z(v_h_zz,j) + (D - 2.) * (-0.5 * h_ZZ * h_WW * d_z(v_h_ww,j));

        // extra geometric term, which goes to zero at z = 0
        if (z > 1e-8)
            states2[j].bssn.c_chris_Z += (D - 2.) * (h_WW * (1 - h_ZZ * states2[j].bssn.h_ww) / z);

        // hard enforce 0 for isotropic data
        if (isotropic)
            states2[j].bssn.c_chris_Z = 0.;
    }

    // Legacy states removed: no mirroring
}

//writes BSSN evolution variables from a particular slice to text file
void BSSNSlice::write_slice(std::string file_name)
{
    int length = states2.size();
    double dr = R / ((double)length - 1.0);
    std::ofstream data_file{file_name};

     if (!data_file)
    {
        cerr << "SliceData.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 0; j < length; j++)
    {
        //if (!active_points[j]) //perform no computations on inactive points
            //continue;

        double theta0 = 0.; //if using CCZ4, write theta at the end
        if (use_CCZ4)
            theta0 = theta[j];

    data_file <<  std::setprecision (16) << (j) * dr << "   " << states2[j].bssn.chi << "    " << states2[j].bssn.h_zz << "    " << states2[j].bssn.h_ww  << "    " << states2[j].bssn.A_zz
        << "   " << states2[j].bssn.A_ww << "    " << states2[j].bssn.K << "    " << states2[j].bssn.c_chris_Z  << "    " << states2[j].csf.phi_re << "    " << states2[j].csf.phi_im
        << "   " << states2[j].csf.K_phi_re << "    " << states2[j].csf.K_phi_im << "    " << states2[j].bssn.alpha << "    " << states2[j].bssn.beta << "    " <<  theta0
        << "   " << states2[j].rsf.psi << "    " << states2[j].rsf.K_psi << endl;
    }
}

//reads slice data from checkpoint file of form checkpoint(time).dat w/o parentheses
void BSSNSlice::read_checkpoint(int time, int n_gridpoints)
{
    states2.resize(n_gridpoints);

    std::ifstream checkpoint_file("checkpoint" + std::to_string(time) + ".dat");
    if (!checkpoint_file.is_open())
    {
        cerr << "Could not open checkpoint file!" << endl;
        exit(1);
    }

    string line;
    int j = 0;
    while (std::getline(checkpoint_file, line))
    {
        std::istringstream iss(line);
    double z;
    // read into temporaries to populate both legacy states and new states2
    double chi,h_zz,h_ww,A_zz,A_ww,K,c_chris_Z,phi_re,phi_im,K_phi_re,K_phi_im,alpha,beta;

        if (use_CCZ4)
        {
            double& t = theta[j];

            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta >> t)
            {
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                j++;
            }
        }
        else
        {
            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta)
            {
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                j++;
            }
        }
    }
}

//determines which refinement level j belongs to. 1 is highest resolution; for each level up the resolution halves.
int BSSNSlice::get_refinement_level(int j, std::vector<int>& refinement_points)
{
    //cout << "Size is " <<  refinement_points.size() << endl;

    if (refinement_points.empty()) //always return 1 if no refinement
        return 1;

    int n_refinements = refinement_points.size();
    int level = 1;
    int k = 0;

    while ( k < n_refinements && j >= refinement_points[k] - pow(2, k))//check this difference -- meant to ensure we can use a stencil at points spaced by 2^(k + 1) safely at j
    {
        level++;
        k++;
    }

    return level;
}

//convergence test for three slices. Resolutions must differ by factor of 2. Also, currently does not support refinement.
void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh)
{
    if (sl.R != sm.R || sm.R != sh.R)
    {
        cout << "ERROR: Convergence test requested on slices of different radii";
        exit(1);
    }

    if (2*sl.states2.size() - 1 != sm.states2.size()  || 2*sm.states2.size() - 1 != sh.states2.size()   )
    {
        cout << "ERROR: Convergence test requested on slices with incompatible sizes";
        exit(1);
    }

    //write to file
    std::ofstream conv_file{ "conv.dat" };

    // If we couldn't open the output file stream for writing
    if (!conv_file)
    {
        // Print an error and exit
        cerr << "Error: conv.dat could not be opened for writing\n";
        exit(1);
    }

    //write change in A between med/low and high/med resolution to file.
    for (unsigned int j = 0; j < sl.states2.size(); j++)
    {
    // For 4th-order convergence the refinement factor scaling is 2^4 = 16
    conv_file << j * sl.R /(sl.states2.size() - 1) << "   " << sm.states2[2*j].csf.phi_re -  sl.states2[j].csf.phi_re  << "    " << 16.*(sh.states2[4*j].csf.phi_re  -  sm.states2[2*j].csf.phi_re ) << endl;
    }
}


//these let us call d_z and d_zz on the current slice without explicitly referencing it
//must remember to update current_slice_ptr appropriately!!!
double Spacetime::d_z(bssn_var var, int index, int order = 1)
{
    return current_slice_ptr->d_z(var, index, order);
}

double Spacetime::d_zz(bssn_var var, int index)
{
    return current_slice_ptr->d_zz(var, index);
}

double Spacetime::V(const double A)
{
    if (!solitonic)
        return mu * mu * A * A + lambda * pow(A,4);

    else
        return mu * mu * A * A * pow((1. - 2. * pow(A / sigma, 2)), 2) + lambda * pow(A,4) ;
}

double Spacetime::dV(const double A)
{
    if (!solitonic)
        return mu * mu + 2 * A * A * lambda;

    else
        return mu * mu - 8. * mu * mu * pow(A / sigma, 2) + 12. * mu * mu * pow(A / sigma, 4) + 2. * A * A * lambda;
}

//asymptotic values of chi and alpha and their derivatives. TODO: generalize D > 4
double Spacetime::chi_asymp(double r)
{
    return(isotropic) ? (pow(1. +  0.5 * M / r, -4.)) : pow(1 - 2. * M / r, -1./3.); //use isotropic/areal Schwarzchilld value for chi as appropriate
}

double Spacetime::alpha_asymp(double r)
{
        return sqrt(1 - 2. * M / R);
}

double Spacetime::d_chi_asymp(double r)
{
    return(isotropic) ? (2 * M * pow(1. +  0.5 * M / r, -5.) / (r * r)) : -(2./3.) * M * pow(1 - 2. * M / r, -4./3.) / (r * r);
}

double Spacetime::d_alpha_asymp(double r)
{
        return (M / (R * R)) / sqrt(1 - 2. * M / r);
}


 // returns RHS as {chi, eta} where eta is the radial derivative of chi
 vector<double> Spacetime::ham_init_rhs(double r, double chi, double eta)
 {
    int k = floor(r / dr) - 1;
    while (k * dr < r)
        k++;

    const BSSNSlice& s = slices[0]; //initial slice
    current_slice_ptr = &slices[0];

    int j0 = bound(k - 1, 0, n_gridpoints - 4);

    //cubic interpolation to get state values off gridpoints (for midstep solving)
    double phi_re = cubic_interp(r, s.states2[j0].csf.phi_re, s.states2[j0 + 1].csf.phi_re, s.states2[j0 + 2].csf.phi_re, s.states2[j0 + 3].csf.phi_re, j0, dr);
    double phi_im = cubic_interp(r, s.states2[j0].csf.phi_im, s.states2[j0 + 1].csf.phi_im, s.states2[j0 + 2].csf.phi_im,s.states2[j0 + 3].csf.phi_im, j0, dr);
    double K_phi_re = cubic_interp(r, s.states2[j0].csf.K_phi_re, s.states2[j0 + 1].csf.K_phi_re, s.states2[j0 + 2].csf.K_phi_re, s.states2[j0 + 3].csf.K_phi_re,j0, dr);
    double K_phi_im = cubic_interp(r, s.states2[j0].csf.K_phi_im, s.states2[j0 + 1].csf.K_phi_im, s.states2[j0 + 2].csf.K_phi_im, s.states2[j0 + 3].csf.K_phi_im,j0, dr);

    //WARNING: these are very similarly named to Spacetime member variables introduced to avoid re-computing derivatives during evolution; do not confuse!
    double dz_phi_re = cubic_interp(r, d_z(v_phi_re, j0), d_z(v_phi_re, j0 + 1), d_z(v_phi_re, j0 + 2), d_z(v_phi_re, j0 + 3),j0,  dr);
    double dz_phi_im = cubic_interp(r, d_z(v_phi_im, j0), d_z(v_phi_im, j0 + 1), d_z(v_phi_im, j0 + 2), d_z(v_phi_im, j0 + 3),j0,  dr);

    double mod_phi = sqrt(phi_re * phi_re + phi_im + phi_im);
    double rho0 = 2. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) +  0.5 * chi * (pow(dz_phi_re,2) + pow(dz_phi_im,2)) + 0.5 * V(mod_phi);
    
    //account for real scalar field if needed
    if (add_real_field)
    {   
        //double psi = cubic_interp(r, s.states2[j0].rsf.psi, s.states2[j0 + 1].rsf.psi, s.states2[j0 + 2].rsf.psi, s.states2[j0 + 3].rsf.psi,j0, dr);
        double dz_psi= cubic_interp(r, d_z(v_psi, j0), d_z(v_psi, j0 + 1), d_z(v_psi, j0 + 2), d_z(v_psi, j0 + 3),j0,  dr);
        double K_psi = cubic_interp(r, s.states2[j0].rsf.K_psi, s.states2[j0 + 1].rsf.K_psi, s.states2[j0 + 2].rsf.K_psi, s.states2[j0 + 3].rsf.K_psi,j0, dr);

        rho0 += 0.5 * chi * dz_psi * dz_psi + 2. * K_psi * K_psi;
    }

    double chi_rhs = eta;
    double eta_rhs = 5. * eta * eta / (4. * chi) + 8. * M_PI * rho0;

    if (r > 0.) eta_rhs -= 2 * eta / r;
    else eta_rhs /= 3.; // in the limit r -> 0, eta / r is replaced by d_r(eta), so net effect is to divide the usual RHS by 3!

    vector<double> rhs = {chi_rhs, eta_rhs};
    return rhs;

 }

 void Spacetime::solve_initial_ham(bool quiet)
 {
    if (!isotropic)
    {
        cout << "WARNING: Spacetime Hamiltonian solver called in non-isotropic coordinates; aborting request" << endl;
        return;
    }

    BSSNSlice& s = slices[0]; //initial slice
    double c1, c2, c3, c4, e1, e2, e3, e4;
    double eta = 0;

    if (!quiet)
        cout << "\n Starting initial Spacetime Hamiltonian constraint solver..." << endl;

   for (int j = 0; j < n_gridpoints - 1; j++) //for (int j = n_gridpoints - 1; j > 0; j--)
    {
        double r = j * dr;

        //RK4 solver for Ham. constraint
        //pretty inefficient but OK for now
        c1 = ham_init_rhs(r, s.states2[j].bssn.chi, eta)[0];
        e1 = ham_init_rhs(r, s.states2[j].bssn.chi, eta)[1];

        c2 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[0];
        e2 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[1];

        c3 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[0];
        e3 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[1];

        c4 = ham_init_rhs(r + dr, s.states2[j].bssn.chi + dr * c3, eta + dr * e3)[0];
        e4 = ham_init_rhs(r + dr, s.states2[j].bssn.chi + dr * c3, eta + dr * e3)[1];

        s.states2[j + 1].bssn.chi = s.states2[j].bssn.chi + (dr / 6.) * (c1 + 2 * c2 + 2 * c3 + c4);
        eta = eta + (dr / 6.) * (e1 + 2 * e2 + 2 * e3 + e4);

        if (isnan(s.states2[j + 1].bssn.chi))
        {
            cout << "ERROR: constraint solver produced nan's on step " << j  << endl;
            return; //exit(1); //should probably make a param about whether this fails
        }
    }

    //initializes alpha = sqrt(chi)-- mostly for formation runs where static lapse solver may have failed.
    for (int j = 0; j < n_gridpoints; j++)
        s.states2[j].bssn.alpha = sqrt(s.states2[j].bssn.chi);


    if (!quiet)
        cout << "Successfully ran constraint solver" << endl;

 }

//compute auxiliary quantities at jth point
void Spacetime::auxiliary_quantities_at_point(BSSNSlice* slice_ptr, int j)
{
    const double n = D - 2.;
    const double z = (j + grid_offset) * dr;

    const double invDm1 = 1.0 / (D - 1.);

    // shorthand versions of the state on desired slice for convenience
    const double chi = max(min_chi, slice_ptr->states2[j].bssn.chi);
    const double h_zz = slice_ptr->states2[j].bssn.h_zz;
    const double h_ww = slice_ptr->states2[j].bssn.h_ww;
    const double c_chris_Z = slice_ptr->states2[j].bssn.c_chris_Z;
    const double phi_re = slice_ptr->states2[j].csf.phi_re;
    const double phi_im = slice_ptr->states2[j].csf.phi_im;
    const double beta = slice_ptr->states2[j].bssn.beta;

    // inverse metric components
    h_ZZ[j] = 1.0 / h_zz;
    h_WW[j] = 1.0 / h_ww;
    const double inv_hzz = h_ZZ[j];
    const double inv_hww = h_WW[j];
    const double inv_chi = 1.0 / chi;

    //Christoffel symbols that are only needed for other aux variables
    double chris_zzz = 0.5 * d_z_h_zz;
    double chris_wwz = 0.5 * d_z_h_ww;
    double chris_zww = -0.5 * d_z_h_ww;
    if (z > min_z)
        chris_zww += (h_zz - h_ww) / z; // extra term that vanishes as z -> 0
    double chris_Wwz = inv_hww * chris_wwz;

    // christoffel symbols stored for evolution equations
    chris_Zzz[j] = inv_hzz * chris_zzz;
    chris_Zww[j] = inv_hzz * chris_zww;
    const double cZzz = chris_Zzz[j];
    const double cZww = chris_Zww[j];

    // auxiliary constraint
    if (z > min_z)
        aux_constraint[j] = (d_z_beta + n * beta / z) * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);
    else
        aux_constraint[j] = (n + 1.) * d_z_beta * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);

    // 3-covariant derivatives of alpha and tracefree parts
    D_zz_alpha[j] = d_zz(v_alpha, j) - cZzz * d_z_alpha + 0.5 * d_z_chi * d_z_alpha * inv_chi;
    
    const double d_alpha_z = (z <= min_z) ? d_zz(v_alpha, j) : (d_z_alpha / z); // regularized alpha_z/z at the origin
    D_ww_alpha[j] = d_alpha_z - (cZww + 0.5 * h_ww * inv_hzz * d_z_chi * inv_chi) * d_z_alpha;

    D_zz_alpha_TF[j] = n * (D_zz_alpha[j] - h_zz * inv_hww * D_ww_alpha[j]) * invDm1;
    D_ww_alpha_TF[j] = (D_ww_alpha[j] - h_ww * inv_hzz * D_zz_alpha[j]) * invDm1;

    // cached derivatives
    const double dzz_chi = d_zz(v_chi, j);
    const double cD_zz_chi = dzz_chi - cZzz * d_z_chi;
    double d_chi_z = (z <= min_z) ? dzz_chi : (d_z_chi / z); // d_z(chi) / z replaced with 2nd deriv at 0
    double c_chris_Z_overz = (z <= min_z) ? d_z(v_c_chris_Z, j) : (c_chris_Z / z); // c_chris_Z / z replaced with 2nd deriv at 0
    const double dzz_h_zz = d_zz(v_h_zz, j);
    const double dzz_h_ww = d_zz(v_h_ww, j);
    double d_h_ww_z = (z <= min_z) ? dzz_h_ww : (d_z_h_ww / z); // d_z(h_ww) / z replaced with 2nd deriv at 0
    const double dzchi_over_chi = d_z_chi * inv_chi;

    // Ricci chi-part
    double R_zz_chi = n * (cD_zz_chi + (0.5 * inv_hww * d_z_h_ww * d_z_chi + d_chi_z) - d_z_chi * dzchi_over_chi) * (0.5 * inv_chi);
    double R_ww_chi = h_ww * inv_hzz * (cD_zz_chi + (2. * D - 5.) * (0.5 * inv_hww * d_z_h_ww * d_z_chi + d_chi_z) - (D - 1.) * d_z_chi * (0.5 * dzchi_over_chi)) * (0.5 * inv_chi);

    // Ricci conformal-part
    double R_zz_c_t1 = (z <= min_z) ? (-0.5 * dzz_h_ww) : ((h_zz - h_ww) / (z * z) - 0.5 * d_z_h_zz / z);
    double h_zw_diff = (z <= min_z) ? (0.5 * (dzz_h_zz - dzz_h_ww)) : ((h_zz - h_ww) / (z * z));

    double R_zz_c = n * inv_hww * (R_zz_c_t1 + chris_Wwz * chris_wwz + 2. * chris_Wwz * chris_zww)
                    - 0.5 * inv_hzz * dzz_h_zz + h_zz * d_z(v_c_chris_Z, j)
                    + c_chris_Z * chris_zzz + 3. * inv_hzz * chris_zzz * cZzz;

    double R_ww_c = -0.5 * inv_hzz * dzz_h_ww + 3. * inv_hzz * chris_Wwz * chris_wwz - 0.5 * n * inv_hww * d_h_ww_z
                    + inv_hww * (cZww * chris_zww + 2. * cZww * chris_wwz)
                    + h_ww * c_chris_Z_overz + c_chris_Z * chris_wwz - inv_hww * h_zw_diff;

    R_zz[j] = R_zz_c + R_zz_chi;
    R_ww[j] = R_ww_c + R_ww_chi;

    if (use_CCZ4) //if CCZ4 is used, add supplementary terms to the (extended) Ricci tensor and fill out theta_Z
    {
        double theta = 0.5 * chi * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);
        theta_Z[j] = theta;
        R_zz[j] += theta * (chi * d_z_h_zz + h_zz * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_zzz * inv_chi;
        R_ww[j] += theta * (chi * d_z_h_ww - h_ww * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_wwz * inv_chi;
    }

    // traceless part of Ricci components
    R_zz_TF[j] = n * (R_zz[j] - h_zz * inv_hww * R_ww[j]) * invDm1;
    R_ww_TF[j] = (R_ww[j] - h_ww * inv_hzz * R_zz[j]) * invDm1;

    // |phi|
    double mod_phi = sqrt(phi_re * phi_re + phi_im * phi_im);
    if (j == 0) test_ctr = mod_phi;

    // matter quantities via ComplexScalarField helpers (reuse member to avoid per-point construction)
    CSF csf = slice_ptr->states2[j].csf;
    const BSSNState& bssn = slice_ptr->states2[j].bssn;
    // Only phi derivatives are used; K_phi derivatives and second derivatives not used in stress-energy
    CSF d_z_csf{d_z_phi_re, d_z_phi_im, 0.0, 0.0};
    CSF d_zz_csf{0.0, 0.0, 0.0, 0.0};

    rho[j] = csf_model.rho(csf, bssn, d_z_csf, d_zz_csf);
    j_z[j] = csf_model.j_z(csf, bssn, d_z_csf, d_zz_csf);
    S_zz[j] = csf_model.S_zz(csf, bssn, d_z_csf, d_zz_csf);
    S_ww[j] = csf_model.S_ww(csf, bssn, d_z_csf, d_zz_csf);

    // If configured, add massless real scalar field contributions
    if (add_real_field)
    {
        RSF rsf = slice_ptr->states2[j].rsf;
        RSF d_z_rsf{d_z(v_psi, j), d_z(v_K_psi, j)};
        RSF d_zz_rsf{0.0, 0.0};

        rho[j] += rsf_model.rho(rsf, bssn, d_z_rsf, d_zz_rsf);
        j_z[j] += rsf_model.j_z(rsf, bssn, d_z_rsf, d_zz_rsf);
        S_zz[j] += rsf_model.S_zz(rsf, bssn, d_z_rsf, d_zz_rsf);
        S_ww[j] += rsf_model.S_ww(rsf, bssn, d_z_rsf, d_zz_rsf);
    }

    // Trace S = gamma^{ij} S_ij = chi * (h^{zz} S_zz + n h^{ww} S_ww)
    S[j] = chi * (inv_hzz * S_zz[j] + n * inv_hww * S_ww[j]);

    S_zz_TF[j] = S_zz[j] - S[j] * h_zz / (chi * (D - 1.));
    S_ww_TF[j] = S_ww[j] - S[j] * h_ww / (chi * (D - 1.));
}

//compute auxiliary quantities on a slice. Should no longer be called in slice_rhs (pointwise evaluator called directly for efficiency to avoid re-computing derivatives)
void Spacetime::compute_auxiliary_quantities(BSSNSlice* slice_ptr, bool derivatives_computed)
{
    //double n = D - 2.;
    dr = R / ((double)n_gridpoints - 1.0);

    for (int j = 0; j < n_gridpoints; j++)
    {
        d_z_chi = d_z(v_chi,j);
        d_z_h_zz = d_z(v_h_zz,j);
        d_z_h_ww = d_z(v_h_ww,j);
        d_z_phi_re = d_z(v_phi_re,j);
        d_z_phi_im = d_z(v_phi_im,j);
        d_z_alpha = d_z(v_alpha,j);
        d_z_beta = d_z(v_beta,j);

        auxiliary_quantities_at_point(slice_ptr, j);
    }
}

//returns a slice corresponding to the RHS of the BSSN evolution equations. Also computes needed auxiliary quantities first
BSSNSlice Spacetime::slice_rhs(BSSNSlice* slice_ptr)
{
    double n = D - 2.;

    //define return slice and size its states array appropriately
    BSSNSlice rhs;
    rhs.R = R;
    rhs.has_BH = slice_ptr->has_BH;
    rhs.refinement_points = slice_ptr->refinement_points;

    rhs.use_CCZ4 = use_CCZ4;
    rhs.states2.resize(n_gridpoints); // composite State vector (initialized to zeros)

    int max_ref = (refinement_levels.size() == 0 ) ? 1 : refinement_levels[refinement_levels.size() - 1]; //refinement level at outermost bdry
    // Precompute 2^(max_ref-1) using bit shift (max_ref >= 1) for hot-loop conditions and stencils
    const int res_fac_outer = (max_ref > 0) ? (1 << (max_ref - 1)) : 1;

    for (int j = 0; j < n_gridpoints - 2; j++)
    {
        if (!active_points[j])//perform no calculations and skip on inactive points
            continue;

        double z = (j + grid_offset) * dr;

        //shorthand versions of the state vars on desired slice for convenience
        const double chi = max(min_chi, slice_ptr->states2[j].bssn.chi);
        const double h_zz = slice_ptr->states2[j].bssn.h_zz;
        const double h_ww = slice_ptr->states2[j].bssn.h_ww;
        const double A_zz = slice_ptr->states2[j].bssn.A_zz;
        const double A_ww = slice_ptr->states2[j].bssn.A_ww;
        const double K = slice_ptr->states2[j].bssn.K;
        const double c_chris_Z = slice_ptr->states2[j].bssn.c_chris_Z;

        const double alpha = slice_ptr->states2[j].bssn.alpha;
        const double beta = slice_ptr->states2[j].bssn.beta;


        if (wave_mode) //converts to wave eq'n solver with beta the time derivative of phi_re for testing purposes
        {
            rhs.states2[j].bssn = (BSSNState){0.,0.,0.,0.,0.,0.,0.,0.,0.};
            rhs.states2[j].csf.phi_re = beta;
            rhs.states2[j].bssn.beta = d_zz(v_phi_re, j);
            continue;
        }

        //store local variables for commonly-used derivatives to avoid unnecessary re-computation
        d_z_chi = d_z(v_chi,j);
        d_z_h_zz = d_z(v_h_zz,j);
        d_z_h_ww = d_z(v_h_ww,j);
        d_z_phi_re =  d_z(v_phi_re,j);
        d_z_phi_im =  d_z(v_phi_im,j);
        d_z_alpha =  d_z(v_alpha,j);
        d_z_beta = d_z(v_beta,j);
        d_z_K = d_z(v_K,j);

        auxiliary_quantities_at_point(slice_ptr, j);

        double beta_z = ((z <= min_z) ? d_z_beta : (beta / z)); //beta / z replaced by z-derivative at z = 0

        // precompute common inverses/constants for this point
        const double invDm1 = 1.0 / (D - 1.);
        const double inv_hzz = h_ZZ[j];
        const double inv_hww = h_WW[j];
        const double cZzz = chris_Zzz[j];
        const double cZww = chris_Zww[j];

        //fill out RHS of the BSSN evolution system
        rhs.states2[j].bssn.chi  = beta * d_z_chi  + 2.* chi * (alpha * K - d_z_beta - n * beta_z ) * invDm1;
        rhs.states2[j].bssn.h_zz = beta * d_z_h_zz + 2. * h_zz * d_z_beta - 2. * h_zz * (d_z_beta + n * beta_z) * invDm1 - 2. * alpha * A_zz;
        rhs.states2[j].bssn.h_ww = beta * d_z_h_ww - 2. * h_ww * (d_z_beta - beta_z) * invDm1 - 2. * alpha * A_ww;

        if (!use_CCZ4) //different evolution eq'n used in CCZ4
        {
            rhs.states2[j].bssn.K = beta * d_z_K - chi * inv_hzz * D_zz_alpha[j] + alpha * inv_hzz * inv_hzz * A_zz * A_zz + alpha * K * K * invDm1
                            + n * inv_hww * (alpha * A_ww * A_ww / h_ww - chi * D_ww_alpha[j]) + 8. * M_PI * alpha * (S[j] + (D - 3.) * rho[j]) / n;
        }

        rhs.states2[j].bssn.A_zz = beta * d_z(v_A_zz, j) + 2. * A_zz * d_z_beta - 2. * A_zz * (d_z_beta + n * beta_z) * invDm1 + alpha * K * A_zz
                                - 2. * alpha * A_zz * inv_hzz * A_zz + chi * (alpha * R_zz_TF[j] - D_zz_alpha_TF[j] - 8. * M_PI * alpha * S_zz_TF[j]);
        rhs.states2[j].bssn.A_ww = beta * d_z(v_A_ww, j) - 2. * A_ww * (d_z_beta -  beta_z) * invDm1 + alpha * A_ww * (K - 2 * inv_hww * A_ww)
                                + chi * (alpha * R_ww_TF[j] - D_ww_alpha_TF[j] - 8. * M_PI * alpha * S_ww_TF[j]);

        rhs.states2[j].bssn.c_chris_Z = beta * d_z(v_c_chris_Z,j) + 2. * c_chris_Z * (d_z_beta + n * beta_z) * invDm1 + inv_hzz * d_zz(v_beta,j) - c_chris_Z * d_z_beta
                                    + (D - 3.) * inv_hzz * d_zz(v_beta,j) * invDm1 - 2. * (D - 2.) * alpha * inv_hzz * d_z_K * invDm1
                                    - A_zz * inv_hzz * inv_hzz * ( (D - 1.) * alpha * d_z_chi / chi + 2 * d_z_alpha)
                                    + 2 * alpha * (cZzz * inv_hzz * inv_hzz * A_zz + n * cZww * inv_hww * inv_hww * A_ww)
                                    - sigma_BSSN * aux_constraint[j] - 16. * M_PI * alpha * j_z[j]* inv_hzz;

        if (z >= min_z) //add terms in d_z(beta) / z - beta /z^2 when z =/= 0
            rhs.states2[j].bssn.c_chris_Z += n * (inv_hww + (D - 3.) * inv_hzz * invDm1) * (d_z_beta / z - beta / (z * z));
    
        //Gauge variable update using moving puncture evolution, unless evolve_shift is off in which case do not update beta
        double f = (shock_gauge) ? (1. + shock_fac / (alpha * alpha) ): (one_log_fac/ alpha); //bona-masso f
        rhs.states2[j].bssn.alpha = beta * d_z_alpha - f * (alpha * alpha) * K;
        rhs.states2[j].bssn.beta =  (evolve_shift)? (beta * d_z_beta + gamma_fac * c_chris_Z - eta * beta): 0.;
       
        // ComplexScalarField RHS for matter variables (uses internal origin regularization and geometric couplings)
        CSF d_z_csf{d_z_phi_re, d_z_phi_im, d_z(v_K_phi_re, j), d_z(v_K_phi_im, j)};
        // Optimization: K_phi second derivatives not used by RHS; pass zeros for those entries
        CSF d_zz_csf{d_zz(v_phi_re, j), d_zz(v_phi_im, j), 0.0, 0.0};

        const CSF csf_rhs = csf_model.rhs_CSF(
            slice_ptr->states2[j].csf,
            slice_ptr->states2[j].bssn,
            d_z_csf,
            d_zz_csf,
            d_z_alpha,
            d_z_chi,
            chris_Zzz[j],
            chris_Zww[j],
            z);

        rhs.states2[j].csf = csf_rhs;

        // RealScalarField RHS for massless real scalar (optional)
        if (add_real_field)
        {
            RSF d_z_rsf{ d_z(v_psi, j), d_z(v_K_psi, j) };
            // Second derivative needed only for psi
            RSF d_zz_rsf{ d_zz(v_psi, j), 0.0 };

            const RSF rsf_rhs = rsf_model.rhs_RSF(
                slice_ptr->states2[j].rsf,
                slice_ptr->states2[j].bssn,
                d_z_rsf,
                d_zz_rsf,
                d_z_alpha,
                d_z_chi,
                chris_Zzz[j],
                chris_Zww[j],
                z);

            rhs.states2[j].rsf = rsf_rhs;
        }
        else
        {
            rhs.states2[j].rsf = RSF{0.0, 0.0};
        }

        //add supplementary terms if CCZ4 evolution is on
        if (use_CCZ4)
        {
            const double theta = slice_ptr->theta[j];
            rhs.theta.resize(n_gridpoints);

            rhs.states2[j].bssn.K = beta * d_z_K - chi * inv_hzz * D_zz_alpha[j]  - n * inv_hww * chi * D_ww_alpha[j]
                            + alpha * (chi * inv_hzz * R_zz[j] + n * chi * inv_hww * R_ww[j] + K * (K - 2. * theta)
                            - (D - 1.) * c1 * (1. + c2) * theta + 8. * M_PI * (S[j] - (D - 1.) * rho[j]) / (D - 2.));

            rhs.states2[j].bssn.A_zz += -2. * alpha * theta * A_zz;
            rhs.states2[j].bssn.A_ww += -2. * alpha * theta * A_ww;

            //note this is now \hat{\Gamma}^z, not \tilde{\Gamma}^z. So we first remove BSSN specific damping then add new terms.
            rhs.states2[j].bssn.c_chris_Z +=   sigma_BSSN * aux_constraint[j] - 2. * d_z_alpha * theta * inv_hzz + 2. * alpha * inv_hzz * d_z(v_theta, j)
                                        - 2. * c1 * alpha * theta_Z[j] / chi - 4. * alpha * K * theta_Z[j] * invDm1 / chi
                                        + sigma_BSSN * theta_Z[j]  * ((3. - D) * d_z_beta + 2. * n * beta_z) / chi;

            rhs.states2[j].bssn.alpha +=  2. * alpha * alpha * f * theta; //4. * alpha * theta; //corrects for K -> K - 2 * theta in 1+log slicing
            //rhs.states2[j].bssn.beta +=  (evolve_shift)? ( -2. * gamma_fac * theta_Z[j] / chi): 0.; //replace theta hat -> theta tilde in CCZ4 moving puncture gauge; seems to produce instability...

            if (!theta_off)
                rhs.theta[j] = beta * d_z(v_theta, j) + 0.5 * alpha * (chi * inv_hzz * R_zz[j] - inv_hzz * inv_hzz * A_zz * A_zz + n * K * K * invDm1
                            + n * (chi * inv_hww * R_ww[j] - A_ww * A_ww / (h_ww * h_ww))  - 2. * K * theta - 2. * theta_Z[j] * d_z_alpha / alpha
                            - c1 * (D + c2 * n) * theta  - 16. * M_PI * rho[j]);
            else
                rhs.theta[j] = 0.;
        }

        // add damping away from edges
        if (damping_factor != 0. && j < n_gridpoints - 3 * res_fac_outer)
            apply_damping_at_point(j, slice_ptr, rhs.states2, rhs.theta);
    }

    //apply Sommerfeld boundary conditions at outer boundary
    sommerfeld_BC(rhs, res_fac_outer, slice_ptr);
    return rhs;
}

//TODO: make these work with MR
double Spacetime::slice_mass(BSSNSlice* slice_ptr)
{
    const double n = D - 2;
    double mass = 0;

    for (int j = 0; j <  0.95 * n_gridpoints - 1; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;
    const double& chi = slice_ptr->states2[j].bssn.chi;

        if (D == 4.)
            mass += -0.25 * dr * (Ham[j]- h_ZZ[j] * chi * R_zz[j] - n * h_WW[j] * chi * R_ww[j]) * z * z * pow(chi, -1.25) ; //-1.25 is spurious! Still haven't figured out why this works...
        else
             mass += -0.125 * pow(M_PI, 0.5 * (D - 3.))  * dr * (Ham[j]- h_ZZ[j] * chi * R_zz[j] - n * h_WW[j] * chi * R_ww[j]) * pow(z, D - 2.) * pow(chi, -1.25 - 0.35 * (D - 4.)) / tgamma(0.5 * (D - 1.));
    }
    return mass;
}

double Spacetime::slice_charge(BSSNSlice* slice_ptr)
{
    const double n = D - 2;
    double charge = 0;

    for (int j = 0; j <  0.95 * n_gridpoints - 1; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;

    const double& chi = slice_ptr->states2[j].bssn.chi;
    const double& phi_re = slice_ptr->states2[j].csf.phi_re;
    const double& phi_im = slice_ptr->states2[j].csf.phi_im;
    const double& K_phi_re = slice_ptr->states2[j].csf.K_phi_re;
    const double& K_phi_im = slice_ptr->states2[j].csf.K_phi_im;

        if (D == 4.)
            charge += -8. * M_PI  * dr * (phi_re * K_phi_im - phi_im * K_phi_re) * z * z * pow(chi, -1.5);
        else
             charge += -4. * pow(M_PI, 0.5 * (D - 1.)) * dr * (phi_re * K_phi_im - phi_im * K_phi_re) * pow(z, n) * pow(chi, 0.5 * (1. - D)) / tgamma(0.5 * (D - 1.));
    }
    return charge;
}

//computes hamiltonian and momentum constraints and conformal metric determinant
void Spacetime::compute_diagnostics (BSSNSlice* slice_ptr)
{
    double n = D - 2.;
    dr = R / ((double)n_gridpoints - 1.0);

    Ham_L2 = 0.;
    Mom_L2 = 0.;

    for (int j = 0; j < n_gridpoints - 2; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;

    //shorthand versions of the BSSN vars on desired slice for convenience
    double chi = slice_ptr->states2[j].bssn.chi;
    const double h_zz = slice_ptr->states2[j].bssn.h_zz;
    const double h_ww = slice_ptr->states2[j].bssn.h_ww;
    const double A_zz = slice_ptr->states2[j].bssn.A_zz;
    const double A_ww = slice_ptr->states2[j].bssn.A_ww;
    const double K = slice_ptr->states2[j].bssn.K;
        //const double beta = slice_ptr->states2[j].bssn.beta;

        if (chi < min_chi)
            chi = min_chi;

        Ham[j] = chi * h_ZZ[j] * R_zz[j] - h_ZZ[j] * h_ZZ[j] * A_zz * A_zz + n * K * K / (D - 1.)
                 + n * (chi * h_WW[j] * R_ww[j] - A_ww * A_ww / (h_ww * h_ww)) - 16. * M_PI * rho[j];


        double A_zzww_diff = (z <= min_z) ? 0. : ((A_zz - A_ww) / z) ; //difference between A_zz and A_ww, vanishing when z = 0 as required

        Mom_Z[j] = - n * d_z(v_K,j) /(D - 1.) + n * h_WW[j] * (A_zzww_diff - - 0.5 * h_WW[j] * A_ww * d_z(v_h_ww,j))
               - n * h_WW[j] * chris_Zww[j] * A_zz + h_ZZ[j] * (d_z(v_A_zz,j) - 2. * chris_Zzz[j] * A_zz - (D - 1.) * A_zz * d_z(v_chi,j) / (2. * chi)) - 8. * M_PI * j_z[j];

        det_h[j] =  h_zz * pow(h_ww,n);

        double start_val = (make_tangherlini ? pow(0.5, 1. / (D - 3.)) : 0.);
        double end_val = R * (only_BS_violation ? (r_99 / R) : 0.9);

        //add contribution to L2 norms of Ham/ Mom constraints; z^2 factor for violation on sphere. Ad hoc cutoff radiui to avoid
        //violation being dominated by non-propagating boundary noise/ BH center spike.
       if (z > start_val && z < end_val)
        {
            Ham_L2 += dr * Ham[j] * Ham[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
            Mom_L2 += dr * Mom_Z[j] * Mom_Z[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
        }
    }

    //take square root to get L^2 norms; normalize by 16 * pi * central energy density
    Ham_L2 = sqrt(Ham_L2) / (16. * M_PI * rho0_init);
    Mom_L2 = sqrt(Mom_L2) / (16. * M_PI * rho0_init);
}

//writes diagnostic quantities to output file. maybe just make appropriate call to write_slice() instead?
void Spacetime:: write_current_slice(std::string file_name)
{
    const BSSNSlice& s = *current_slice_ptr;
    dr = R / ((double)n_gridpoints - 1.0);

    std::ofstream data_file{file_name};

     if (!data_file)
    {
        cerr << "SliceData.dat could not be opened for writing!\n";
        exit(1);
    }

    double theta0 = 0;

    for (int j = 0; j < n_gridpoints; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

        if (use_CCZ4)
            theta0 = s.theta[j];

    data_file <<  std::setprecision (16) << (j + grid_offset) * dr << "   " << s.states2[j].bssn.chi << "    " << s.states2[j].bssn.h_zz << "    " << s.states2[j].bssn.h_ww  << "    " << s.states2[j].bssn.A_zz
        << "   " << s.states2[j].bssn.A_ww << "    " << s.states2[j].bssn.K << "    " << s.states2[j].bssn.c_chris_Z  << "    " << s.states2[j].csf.phi_re << "    " << s.states2[j].csf.phi_im
        << "   " << s.states2[j].csf.K_phi_re << "    " << s.states2[j].csf.K_phi_im << "    " << s.states2[j].bssn.alpha << "    " << s.states2[j].bssn.beta << "    " << theta0 
        << "   " << s.states2[j].rsf.psi << "    " << s.states2[j].rsf.K_psi << endl;
    }
}

//writes diagnostic quantities to output file
void Spacetime:: write_diagnostics()
{
    int length = (current_slice_ptr->states2).size();

    //stick in desired auxiliary variable here to plot its value in last position
    vector<double>& aux_test = rho;

    double dr = R / ((double)length - 1.0);
    std::ofstream data_file{"Diagnostics.dat"};

     if (!data_file)
    {
        std::cerr << "Diagnostics.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 0; j < length - 2; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    const double& phi_re = current_slice_ptr->states2[j].csf.phi_re;
    const double& phi_im = current_slice_ptr->states2[j].csf.phi_im;
    double A = sqrt(phi_re * phi_re + phi_im * phi_im);

    data_file << std::setprecision (10) <<  (j + grid_offset) * dr << "   " << Ham[j] << "    " << Mom_Z[j]<< "    " << det_h[j]  << "    " << aux_test[j]
        << "    " << A << "    "  << d_zz(v_chi, j) << "    "  << d_zz(v_alpha, j)  << "    " << d_z(v_phi_re,j) << endl;
    }

    //cout << "Wrote diagnostics" << endl;
}

//read additional parameters from BSparams.par
void Spacetime::read_parameters(bool quiet)
{
    std::ifstream params{ "BSParams.par" };

    if (!params)
    {
        std::cerr << "Could not open BSParams.par\n";
        abort();
    }

    string current_line{};

    while (getline(params, current_line))
    {
        fill_parameter(current_line, "min_z = ", min_z, quiet);
        fill_parameter(current_line, "min_chi = ", min_chi, quiet);
        fill_parameter(current_line, "sigma_BSSN = ", sigma_BSSN, quiet);
        fill_parameter(current_line, "max_stored_slices = ", max_stored_slices, quiet);
        fill_parameter(current_line, "eta = ", eta, quiet);
        fill_parameter(current_line, "damping_factor = ", damping_factor, quiet);
        fill_parameter(current_line, "write_interval = ", write_interval, quiet);
        fill_parameter(current_line, "write_CN_interval = ", write_CN_interval, quiet);
        fill_parameter(current_line, "BS_resolution_factor = ", BS_resolution_factor, quiet);
        fill_parameter(current_line, "evolve_shift = ", evolve_shift, quiet);
        fill_parameter(current_line, "make_tangherlini = ", make_tangherlini, quiet);
        fill_parameter(current_line, "store_A0 = ", store_A0, quiet);
        fill_parameter(current_line, "run_quietly = ", run_quietly, quiet);
        fill_parameter(current_line, "start_time = ", start_time, quiet);
        fill_parameter(current_line, "checkpoint_time = ", checkpoint_time, quiet);
        fill_parameter(current_line, "read_thinshell = ", read_thinshell, quiet);
        fill_parameter(current_line, "cutoff_frac = ", cutoff_frac, quiet);
        fill_parameter(current_line, "only_BS_violation = ", only_BS_violation, quiet);
        fill_parameter(current_line, "run_spacetime_solver = ", run_spacetime_solver, quiet);
        fill_parameter(current_line, "gamma_fac = ", gamma_fac, quiet);
        fill_parameter(current_line, "spatially_varying_BC = ", spatially_varying_BC, quiet);
        fill_parameter(current_line, "use_CCZ4 = ", use_CCZ4, quiet);
        fill_parameter(current_line, "c1 = ", c1, quiet);
        fill_parameter(current_line, "c2 = ", c2, quiet);
        fill_parameter(current_line, "theta_off = ", theta_off, quiet);
        fill_parameter(current_line, "shock_gauge = ", shock_gauge, quiet);
        fill_parameter(current_line, "shock_fac = ", shock_fac, quiet);
        fill_parameter(current_line, "shock_gauge = ", shock_gauge, quiet);
        fill_parameter(current_line, "one_log_fac = ", one_log_fac, quiet);
        fill_parameter(current_line, "stop_on_migrate = ", stop_on_migrate, quiet);
        fill_parameter(current_line, "cell_centered = ", cell_centered, quiet);
        fill_parameter(current_line, "add_real_field = ", add_real_field, quiet);
        fill_parameter(current_line, "real_amp = ", real_amp, quiet);
        fill_parameter(current_line, "real_sigma = ", real_sigma, quiet);
        fill_parameter(current_line, "real_center = ", real_center, quiet);
        fill_parameter(current_line, "critical_study = ", critical_study, quiet);
        fill_parameter(current_line, "critical_state = ", critical_state, quiet);
        fill_parameter(current_line, "critical_eps = ", critical_eps, quiet);
        fill_parameter(current_line, "lapse_thresh = ", lapse_thresh, quiet);
        fill_parameter(current_line, "hi_guess = ", hi_guess, quiet);
        fill_parameter(current_line, "lo_guess = ", lo_guess, quiet);
        fill_parameter(current_line, "subcritical_time = ", subcritical_time, quiet);
        fill_parameter(current_line, "do_ah_search = ", do_ah_search, quiet);

        fill_param_array(current_line, "refinement_points = ", refinement_points, quiet);

    }

    cout << sigma_BSSN << eta << endl;
}

//halves number of gridpoints while keeping R. Should only be called before evolving!
void Spacetime::halve_resolution()
{
    int n_old = n_gridpoints;

    // Downsample the full composite state (geometry + matter)
    vector<State> old_states(n_old);

    //cut off one point at end to get odd # of gridpoints if needed
    if (n_gridpoints % 2 == 0)
        {
            n_gridpoints -= 1;
            R *= (n_gridpoints) / (n_old);
        }

    for (int j = 0; j < n_gridpoints; j++) // fill old array with states
        old_states[j] = slices[0].states2[j];

    n_gridpoints = (n_gridpoints + 1) / 2;

    // Write back every other point and then shrink
    for (int k = 0; k < n_gridpoints; k++)
        slices[0].states2[k] = old_states[2 * k];
    slices[0].states2.resize(n_gridpoints);

    dr = R / (n_gridpoints - 1.0);
}

void Spacetime::fill_active_points()
{
    //first deactivate all points in case of empty refinement_levels
    for (int j = 0; j < n_gridpoints; j++)
        active_points[j] = 0;

    int step_size = 1; //amount we step forward by in activating points
    unsigned int refinement_layers_crossed = 0;

    //now fill in points, halving the number of points filled in every time we pass a refinement layer
    for (int k = 0; k < n_gridpoints; k += step_size)
        {
            active_points[k] = 1;

            //when we cross refinement layer, double step size unless we're already on the last one
            if (refinement_layers_crossed < refinement_points.size() && k > refinement_points[refinement_layers_crossed]  )
            {
                refinement_layers_crossed++;
                step_size *= 2;
            }
        }
    //for (int j = 0; j < n_gridpoints; j++)
        //cout << active_points[j] << endl;
}

void Spacetime::fill_refinement_levels()
{
    refinement_levels.resize(n_gridpoints);
    for (int j = 0; j < n_gridpoints; j++)
        refinement_levels[j] = slices[0].get_refinement_level(j, refinement_points);
}
//attempts to remove the noise that accumulates at refinement boundaries
void Spacetime::kill_refinement_noise()
{
    int start_point = 0;
    for (unsigned int level = 0; level < refinement_points.size(); level++)
    {
        int step = pow(2, level);

        int j4 = start_point;
        while (j4 + step < n_gridpoints && active_points[j4 + step])
            j4 += step;

        int j3 = j4 - step;
        int j1 = j4 - 3 * step;

        // interpolate at mid-refinement boundary using cubic Lagrange polynomial
        // Nodes around j1 are at j1-3*step, j1-step, j1+step, j1+3*step
        // For uniform spacing, the weights at the center (x=0 relative to j1) are [-1/16, 9/16, 9/16, -1/16]
        BSSNSlice& s1 = *current_slice_ptr;

        const State& st1 = s1.states2[j1 - 3 * step];
        const State& st2 = s1.states2[j1 - step];
        const State& st3 = s1.states2[j1 + step];
        const State& st4 = s1.states2[j1 + 3 * step];

        s1.states2[j1] = (-1.0/16.0) * st1 + (9.0/16.0) * st2 + (9.0/16.0) * st3 + (-1.0/16.0) * st4;

        // Repeat for j3 with the same symmetric stencil around j3
        const State& tt1 = s1.states2[j3 - 3 * step];
        const State& tt2 = s1.states2[j3 - step];
        const State& tt3 = s1.states2[j3 + step];
        const State& tt4 = s1.states2[j3 + 3 * step];

        s1.states2[j3] = (-1.0/16.0) * tt1 + (9.0/16.0) * tt2 + (9.0/16.0) * tt3 + (-1.0/16.0) * tt4;

        //s1.states2[j1] = 0.5 * (s1.states2[j1 - step] + s1.states2[j1 + step]);
        //s1.states2[j3] = 0.5 * (s1.states2[j3 - step] + s1.states2[j3 + step]);

        start_point = j4;
    }
}

//corrects the initial field momentum in such a way that constraints are satisfied. Must have filled in initial slice already.
//seems to only work for negative initial perturbation, at least for early tests... now deprecated
void Spacetime::fix_initial_field_mom()
{
    compute_auxiliary_quantities(&slices[0], 0);
    compute_diagnostics(&slices[0]);
    bool failed = 0;

   // vector<double> old_momenta(n_gridpoints);

    for (int j = 0; j < n_gridpoints; j++)
    {
    //old_momenta[j] = slices[0].states2[j].csf.K_phi_im;
    double correction_sq = slices[0].states2[j].csf.K_phi_im * slices[0].states2[j].csf.K_phi_im + Ham[j] / (32. * M_PI);

        if (correction_sq < 0. && failed == 0)
            {
                failed = 1;
                cout << "WARNING: trick to adjust initial field momentum may have failed starting on gridpoint " << j << endl;
            }

       if (correction_sq >= 0.)
           slices[0].states2[j].csf.K_phi_im = - sqrt(correction_sq);
    }
}

void Spacetime::compute_dtK(int time_index)
{
    if (time_index == 0)
    {
        cout<< "WARNING: called compute_dtK on timestep 0" << endl;
        return;
    }

    BSSNSlice& slice_current = slices[time_index];
    BSSNSlice& slice_last = slices[time_index - 1];

    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j + 0.5) * dr;

        if (z < R * (only_BS_violation ? (r_99 / R) : 0.9))
        {
            double local_dtK =  (slice_current.states2[j].bssn.K - slice_last.states2[j].bssn.K ) / dt;
            dtK_L2 += dr * local_dtK * local_dtK * z * z   * pow(slice_current.states2[j].bssn.chi, -1.5);
        }
    }
    dtK_L2 = sqrt(dtK_L2);
    dtK_L2 = abs(dtK_L2 - 1);
}

void Spacetime::prepare_ham_solve()  //EXPERIMENTAL: try to return to pure isotropic coords + solve Ham constraint mid-run
{

}

//helper function to set all temp array sizes to n_gridpoints and initialize some necessary variables
void Spacetime::resize_temp_arrays()
{
    h_ZZ.resize(n_gridpoints);
    h_WW.resize(n_gridpoints);

    chris_Zww.resize(n_gridpoints);
    chris_Zzz.resize(n_gridpoints);
    aux_constraint.resize(n_gridpoints);
    D_zz_alpha.resize(n_gridpoints);
    D_ww_alpha.resize(n_gridpoints);
    D_zz_alpha_TF.resize(n_gridpoints);
    D_ww_alpha_TF.resize(n_gridpoints);
    R_zz.resize(n_gridpoints);
    R_ww.resize(n_gridpoints);
    R_zz_TF.resize(n_gridpoints);
    R_ww_TF.resize(n_gridpoints);

    rho.resize(n_gridpoints);
    j_z.resize(n_gridpoints);
    S_zz.resize(n_gridpoints);
    S_ww.resize(n_gridpoints);
    S.resize(n_gridpoints);
    S_zz_TF.resize(n_gridpoints);
    S_ww_TF.resize(n_gridpoints);

    Ham.resize(n_gridpoints);
    Mom_Z.resize(n_gridpoints);
    det_h.resize(n_gridpoints);

    if (use_CCZ4)
    {
        theta_Z.resize(n_gridpoints);
    }

    dtK_L2 = 0;

    last_active_j = n_gridpoints - 1; //find last active gridpoint
    while (!active_points[last_active_j])
        last_active_j--;
}

//add perturbation directly to spacetime field profile, thus skipping BS step
void Spacetime::add_spacetime_pert(double a, double k, double center)
{
    BSSNSlice& s = *current_slice_ptr;
    double k2 = k * k;
    double dr = R / ((double)n_gridpoints - 1.0);

    for (int j = 0; j < n_gridpoints; j++)
    {
        double phase = std::arg( std::complex<double>(s.states2[j].csf.phi_re, s.states2[j].csf.phi_im));
        double r = j * dr;

        s.states2[j].csf.phi_re += cos(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states2[j].csf.phi_im += sin(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states2[j].csf.K_phi_re += 0.5 * sin(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states2[j].bssn.alpha);
        s.states2[j].csf.K_phi_im -= 0.5 * cos(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states2[j].bssn.alpha);
    }

    solve_initial_ham(run_quietly);
}

// Add a massless real scalar Gaussian to the provided slice.
// psi(z) = real_amp * exp(-0.5 * ((z - real_center)/real_sigma)^2), K_psi(z) = 0
void Spacetime::add_real_gaussian(BSSNSlice& slice)
{
    // Guard: if disabled or degenerate sigma, do nothing
    if (!add_real_field || real_sigma <= 0.0 || real_amp == 0.0)
        return;

    const int N = n_gridpoints; // use current spacetime grid size
    const double dr_local = slice.R / ((double)N - 1.0);
    const double offset = 0.5 * (double)cell_centered; // consistent with grid_offset

    const double inv_sigma2 = 1.0 / (real_sigma * real_sigma);

    for (int j = 0; j < N; ++j)
    {
        const double z = (j + offset) * dr_local;
        const double dz = z - real_center;
        const double gauss = real_amp * exp(-0.5 * dz * dz * inv_sigma2);

        slice.states2[j].rsf.psi = gauss;
        slice.states2[j].rsf.K_psi= 0.0;
    }
}

//read data from BosonStar to spacetime and construct initial time slice
void Spacetime::initialize(BosonStar& boson_star, bool skip_read)
{
    wave_mode = 0; //make 1 for MR testing purposes only
    ah_radius = -1.;
    D = boson_star.D;

    //inherit parameters from BS
    n_gridpoints = boson_star.n_gridpoints;
    R = boson_star.R;
    courant_factor = boson_star.courant_factor;
    stop_time = boson_star.stop_time;
    mu = boson_star.mu;
    lambda = boson_star.lambda;
    sigma = boson_star.sigma;
    solitonic = boson_star.solitonic;
    omega = boson_star.omega;
    isotropic = boson_star.isotropic;
    M = boson_star.M;
    r_99 = boson_star.r_99;
    BS_perturbed = boson_star.perturb;

    if (!skip_read)
    {
        // Initialize defaults prior to reading params for backwards compatibility
        gamma_fac = 0.75; // 3/4 by default
        spatially_varying_BC = 1;
        one_log_fac = 2.0;
    }
    critical_state = 2;

    if (!skip_read)
    {   
        // Default subcritical_time if not provided
        subcritical_time = stop_time;
        refinement_points = {};
        read_parameters();
    }

    if (critical_study)
        do_ah_search = 1;
    

    if (D != 4.)
        spatially_varying_BC = 0; //only try this for D = 4 for now

    //cout << refinement_points.size() << endl;
    if (!refinement_points.empty() && refinement_points[0] <= 0) // disable refinement if sentinel present
        refinement_points.clear();

    if (read_thinshell)
        BS_resolution_factor = 1;

    if ((BS_resolution_factor & (BS_resolution_factor - 1)) != 0 || BS_resolution_factor <= 0)
        {
            cerr << "ERROR: BS_resolution_factor must be a power of 2" << endl;
            exit(1);
        }

    active_points.resize(n_gridpoints);
    fill_active_points();
    fill_refinement_levels();

    // Cache original BosonStar grid so we can restore after temporary high-res solve
    const int bs_n_gridpoints_orig = boson_star.n_gridpoints;
    bool bumped_bs_resolution = false;

    //solve BS at higher resolution and read in data to first slice, if not starting in other mode (checkpoint / thinshell read / gaussian start)
    if (BS_resolution_factor > 1 && start_time == 0 && !read_thinshell && !boson_star.gaussian_start)
    {
        boson_star.n_gridpoints = boson_star.n_gridpoints * BS_resolution_factor - BS_resolution_factor + 1;
        
        if (!run_quietly)
            cout << " \nRe-solving BS at higher resolution with " << boson_star.n_gridpoints <<  endl;

        boson_star.solve(1);
        boson_star.write_field();
        boson_star.fill_isotropic_arrays();
        boson_star.write_isotropic();

        //using higher-res omega seems like an improvement
        omega = boson_star.omega;
        bumped_bs_resolution = true;
    }

    if (read_thinshell && start_time == 0)
    {
        boson_star.read_thinshell();
        boson_star.write_field();

        if (boson_star.isotropic)
        {
            boson_star.fill_isotropic_arrays();
            boson_star.write_isotropic();
        }

        R = boson_star.R;
        n_gridpoints = boson_star.n_gridpoints;
        omega = boson_star.omega;

        cout << "R = " << boson_star.R << endl;
    }

    //add perturbation to "standard" BS
    if (!read_thinshell && !boson_star.gaussian_start && boson_star.perturb)
    {
        boson_star.add_perturbation(boson_star.perturb_amp, boson_star.perturb_spread, boson_star.perturb_center);
        boson_star.fill_given_A(omega);
    }

    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    int num_timesteps = ceil(stop_time / dt);

    slices.resize(std::min(num_timesteps + 1, max_stored_slices));
    slices[0].active_points = active_points;
    slices[0].refinement_points = refinement_points;
    slices[0].use_CCZ4 = use_CCZ4;

    //cout << "About to read" << endl;

    if (start_time == 0)
        slices[0].read_BS_data(boson_star, BS_resolution_factor, isotropic, cell_centered);

    // Restore BosonStar gridpoints so repeated initialize() calls (e.g., during critical study)
    // don't keep multiplying the resolution.
    if (bumped_bs_resolution)
        boson_star.n_gridpoints = bs_n_gridpoints_orig;

    grid_offset = 0.5 * (double)cell_centered;

    if (!run_quietly)
        cout << "Read BS data" << endl;

    // Optionally add a real scalar field Gaussian before evolution
    if (add_real_field && start_time == 0)
        add_real_gaussian(slices[0]);

    //solve Hamiltonian constraint directly in isotropic coords for initial slice if requested
    if (run_spacetime_solver && start_time == 0)
        solve_initial_ham(run_quietly);

    //cut off outermost 2 gripoints, where the christoffel symbols will be generally polluted by garbage due to not having data to take derivatives there. Might be better to just extrapolate long term.
    n_gridpoints -= 2 ;

    R *= (n_gridpoints - 1.0) / (n_gridpoints + 1.0); //also need to rescale R to avoid stretching solution
    slices[0].R = R;

    int n_old = n_gridpoints;
    n_gridpoints = round(cutoff_frac * n_gridpoints); //shrink domain by cutoff_frac, ideally to remove detritus in H
    R = (R * (n_gridpoints - 1.0)) / (n_old - 1.0); //TODO: think about this

    // Keep composite state vector in sync with the trimmed grid size
    slices[0].states2.resize(n_gridpoints);
    slices[0].R = R;

    if (use_CCZ4)
        slices[0].theta.resize(n_gridpoints);

    if (start_time > 0)
        slices[0].read_checkpoint(start_time, n_gridpoints);

    //n_gridpoints and R should not change after this point!!!

    //resize all auxiliary/diagnostic arrays as appropriate
    resize_temp_arrays();

    //compute auxiliary/diagnostic quantities on initial slice
    current_slice_ptr = &slices[0];

    if (make_tangherlini)
        M = slices[0].make_tangherlini(1., min_chi, D);

    if (make_tangherlini || slices[0].states2[0].bssn.chi < 10 * min_chi)
        slices[0].has_BH = 1;
    else slices[0].has_BH = 0;

    compute_auxiliary_quantities(current_slice_ptr);
    rho0_init = make_tangherlini ? 1. : rho[0];
    compute_diagnostics(current_slice_ptr);
    dtK_L2 = 0;

    double mass0 = slice_mass(&slices[0]);
    double charge0 = slice_charge(&slices[0]);
    M = mass0;

    std::ofstream dynamical_file{"dynamical_constants.dat"};
    dynamical_file << mass0 << "    " << charge0 << "    " << mass0 - mu * charge0;

    if (!run_quietly)
        cout << "Dynamical N =" << charge0 << ", M = "  << mass0 << " and binding energy is " << mass0 - mu * charge0  << endl;
}


//Time evolution! Must have initialized using a constructed BosonStar first.
void Spacetime::evolve()
{
    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    
    if (!run_quietly)
        cout << "dr = " << dr << ", dt = " << dt <<  "   " << endl;

    int num_timesteps = ceil(stop_time / dt);
    int last_checkpoint_time = 0;

    //write constraint norms at each timestep to file
    std::ofstream constraints_file{"constraint_norms.dat"};
    if (!constraints_file)
    {
        std::cerr << "constraint_norms.dat could not be opened for writing!\n";
        exit(1);
    }
    double A0 = sqrt (pow(slices[0].states2[0].csf.phi_re,2) + pow(slices[0].states2[0].csf.phi_im,2));

    if (!run_quietly)
        cout <<" \n Will evolve with up to " << num_timesteps << " time steps \n" << endl;
    

    //s_i are returned RHS's, t represents temporary RHS + current_slice quantities that must be stored so derivatives can be accessed
    BSSNSlice s1, s2, s3, s4, t1, t2, t3;

    for (int time_step = 0; time_step < num_timesteps; time_step++)
    {
        double t = start_time + time_step * dt;

        //fill out array until we've reached maximum number of stored slices, then update last element + rotate at end.
        int n = (time_step > max_stored_slices - 2) ? (max_stored_slices - 2) : time_step;

        double phase_ctr, phase_last;
        if (time_step > 0) phase_last = phase_ctr;
        phase_ctr = std::arg( std::complex<double>(slices[n].states2[0].csf.phi_re, slices[n].states2[0].csf.phi_im));

        double omega_approx = (phase_ctr - phase_last) / dt;

        if (omega_approx < -M_PI / dt) omega_approx += 2 * M_PI / dt; //corrects jumps due to branch cuts
        
        //double phase_diff= phase_ctr - std::fmod(omega * dt * time_step, M_PI ); //difference between actual and expected phase at center of BS
        //if (phase_diff < 0.) phase_diff += M_PI;

        if (n > 0) compute_dtK(n);
        M = slice_mass(current_slice_ptr); // maybe remove if causes bad bdry oscillations?

        double Ap0 = sqrt (pow(slices[n].states2[0].csf.phi_re,2) + pow(slices[n].states2[0].csf.phi_im,2));
        double Ap1 = sqrt (pow(slices[n].states2[1].csf.phi_re,2) + pow(slices[n].states2[1].csf.phi_im,2));
        double A_ctr =  (cell_centered) ? (Ap0- 0.5 * dr * (Ap1 - Ap0)) : Ap0; //just linearly extrapolate to r = 0 when cell-centered

        if (time_step % write_CN_interval == 0) //write time-dependent diagnostics to constraint_norms.dat
            constraints_file << std::setprecision (10) << start_time + dt * time_step << "   " << Ham_L2  
            << "   " << Mom_L2 <<  "   " << slices[n].states2[0].bssn.chi << "   "
            << A_ctr << "   "  << ah_radius << "   " <<  M << "   "  << slice_charge(current_slice_ptr)
            << "   " << dtK_L2 << "   " << omega_approx << "   " << slices[n].states2[0].bssn.alpha <<  endl;

        slices[n + 1].states2.resize(n_gridpoints);
        slices[n + 1].R = R;
        slices[n + 1].refinement_points = refinement_points;
        slices[n + 1].active_points = active_points;

        if (make_tangherlini || slices[n].states2[0].bssn.chi < 10 * min_chi)
            slices[n + 1].has_BH = 1;

        if (use_CCZ4)
            slices[n + 1].theta.resize(n_gridpoints);

        //evaluate intermediate RK4 quantities
        current_slice_ptr = &slices[n];

        s1 = slice_rhs(current_slice_ptr);
        t1 = slices[n]  + (0.5 * dt) * s1;

        compute_diagnostics(current_slice_ptr); //do this here so current_slice_ptr is in right place and auxiliary quantities computed

        current_slice_ptr = &t1; //must update current_slice_ptr before calling slice_rhs or derivatives will not work properly! Should consider better approach...
        s2 = slice_rhs(current_slice_ptr);
        t2 = slices[n] + (0.5 * dt) * s2;

        current_slice_ptr = &t2;
        s3 = slice_rhs(current_slice_ptr);
        t3 = slices[n] + dt * s3;

        current_slice_ptr = &t3;
        s4 = slice_rhs(current_slice_ptr);

        //update slice
        slices[n + 1] = slices[n] + (dt / 6.) * (s1 + 2. * s2 + 2. * s3 + s4);

        //enforce that A is traceless
        current_slice_ptr = &slices[n + 1];
        kill_refinement_noise(); //interpolate over refinement boundary noise before making A traceless
        make_A_traceless(current_slice_ptr);

        //enforce minimum chi
        for (State& st: slices[n + 1].states2)
            {if (st.bssn.chi < min_chi) st.bssn.chi = min_chi; }
        
        if (do_ah_search)
        {   
            double ah_radius_old = ah_radius;
            ah_radius = apparent_horizon_search(current_slice_ptr);
            if ( ah_radius_old < 0. && ah_radius > 0. && !critical_study)
                cout << "Apparent horizon found at radius " << ah_radius << " at time " << t + dt << endl;
        }

        if (time_step % write_interval == 0)
        {
            //current_slice_ptr->write_slice();
            write_current_slice();
            write_diagnostics();
        }

        //write checkpoint files
        if ((int)std::floor(t) % checkpoint_time == 0 && (int)std::floor(t) > last_checkpoint_time)
        {
            current_slice_ptr->write_slice("checkpoint" + std::to_string((int)std::floor(t)) + ".dat");
            last_checkpoint_time = (int)std::floor(t);
            cout << "Wrote checkpoint at t = " << last_checkpoint_time << endl;
        }

        //cycles slice array back by one so that last entry can be overwritten
        if (time_step >= max_stored_slices - 2)
            rotate(slices.begin(), slices.begin() + 1, slices.end());

        if ((time_step + 1) % 10 == 0 && !run_quietly) cout << "Time step " << time_step + 1 << " complete! t = " << t << endl;

        if (isnan(slices[n + 1].states2[0].bssn.chi))
        {
            cout << "Central chi became nan on step " << time_step << endl;
            exit(1);
        }

        //if stop_on_migrate enabled, exit when central amplitude changes by more than 10%
        if (stop_on_migrate && abs(A0 - test_ctr) > 0.1 * A0)
        {
            cout << "Migration occurred on time step  " << time_step <<" at time t = " << time_step * dt << endl;
            exit(1);
        }

        if (critical_study && ah_radius > 0.) //use AH formation as supercritical indicator 
        //(critical_study && slices[n + 1].states2[0].bssn.alpha < lapse_thresh)
        {
            cout << "Supercritical at time " << t << endl;
            critical_state = 1;
            critical_time = t;
            break;
        }
        //wait some time to avoid initial transients, then use net central amplitude decrease as subcritical indicator
        if (critical_study && A_ctr < A0 && t > 10.0) 
        {
            cout << "Subcritical at time " << t << endl;
            critical_state = 0;
            critical_time = t;
            break;
        }

        // Additional subcritical criterion: if no AH by subcritical_time, assume subcritical
        if (critical_study && t >= subcritical_time)
        {
            if (ah_radius < 0.0)
            {
                cout << "Subcritical by time criterion at time " << t << endl;
                critical_state = 0;
                critical_time = t;
                break;
            }
        }
    }
}

// Tune a parameter by bisection to bracket the critical threshold.
// hi_guess must produce supercritical (critical_state = 1), lo_guess subcritical (critical_state = 0).
// Will repeatedly initialize and evolve, restoring tuning_param after initialize in (likely) case it's a Spacetime member overwritten by params.
void Spacetime::tune_to_critical(double& tuning_param, double hi_guess, double lo_guess, BosonStar* bs)
{
    std::ofstream subcritical_file{"subcritical.dat"};
    std::ofstream supercritical_file{"supercritical.dat"};

    subcritical_file << "#tuning_param   critical_time" << endl;
    supercritical_file << "#tuning_param   critical_time" << endl;

    if (bs == nullptr)
    {
        cerr << "ERROR: tune_to_critical called with null BosonStar pointer" << endl;
        exit(1);
    }

    // Ensure we record critical outcomes during evolve()
    critical_study = true;

    auto run_with_param = [&](double val) -> int {
        // Assign target value, initialize (which may overwrite), then restore
        tuning_param = val;
        double tmp = tuning_param;
        cout << "Running with tuning_param = " << std::setprecision(16) << tuning_param << endl;

        
        // Skip re-reading parameters so tuned value persists into initial data seeding
        initialize(*bs, /*skip_read=*/true);
        tuning_param = tmp;
        evolve();
        if (critical_state != 0 && critical_state != 1)
        {
            cerr << "ERROR: evolve() ended with critical_state = " << critical_state
                 << ", expected 0 (subcritical) or 1 (supercritical)." << endl;
            exit(1);
        }
        return critical_state;
    };

    // Verify high guess is supercritical
    int hi_state = run_with_param(hi_guess);
    if (hi_state != 1)
    {
        cerr << "ERROR: hi_guess = " << std::setprecision(16) << hi_guess
             << " did not produce supercritical outcome (critical_state = 1). Got critical_state = "
             << hi_state << "." << endl;
        exit(1);
    }
    supercritical_file << std::setprecision(16)
                               << hi_guess<< "   " << critical_time << endl;

    // Verify low guess is subcritical
    int lo_state = run_with_param(lo_guess);
    if (lo_state != 0)
    {
        cerr << "ERROR: lo_guess = " << std::setprecision(16) << lo_guess
             << " did not produce subcritical outcome (critical_state = 0). Got critical_state = "
             << lo_state << "." << endl;
        exit(1);
    }

    subcritical_file << std::setprecision(16)
                             << lo_guess << "   " << critical_time << endl;
    double hi = hi_guess;
    double lo = lo_guess;

    // Bisection until interval smaller than tolerance
    while (fabs(hi - lo) / lo > critical_eps)
    {
        double mid = 0.5 * (hi + lo);
        int state = run_with_param(mid);
        if (state == 1)
        {
            hi = mid;   // supercritical at mid -> reduce upper bound
            supercritical_file << std::setprecision(16)
                               << mid << "   " << critical_time << endl;
        }
        else if (state == 0)
        {
            lo = mid;   // subcritical at mid -> raise lower bound
            subcritical_file << std::setprecision(16)
                             << mid << "   " << critical_time << endl;
        }
        else
        {
            cerr << "ERROR: Unexpected critical_state during bisection: " << state << endl;
            exit(1);
        }
    }

    tuning_param = 0.5 * (hi + lo);
    cout << std::setprecision(16)
         << "Critical tuning complete. tuning_param = " << tuning_param
         << ", bracket [" << lo << ", " << hi << "] with tol = " << critical_eps << endl;
    exit(0);
}

// --------------------------
// Apparent horizon search
// --------------------------

// Placeholder AH indicator at gridpoint j. This will be replaced with the actual expression
// using BSSNState fields (e.g., expansion of outgoing null geodesics).
double Spacetime::ah_indicator_at(const BSSNSlice* slice_ptr, int j) 
{
    // Access BSSN variables to avoid unused-variable warnings and as a cue for future implementation.
    const BSSNState& s = slice_ptr->states2[j].bssn;
    double z = (j + grid_offset) * dr;
    if (z == 0.0)
        z = 1e-12; // avoid division by zero at origin

    double ah_indicator = sqrt(s.chi / s.h_zz) *
                          (2. / z + (d_z(v_h_ww, j) - s.h_ww * d_z(v_chi, j) / s.chi) / s.h_ww)
                          - (D - 2.) * (s.A_ww + s.K * s.h_ww / (D - 1.));

    return ah_indicator;
}

// Scan slice for sign changes of the AH indicator and return outermost root radius. Returns -1.0 if no AH is found.
double Spacetime::apparent_horizon_search(const BSSNSlice* slice_ptr)
{
    if (slice_ptr == nullptr || slice_ptr->states2.empty())
        return -1.0;

    const int N = static_cast<int>(slice_ptr->states2.size());
    const double dr_local = slice_ptr->R / (static_cast<double>(N) - 1.0);
    const double offset = 0.5 * static_cast<double>(cell_centered); // match grid_offset convention

    double outermost_r = -1.0;

    double f_prev = ah_indicator_at(slice_ptr, 0);
    double r_prev = (0 + offset) * dr_local;

    for (int j = 1; j < N; ++j)
    {
        double f_curr = ah_indicator_at(slice_ptr, j);
        double r_curr = (j + offset) * dr_local;

        // Exact zero at a node counts as a root
        if (f_curr == 0.0)
            outermost_r = r_curr;

        // Look for sign changes indicating a crossing
        if ((f_prev > 0.0 && f_curr < 0.0) || (f_prev < 0.0 && f_curr > 0.0))
        {
            // Linear interpolation to estimate root location within [r_prev, r_curr]
            double r_zero = r_prev - f_prev * (r_curr - r_prev) / (f_curr - f_prev);
            if (r_zero > outermost_r)
                outermost_r = r_zero;
        }

        f_prev = f_curr;
        r_prev = r_curr;
    }
    return outermost_r;
}

#endif /*EVOLUTIONVARIABLES*/