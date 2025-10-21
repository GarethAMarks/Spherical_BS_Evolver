#ifndef EVOLUTIONVARIABLES_CPP_
#define EVOLUTIONVARIABLES_CPP_

#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include "EvolutionVariables.h"
#include "mathutils.h"
#include <iomanip>
#include<algorithm>
#include<filesystem>
#include <complex.h>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::rotate;


//overloads for addition/ scalar multiplication of BSSNState sets and the slice arrays containing them
BSSNState operator+(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi + s2.chi, s1.h_zz + s2.h_zz, s1.h_ww + s2.h_ww, s1.A_zz + s2.A_zz, s1.A_ww + s2.A_ww, s1.K + s2.K, s1.c_chris_Z + s2.c_chris_Z,

    s1.phi_re + s2.phi_re, s1.phi_im + s2.phi_im, s1.K_phi_re + s2.K_phi_re, s1.K_phi_im + s2.K_phi_im, s1.alpha + s2.alpha, s1.beta + s2.beta};
}

BSSNState operator-(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi - s2.chi, s1.h_zz - s2.h_zz, s1.h_ww - s2.h_ww, s1.A_zz - s2.A_zz, s1.A_ww - s2.A_ww, s1.K - s2.K, s1.c_chris_Z - s2.c_chris_Z,

    s1.phi_re - s2.phi_re, s1.phi_im - s2.phi_im, s1.K_phi_re - s2.K_phi_re, s1.K_phi_im - s2.K_phi_im, s1.alpha - s2.alpha, s1.beta - s2.beta};
}

//termwise multiplication for convenient use with characteristic speeds etc.
BSSNState operator*(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi * s2.chi, s1.h_zz * s2.h_zz, s1.h_ww * s2.h_ww, s1.A_zz * s2.A_zz, s1.A_ww * s2.A_ww, s1.K * s2.K, s1.c_chris_Z * s2.c_chris_Z,

    s1.phi_re * s2.phi_re, s1.phi_im * s2.phi_im, s1.K_phi_re * s2.K_phi_re, s1.K_phi_im * s2.K_phi_im, s1.alpha * s2.alpha, s1.beta * s2.beta};
}

BSSNState operator*(double c, const BSSNState& s)
{
    return (BSSNState){c * s.chi, c * s.h_zz, c * s.h_ww , c * s.A_zz, c * s.A_ww, c * s.K, c * s.c_chris_Z,

    c * s.phi_re, c * s.phi_im, c * s.K_phi_re, c * s.K_phi_im, c * s.alpha, c * s.beta};
}

inline BSSNState operator/(const BSSNState& s, double c)
{
    return (1. / c) * s;
}

BSSNSlice operator+(const BSSNSlice& slice1, const BSSNSlice& slice2)
{
    BSSNSlice return_slice;
    unsigned int length = slice1.states.size();

    if (length != slice2.states.size() || slice1.R != slice2.R)
    {
        cerr << "ERROR: attempted to add data from slices of different sizes!" << endl;
        exit(1);
    }

    return_slice.states.resize(length);

    if (slice1.use_CCZ4)
        return_slice.theta.resize(length);


    return_slice.R = slice1.R;
    return_slice.has_BH = slice1.has_BH;
    return_slice.use_CCZ4 = slice1.use_CCZ4;
    return_slice.refinement_points = slice1.refinement_points;

    for (unsigned int j = 0; j < length; j++)
    {
        return_slice.states[j] = slice1.states[j] + slice2.states[j];

        if (slice1.use_CCZ4)
            return_slice.theta[j] = slice1.theta[j] + slice2.theta[j];
    }

    return return_slice;
}

BSSNSlice operator*(double c, const BSSNSlice& slice)
{
    BSSNSlice return_slice;

    int length = slice.states.size();
    return_slice.states.resize(length);

    if (slice.use_CCZ4)
        return_slice.theta.resize(length);

    return_slice.R = slice.R;
    return_slice.has_BH = slice.has_BH;
    return_slice.refinement_points = slice.refinement_points;
    return_slice.use_CCZ4 = slice.use_CCZ4;

    for (int j = 0; j < length; j++)
    {
        return_slice.states[j] = c * slice.states[j];

        if (slice.use_CCZ4)
            return_slice.theta[j] = c * slice.theta[j];
    }

    return return_slice;
}


//radial partial derivative of a given bssn var at index. Order is an optional argument that allows higher z-derivatives to be taken, default is 1.
//should add chacks on refinement levels
double BSSNSlice::d_z(bssn_var var, int index, int order = 1 )
{
    int n_gridpoints = states.size();
    // Use cell-centred grid spacing to match read_BS_data (innermost point at z = 0.5*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    int ref_level = get_refinement_level(index, refinement_points); //refinement level; starts at 1 for no refinement and halves every time.

    int res_fac = pow(2, ref_level - 1);

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
                    F[j] = states[J[j]].chi;
                    break;
                case v_h_zz:
                    F[j] = states[J[j]].h_zz;
                    break;
                case v_h_ww:
                    F[j] = states[J[j]].h_ww;
                    break;
                case v_A_zz:
                    F[j] = states[J[j]].A_zz;
                    break;
                case v_A_ww:
                    F[j] = states[J[j]].A_ww;
                    break;
                case v_K:
                    F[j] = states[J[j]].K;
                    break;
                case v_c_chris_Z:
                    F[j] = states[J[j]].c_chris_Z;
                    break;
                case v_phi_re:
                    F[j] = states[J[j]].phi_re;
                    break;
                case v_phi_im:
                    F[j] = states[J[j]].phi_im;
                    break;
                case v_K_phi_re:
                    F[j] = states[J[j]].K_phi_re;
                    break;
                case v_K_phi_im:
                    F[j] = states[J[j]].K_phi_im;
                    break;
                case v_alpha:
                    F[j] = states[J[j]].alpha;
                    break;
                case v_beta:
                    F[j] = states[J[j]].beta;
                    break;
                case v_theta:
                    F[j] = theta[J[j]];
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
    int n_gridpoints = states.size();

    // Use cell-centred grid spacing to match slice convention (z = (j+0.5)*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    double psi_power = -4. / (D - 3.); //12. / ((D - 3.) * (1. - D));

    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = ( (j + 0.5) * dr);

        //conformal factor in isotropic convention
        double psi = 1. + pow(m, D - 3.) / (4 * pow(r, D - 3.));

        //conformal factor and conformally rescaled metric components are all powers of X in this gauge
        states[j].chi = max(min_chi, pow(psi, psi_power));
        states[j].h_zz = 1.;
        states[j].h_ww = 1.;

        states[j].A_zz = 0.;
        states[j].A_ww = 0.;
        states[j].K = 0.;
        states[j].c_chris_Z = 0.;
        states[j].phi_re = 0.;
        states[j].phi_im = 0.;
        states[j].K_phi_re = 0.;
        states[j].K_phi_im = 0.;

        states[j].alpha = pow(states[j].chi, 0.5) /*abs((1 - 0.5 * m / r) / (1 + 0.5 * m / r))*/ /*sqrt(abs((4 * pow(r, D - 3.) - pow(m, D - 3.)) / (4 * pow(r, D - 3.) + pow(m, D - 3.))))*/;
        states[j].beta = 0.;

    }
    //states[0].alpha = 1e-4;
    //cout << states[0].alpha << endl;
    return m;
}

//reads data from a boson_star object we have already solved for into initial BSSN slice.
void BSSNSlice::read_BS_data (BosonStar& boson_star, int BS_resolution_factor, bool isotropic, bool cell_centered)
{
    //prepare to fill states array

    if ((boson_star.n_gridpoints + BS_resolution_factor - 1) % BS_resolution_factor != 0 )
        cout << "WARNING: Incompatible resolution factor used in read_BS_data" << endl;

    int n_gridpoints = (boson_star.n_gridpoints + BS_resolution_factor - 1) / BS_resolution_factor;
    states.resize(n_gridpoints);

    if (use_CCZ4)
       theta.resize(n_gridpoints);

    R = boson_star.R;
    // Use a cell-centred z-grid so the innermost gridpoint is at z = 0.5 * dr
    // (previously grid was node-centred with z = j * dr and dr = R/(N-1)).
    double dr = R / ((double)n_gridpoints - 1.0);

    double D = boson_star.D;

    //
    // Precompute BosonStar spacing and sizes for interpolation
    int N_bs = boson_star.n_gridpoints;
    double dr_bs = boson_star.R / (double)(max(1, N_bs - 1));

    for (int j = 0; j < n_gridpoints; j++)
    {
        // cell- or vertex-centered coordinate
        double z = (j + 0.5 * (double)cell_centered) * dr;

        // find nearest lower node index on BS grid and set j0 for cubic_interp
        int idx = static_cast<int>(floor(z / dr_bs));
        int j0 = idx - 1; // cubic_interp expects r between (j0+1)*dr and (j0+2)*dr
        j0 = bound(j0, 0, max(0, N_bs - 4));

        // Interpolate X (areal conformal factor) or psi (isotropic conformal factor)
        double X_interp = 1.;
        double psi_interp = 1.;
        if (isotropic)
        {
            // interpolate isotropic psi (conformal factor) and A
            psi_interp = cubic_interp(z, boson_star.psi_iso_array[j0], boson_star.psi_iso_array[j0 + 1], boson_star.psi_iso_array[j0 + 2], boson_star.psi_iso_array[j0 + 3], j0, dr_bs);
            states[j].chi = pow(psi_interp, -4.);
            states[j].h_zz = 1.; //isotropic gauge is conformally flat
            states[j].h_ww = 1.;
        }
        else
        {
            X_interp = cubic_interp(z, boson_star.state[j0].X, boson_star.state[j0 + 1].X, boson_star.state[j0 + 2].X, boson_star.state[j0 + 3].X, j0, dr_bs);
            states[j].chi = pow(X_interp, -2. / (D - 1.));
            states[j].h_zz = pow(X_interp, (2. * D - 4.) / (D - 1.));
            states[j].h_ww = pow(X_interp, -2. / (D - 1.));
        }

        //for time-independent static, spherically symmetric boson stars K_ij = 0-- may need to do more work when we consider perturbations.
        states[j].A_zz = 0.;
        states[j].A_ww = 0.;
        states[j].K = 0.;

        // start at t = 0 so the scalar field is real regardless of phase
        double A_interp = 0.;
        double phi_interp = 0.;
        if (isotropic)
        {
            A_interp = cubic_interp(z, boson_star.A_iso_array[j0], boson_star.A_iso_array[j0 + 1], boson_star.A_iso_array[j0 + 2], boson_star.A_iso_array[j0 + 3], j0, dr_bs);
            phi_interp = cubic_interp(z, boson_star.phi_iso_array[j0], boson_star.phi_iso_array[j0 + 1], boson_star.phi_iso_array[j0 + 2], boson_star.phi_iso_array[j0 + 3], j0, dr_bs);
            states[j].phi_re = A_interp;
            states[j].phi_im = 0.;
            states[j].alpha = exp(phi_interp);
            states[j].beta = 0.;
        }
        else
        {
            A_interp = cubic_interp(z, boson_star.state[j0].A, boson_star.state[j0 + 1].A, boson_star.state[j0 + 2].A, boson_star.state[j0 + 3].A, j0, dr_bs);
            phi_interp = cubic_interp(z, boson_star.state[j0].phi, boson_star.state[j0 + 1].phi, boson_star.state[j0 + 2].phi, boson_star.state[j0 + 3].phi, j0, dr_bs);
            states[j].phi_re = A_interp;
            states[j].phi_im = 0.;
            states[j].alpha = exp(phi_interp);
            states[j].beta = 0.;
        }

        //starting BS real means its momentum is imaginary (with 0 starting shift)
        states[j].K_phi_re = 0.;

        double pert_correction = 0.;
        if (boson_star.perturb)
        {
            double pert_interp = 0.;
            if (isotropic)
                pert_interp = cubic_interp(z, boson_star.pert_iso_array[j0], boson_star.pert_iso_array[j0 + 1], boson_star.pert_iso_array[j0 + 2], boson_star.pert_iso_array[j0 + 3], j0, dr_bs);
            else
                pert_interp = cubic_interp(z, boson_star.pert_array[j0], boson_star.pert_array[j0 + 1], boson_star.pert_array[j0 + 2], boson_star.pert_array[j0 + 3], j0, dr_bs);

            pert_correction = 1. * (boson_star.omega / states[j].alpha) * pert_interp; //conserves Noether charge at 1st order in perturbation
        }

        states[j].K_phi_im = - boson_star.omega * states[j].phi_re / (2. * states[j].alpha) + pert_correction;

        if (use_CCZ4)
            theta[j] = 0.;
    }

    //separate loop to fill in values for the contracted Christoffel symbols as they require derivatives of h (hence info off j value)
    // separate loop to fill in values for the contracted Christoffel symbols as they require
    // derivatives of h (hence info off j value). With a cell-centred grid z = (j+0.5)*dr
    // so z>0 for all j and the extra geometric term can be evaluated safely for every j.
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j) * dr;

        // inverse metric components, for convenience
        double h_ZZ = 1 / states[j].h_zz;
        double h_WW = 1 / states[j].h_ww;

        states[j].c_chris_Z = 0.5 * h_ZZ * h_ZZ * d_z(v_h_zz,j) + (D - 2.) * (-0.5 * h_ZZ * h_WW * d_z(v_h_ww,j));

        // extra geometric term, which goes to zero at z = 0
        if (z > 1e-8)
            states[j].c_chris_Z += (D - 2.) * (h_WW * (1 - h_ZZ * states[j].h_ww) / z);

        // hard enforce 0 for isotropic data
        if (isotropic)
            states[j].c_chris_Z = 0;
    }
}

//writes BSSN evolution variables from a particular slice to text file
void BSSNSlice::write_slice(std::string file_name)
{
    int length = states.size();
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


    data_file <<  std::setprecision (16) << (j) * dr << "   " << states[j].chi << "    " << states[j].h_zz << "    " << states[j].h_ww  << "    " << states[j].A_zz
        << "   " << states[j].A_ww << "    " << states[j].K << "    " << states[j].c_chris_Z  << "    " << states[j].phi_re << "    " << states[j].phi_im
        << "   " << states[j].K_phi_re << "    " << states[j].K_phi_im << "    " << states[j].alpha << "    " << states[j].beta << "    " << theta0 << endl;
    }

}

//reads slice data from checkpoint file of form checkpoint(time).dat w/o parentheses
void BSSNSlice::read_checkpoint(int time, int n_gridpoints)
{
    states.resize(n_gridpoints);

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
        BSSNState& s = states[j];

        if (use_CCZ4)
        {
            double& t = theta[j];

            if(iss >> z >> s.chi >> s.h_zz >> s.h_ww >> s.A_zz >> s.A_ww >> s.K >> s.c_chris_Z >> s.phi_re >> s.phi_im >> s.K_phi_re >> s.K_phi_im >> s.alpha >> s.beta >> t)
                j++;
        }
        else
        {
            if(iss >> z >> s.chi >> s.h_zz >> s.h_ww >> s.A_zz >> s.A_ww >> s.K >> s.c_chris_Z >> s.phi_re >> s.phi_im >> s.K_phi_re >> s.K_phi_im >> s.alpha >> s.beta)
                j++;
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

    if (2*sl.states.size() - 1 != sm.states.size()  || 2*sm.states.size() - 1 != sh.states.size()   )
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
    for (unsigned int j = 0; j < sl.states.size(); j++)
    {
    // For 4th-order convergence the refinement factor scaling is 2^4 = 16
    conv_file << j * sl.R /(sl.states.size() - 1) << "   " << sm.states[2*j].phi_re -  sl.states[j].phi_re  << "    " << 16.*(sh.states[4*j].phi_re  -  sm.states[2*j].phi_re ) << endl;
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
    /*if (!make_tangherlini)
        return 1;
    else*/
        return sqrt(1 - 2. * M / R);

}

double Spacetime::d_chi_asymp(double r)
{
    return(isotropic) ? (2 * M * pow(1. +  0.5 * M / r, -5.) / (r * r)) : -(2./3.) * M * pow(1 - 2. * M / r, -4./3.) / (r * r);
}

double Spacetime::d_alpha_asymp(double r)
{
    /*if (!make_tangherlini)
        return 0;
    else*/
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
    double phi_re = cubic_interp(r, s.states[j0].phi_re, s.states[j0 + 1].phi_re, s.states[j0 + 2].phi_re, s.states[j0 + 3].phi_re, j0, dr);
    double phi_im = cubic_interp(r, s.states[j0].phi_im, s.states[j0 + 1].phi_im, s.states[j0 + 2].phi_im,s.states[j0 + 3].phi_im, j0, dr);
    double K_phi_re = cubic_interp(r, s.states[j0].K_phi_re, s.states[j0 + 1].K_phi_re, s.states[j0 + 2].K_phi_re, s.states[j0 + 3].K_phi_re,j0, dr);
    double K_phi_im = cubic_interp(r, s.states[j0].K_phi_im, s.states[j0 + 1].K_phi_im, s.states[j0 + 2].K_phi_im, s.states[j0 + 3].K_phi_im,j0, dr);

    //WARNING: these are very similarly named to Spacetime member variables introduced to avoid re-computing derivatives during evolution; do not confuse!
    double dz_phi_re = cubic_interp(r, d_z(v_phi_re, j0), d_z(v_phi_re, j0 + 1), d_z(v_phi_re, j0 + 2), d_z(v_phi_re, j0 + 3),j0,  dr);
    double dz_phi_im = cubic_interp(r, d_z(v_phi_im, j0), d_z(v_phi_im, j0 + 1), d_z(v_phi_im, j0 + 2), d_z(v_phi_im, j0 + 3),j0,  dr);

   /* if (j0 == 0)
    {
        phi_re = s.states[0].phi_re * (dr - r) / dr + s.states[1].phi_re * r / dr;
        phi_im = s.states[0].phi_im * (dr - r) / dr + s.states[1].phi_im * r / dr;

        K_phi_re = s.states[0].K_phi_re * (dr - r) / dr + s.states[1].K_phi_re * r / dr;
        K_phi_im = s.states[0].K_phi_im * (dr - r) / dr + s.states[1].K_phi_im * r / dr;

        dz_phi_re = d_z(v_phi_re, 0) * (dr - r) / dr + d_z(v_phi_re, 1) * r / dr;
        dz_phi_im = d_z(v_phi_im, 0) * (dr - r) / dr + d_z(v_phi_im, 1) * r / dr;
    }*/

    double mod_phi = sqrt(phi_re * phi_re + phi_im + phi_im);
    double rho0 = 2. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) +  0.5 * chi * (pow(dz_phi_re,2) + pow(dz_phi_im,2)) + 0.5 * V(mod_phi);

    double chi_rhs = eta;
    double eta_rhs = 5. * eta * eta / (4. * chi) + 8. * M_PI * rho0;

    if (r > 0.) eta_rhs -= 2 * eta / r;
    else eta_rhs /= 3.; // in the limit r -> 0, eta / r is replaced by d_r(eta), so net effect is to divide the usual RHS by 3!

    vector<double> rhs = {chi_rhs, eta_rhs};
    return rhs;

 }

 void Spacetime::solve_initial_ham()
 {
    if (!isotropic)
    {
        cout << "WARNING: Spacetime Hamiltonian solver called in non-isotropic coordinates; aborting request" << endl;
        return;
    }

    BSSNSlice& s = slices[0]; //initial slice
    double c1, c2, c3, c4, e1, e2, e3, e4;
    double eta = 0;

    cout << "\n Starting initial Spacetime Hamiltonian constraint solver..." << endl;

   for (int j = 0; j < n_gridpoints - 1; j++) //for (int j = n_gridpoints - 1; j > 0; j--)
    {
        double r = j * dr;

        //RK4 solver for Ham. constraint
        //pretty inefficient but OK for now
        c1 = ham_init_rhs(r, s.states[j].chi, eta)[0];
        e1 = ham_init_rhs(r, s.states[j].chi, eta)[1];

        c2 = ham_init_rhs(r + dr / 2., s.states[j].chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[0];
        e2 = ham_init_rhs(r + dr / 2., s.states[j].chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[1];

        c3 = ham_init_rhs(r + dr / 2., s.states[j].chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[0];
        e3 = ham_init_rhs(r + dr / 2., s.states[j].chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[1];

        c4 = ham_init_rhs(r + dr, s.states[j].chi + dr * c3, eta + dr * e3)[0];
        e4 = ham_init_rhs(r + dr, s.states[j].chi + dr * c3, eta + dr * e3)[1];

        s.states[j + 1].chi = s.states[j].chi + (dr / 6.) * (c1 + 2 * c2 + 2 * c3 + c4);
        eta = eta + (dr / 6.) * (e1 + 2 * e2 + 2 * e3 + e4);

        if (isnan(s.states[j + 1].chi))
        {
            cout << "ERROR: constraint solver produced nan's on step " << j  << endl;
            return; //exit(1); //should probably make a param about whether this fails
        }
    }

    //initializes alpha = sqrt(chi)-- mostly for formation runs where static lapse solver may have failed.
    for (int j = 0; j < n_gridpoints; j++)
        s.states[j].alpha = sqrt(s.states[j].chi);


    cout << "Successfully ran constraint solver" << endl;

 }

//compute auxiliary quantities at jth point
void Spacetime::auxiliary_quantities_at_point(BSSNSlice* slice_ptr, int j)
{
    double n = D - 2.;
    double z = (j + grid_offset) * dr;

    //shorthand versions of the BSSN vars on desired slice for convenience
    const double chi = max(min_chi, slice_ptr->states[j].chi);
    const double h_zz = slice_ptr->states[j].h_zz;
    const double h_ww = slice_ptr->states[j].h_ww;
    const double c_chris_Z = slice_ptr->states[j].c_chris_Z;
    const double phi_re = slice_ptr->states[j].phi_re;
    const double phi_im = slice_ptr->states[j].phi_im;
    const double K_phi_re = slice_ptr->states[j].K_phi_re;
    const double K_phi_im = slice_ptr->states[j].K_phi_im;
    const double beta = slice_ptr->states[j].beta;

    //inverse metric components
    h_WW[j] = 1 / h_ww;
    h_ZZ[j] = 1 / h_zz;

    //Christoffel symbols that are only needed for other aux variables
    double chris_zzz = 0.5 * d_z_h_zz;
    double chris_wwz = 0.5 * d_z_h_ww;
    double chris_zww = -0.5 * d_z_h_ww;
    if (z > min_z)
        chris_zww +=  (h_zz - h_ww) / z; //extra term that vanishes as z -> 0
    double chris_Wwz = h_WW[j] * chris_wwz;

    //christoffel symbols that are stored for evolution equations
    chris_Zww[j] = h_ZZ[j] * chris_zww;
    chris_Zzz[j] = h_ZZ[j] * chris_zzz;

    //auxiliary constraint, accounting for limit beta / z -> d_z(beta) as z -> 0
    if (z > min_z)
        aux_constraint[j] = (d_z_beta + n * beta / z) * (c_chris_Z - h_ZZ[j] * chris_Zzz[j] - n * h_WW[j] * chris_Zww[j] );
    else
        aux_constraint[j] = (n + 1.) * d_z_beta  * (c_chris_Z - h_ZZ[j] * chris_Zzz[j] - n * h_WW[j] * chris_Zww[j] );

    //3-covariant derivatives of alpha and tracefree parts
    double d_alpha_z = ((z <= min_z) ? d_zz(v_alpha,j) : (d_z_alpha / z)); //d_z(alpha) / z replaced with 2nd deriv at 0

    D_zz_alpha[j] = d_zz(v_alpha, j) - chris_Zzz[j] * d_z_alpha + d_z_chi * d_z_alpha / (2. * chi);
    D_ww_alpha[j] = d_alpha_z - (chris_Zww[j] + 0.5 * h_ww * h_ZZ[j] * d_z_chi / chi ) * d_z_alpha;

    D_zz_alpha_TF[j] = n * (D_zz_alpha[j] - h_zz * h_WW[j] * D_ww_alpha[j])  / (D - 1.);
    D_ww_alpha_TF[j] =  (D_ww_alpha[j] - h_ww * h_ZZ[j] * D_zz_alpha[j])  / (D - 1.);

    //2nd conformal derivative of chi wrt z\/
    double cD_zz_chi = d_zz(v_chi,j) - chris_Zzz[j] * d_z_chi;

    double d_chi_z = ((z <= min_z) ? d_zz(v_chi,j) : (d_z_chi / z)); //d_z(chi) / z replaced with 2nd deriv at 0
    double c_chris_Z_overz = ((z <= min_z) ? d_z(v_c_chris_Z,j) : (c_chris_Z / z)); // c_chris_Z / z replaced with 2nd deriv at 0
    double d_h_ww_z = ((z <= min_z) ? d_zz(v_h_ww, j) : (d_z_h_ww / z)); //d_z(h_ww) / z replaced with with 2nd deriv at 0

    //conformal + chi parts of the Ricci tensor components
    double R_zz_chi = n * ( cD_zz_chi + (0.5 * h_WW[j] * d_z_h_ww * d_z_chi + d_chi_z ) - d_z_chi * d_z_chi / chi   ) / (2. * chi);

    double R_ww_chi = h_ww * h_ZZ[j] * (cD_zz_chi + (2. * D - 5.) * (0.5 * h_WW[j] * d_z_h_ww * d_z_chi + d_chi_z) - (D - 1.) * d_z_chi * d_z_chi / (2. * chi) ) / (2. * chi);

    double R_zz_c_t1 = ((z <= min_z) ?  ( -0.5 * d_zz(v_h_ww,j) ): ((h_zz - h_ww) / (z * z) - 0.5 * d_z_h_zz / z)  ); // first bracketed term in R_zz_c
    double h_zw_diff = ((z <= min_z) ?  ( 0.5 * (d_zz(v_h_zz,j) - d_zz(v_h_ww,j) )) : ((h_zz - h_ww) / (z * z)) ); // (h_zz -h_ww) / z^2 replaced by half 2nd deriv difference at z = 0

        /*double R_zz_c = n * h_WW[j] * (R_zz_c_t1 - 0.25 * h_WW[j] * d_z_h_ww * d_z_h_ww )
                        - 0.5 * h_ZZ[j] * d_zz(v_h_zz,j) + h_zz * d_z(v_c_chris_Z,j) + c_chris_Z * chris_zzz + 3. * h_ZZ[j] * chris_zzz * chris_Zzz[j];
        if (z > min_z)
            R_zz_c += n * h_WW[j] * (h_WW[j]*h_zz - 1.) * d_z_h_ww / z;*/

    double R_zz_c = n * h_WW[j] * (R_zz_c_t1 + chris_Wwz * chris_wwz + 2 * chris_Wwz * chris_zww)
                        - 0.5 * h_ZZ[j] * d_zz(v_h_zz,j) + h_zz * d_z(v_c_chris_Z,j) + c_chris_Z * chris_zzz + 3. * h_ZZ[j] * chris_zzz * chris_Zzz[j];

   /* double R_ww_c = -0.5 * h_ZZ[j] * d_zz(v_h_ww,j) + 0.5 * h_WW[j] * h_ZZ[j] * d_z(v_h_ww, j) * d_z(v_h_ww, j) - 0.5 * n * h_WW[j] * d_h_ww_z
                        + h_ww * c_chris_Z_overz + 0.5 * c_chris_Z * d_z_h_ww - h_ZZ[j] * h_zw_diff;*/

    double R_ww_c = -0.5 * h_ZZ[j] * d_zz(v_h_ww,j) + 3 * h_ZZ[j] * chris_Wwz * chris_wwz - 0.5 * n * h_WW[j] * d_h_ww_z
                        + h_WW[j] * (chris_Zww[j] * chris_zww + 2. * chris_Zww[j] * chris_wwz) + h_ww * c_chris_Z_overz +  c_chris_Z * chris_wwz - h_WW[j] * h_zw_diff;


    R_zz[j] = R_zz_c + R_zz_chi;// d_z(v_chi,j);
    R_ww[j] = R_ww_c + R_ww_chi;

    if (use_CCZ4) //if CCZ4 is used, add supplementary terms to the (extended) Ricci tensor and fill out theta_Z
    {
        theta_Z[j] = 0.5 * chi * (c_chris_Z - h_ZZ[j] * chris_Zzz[j] - n * h_WW[j] * chris_Zww[j]);
        R_zz[j] += theta_Z[j] * (chi * d_z_h_zz + h_zz * d_z_chi) / (chi * chi) - 2. * theta_Z[j] * chris_zzz / chi; //R^Z_zz with correction to c_chris_Z (tilde->hat)
        R_ww[j] += theta_Z[j] * (chi * d_z_h_ww - h_ww * d_z_chi) / (chi * chi) - 2. * theta_Z[j] * chris_wwz / chi; // R^Z_ww with correction to c_chris_Z (tilde->hat)
    }

    //traceless part of Ricci components
    R_zz_TF[j] = n * (R_zz[j] - h_zz * h_WW[j] * R_ww[j]) / (D - 1.);
    R_ww_TF[j] = (R_ww[j] - h_ww * h_ZZ[j] * R_zz[j] ) / (D - 1.);

    //|phi|
    double mod_phi = sqrt(phi_re * phi_re + phi_im * phi_im);

    if (j == 0) test_ctr = mod_phi;

    //matter quantities: some may only be valid in 3+1; check if we want to do D != 4
    rho[j] = 2. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) +  0.5 * h_ZZ[j] * chi * (pow(d_z_phi_re,2) + pow(d_z_phi_im,2)) + 0.5 * V(mod_phi);
    j_z[j] = 2. * (K_phi_re * d_z_phi_re + K_phi_im * d_z_phi_im );

    S_zz[j] = 0.5 * (pow(d_z_phi_re,2) + pow(d_z_phi_im,2)) - 0.5 * h_zz * (V(mod_phi) - 4. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) ) / chi;
    S_ww[j] = - 0.5 * h_ww * (h_ZZ[j] * (pow(d_z_phi_re,2) + pow(d_z_phi_im,2)) + (V(mod_phi) - 4. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re)) / chi );

    //convenient to use identity involving rho in 3+1; otherwise compute trace normally
    if (D == 4.)
        S[j] = 8. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) - V(mod_phi) - rho[j];
    else
        S[j] = 0.5 * (3. - D) * h_ZZ[j]  * chi * (pow(d_z_phi_re,2) + pow(d_z_phi_im,2))  - 0.5 * (D - 1.) * (V(mod_phi) - 4. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re));

    S_zz_TF[j] = S_zz[j] - S[j] * h_zz /( chi * (D - 1.));
    S_ww_TF[j] = S_ww[j] - S[j] * h_ww/ ( chi * (D - 1.));

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
    rhs.states.resize(n_gridpoints);

    int max_ref = (refinement_levels.size() == 0 ) ? 1 : refinement_levels[refinement_levels.size() - 1]; //refinement level at outermost bdry

    for (int j = 0; j < n_gridpoints - 2; j++)
    {
        if (!active_points[j])//perform no calculations and skip on inactive points
            continue;//rhs.states[j] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}

    double z = (j + grid_offset) * dr;

        //shorthand versions of the BSSN vars on desired slice for convenience
        const double chi = max(min_chi, slice_ptr->states[j].chi);
        const double h_zz = slice_ptr->states[j].h_zz;
        const double h_ww = slice_ptr->states[j].h_ww;
        const double A_zz = slice_ptr->states[j].A_zz;
        const double A_ww = slice_ptr->states[j].A_ww;
        const double K = slice_ptr->states[j].K;
        const double c_chris_Z = slice_ptr->states[j].c_chris_Z;
        const double phi_re = slice_ptr->states[j].phi_re;
        const double phi_im = slice_ptr->states[j].phi_im;
        const double K_phi_re = slice_ptr->states[j].K_phi_re;
        const double K_phi_im = slice_ptr->states[j].K_phi_im;
        const double alpha = slice_ptr->states[j].alpha;
        const double beta = slice_ptr->states[j].beta;


        if (wave_mode) //converts to wave eq'n solver with beta the time derivative of phi_re for testing purposes
        {
            rhs.states[j] = (BSSNState){0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            rhs.states[j].phi_re = beta;
            rhs.states[j].beta = d_zz(v_phi_re, j);
            continue;
            //return rhs;
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

        //fill out RHS of the BSSN evolution system
        rhs.states[j].chi  = beta * d_z_chi  + 2.* chi * (alpha * K - d_z_beta - n * beta_z ) / (D - 1.);

        rhs.states[j].h_zz = beta * d_z_h_zz + 2. * h_zz * d_z_beta - 2. * h_zz * (d_z_beta + n * beta_z) / (D - 1.) - 2. * alpha * A_zz;
        rhs.states[j].h_ww = beta * d_z_h_ww - 2. * h_ww * (d_z_beta - beta_z) / (D - 1.) - 2. * alpha * A_ww;

        if (!use_CCZ4) //different evolution eq'n used in CCZ4
        {rhs.states[j].K = beta * d_z_K - chi * h_ZZ[j] * D_zz_alpha[j] + alpha * h_ZZ[j] * h_ZZ[j] * A_zz * A_zz + alpha * K * K / (D - 1.)
                        + n * h_WW[j] * (alpha * A_ww * A_ww / h_ww - chi * D_ww_alpha[j]) + 8. * M_PI * alpha * (S[j] + (D - 3.) * rho[j]) / n;}

        rhs.states[j].A_zz = beta * d_z(v_A_zz, j) + 2. * A_zz * d_z_beta - 2. * A_zz * (d_z_beta + n * beta_z) / (D - 1.) + alpha * K * A_zz
                            - 2. * alpha * A_zz * h_ZZ[j] * A_zz + chi * (alpha * R_zz_TF[j] - D_zz_alpha_TF[j] - 8. * M_PI * alpha * S_zz_TF[j]);
        rhs.states[j].A_ww = beta * d_z(v_A_ww, j) - 2. * A_ww * (d_z_beta -  beta_z) / (D - 1.) + alpha * A_ww * (K - 2 * h_WW[j] * A_ww)
                            + chi * (alpha * R_ww_TF[j] - D_ww_alpha_TF[j] - 8. * M_PI * alpha * S_ww_TF[j]);

        rhs.states[j].c_chris_Z = beta * d_z(v_c_chris_Z,j) + 2. * c_chris_Z * (d_z_beta + n * beta_z) / (D - 1.) + h_ZZ[j] * d_zz(v_beta,j) - c_chris_Z * d_z_beta
                                + (D - 3.) * h_ZZ[j] * d_zz(v_beta,j) / (D - 1.) - 2. * (D - 2.) * alpha * h_ZZ[j] * d_z_K / (D - 1.)
                                - A_zz * h_ZZ[j] * h_ZZ[j] * ( (D - 1.) * alpha * d_z_chi / chi + 2 * d_z_alpha)
                                + 2 * alpha * (chris_Zzz[j] * h_ZZ[j] * h_ZZ[j] * A_zz + n * chris_Zww[j] * h_WW[j] * h_WW[j] * A_ww)
                                - sigma_BSSN * aux_constraint[j] - 16. * M_PI * alpha * j_z[j]* h_ZZ[j];

        if (z >= min_z) //add terms in d_z(beta) / z - beta /z^2 when z =/= 0
            rhs.states[j].c_chris_Z += n * (h_WW[j] + (D - 3.) * h_ZZ[j] / (D - 1.)) * (d_z_beta / z - beta / (z * z));

        //|phi|
        double mod_phi = sqrt(phi_re * phi_re + phi_im * phi_im);

        rhs.states[j].phi_re = beta * d_z_phi_re - 2. * alpha * K_phi_re;
        rhs.states[j].phi_im = beta * d_z_phi_im - 2. * alpha * K_phi_im;

        double d_phi_re_z = ((z <= min_z) ? d_zz(v_phi_re,j) : d_z_phi_re / z); //dphi(z) / z replaced by 2nd deriv at 0
        double d_phi_im_z = ((z <= min_z) ? d_zz(v_phi_im,j) : d_z_phi_im / z);

        //conformal 2nd covariant derivatives of scalar field
        double cD_zz_phi_re = d_zz(v_phi_re, j ) - chris_Zzz[j] * d_z_phi_re ;
        double cD_zz_phi_im = d_zz(v_phi_im, j ) - chris_Zzz[j] * d_z_phi_im ;

        double cD_ww_phi_re = d_phi_re_z - chris_Zww[j] * d_z_phi_re;
        double cD_ww_phi_im = d_phi_im_z - chris_Zww[j] * d_z_phi_im;


        rhs.states[j].K_phi_re = beta * d_z(v_K_phi_re,j) + alpha * K * K_phi_re + 0.5 * alpha * phi_re * dV(mod_phi)
                               - 0.5 * chi * (h_ZZ[j] * d_z_alpha * d_z_phi_re + alpha * (h_ZZ[j] * cD_zz_phi_re + n * h_WW[j] * cD_ww_phi_re))
                               + 0.25 * (D - 3.) * alpha * h_ZZ[j] * d_z_chi * d_z_phi_re;

        rhs.states[j].K_phi_im =  beta * d_z(v_K_phi_im,j) + alpha * K * K_phi_im + 0.5 * alpha * phi_im * dV(mod_phi)
                               - 0.5 * chi * (h_ZZ[j] * d_z_alpha * d_z_phi_im + alpha * (h_ZZ[j] * cD_zz_phi_im + n * h_WW[j] * cD_ww_phi_im))
                               + 0.25 * (D - 3.) * alpha * h_ZZ[j] * d_z_chi * d_z_phi_im;

        //Gauge variable update using moving puncture evolution, unless evolve_shift is off in which case do not update beta

        double f = (shock_gauge) ? (1. + shock_fac / (alpha * alpha) ): (one_log_fac/ alpha); //bona-masso f
        //cout<< one_log_fac << endl;

        rhs.states[j].alpha = beta * d_z_alpha - f * pow(alpha, 2.) * K;
        rhs.states[j].beta =  (evolve_shift)? (beta * d_z_beta + gamma_fac * c_chris_Z - eta * beta): 0.;

        //add supplementary terms if CCZ4 evolution is on
        if (use_CCZ4)
        {
            const double theta = slice_ptr->theta[j];
            rhs.theta.resize(n_gridpoints);

            //rhs.states[j].K += alpha * (-2. * K * theta - (D - 1.) * c1 * (1. + c2) * theta); //FIX ME!!!

            //should add correction to above for BSSN to avoid recomputing
            rhs.states[j].K = beta * d_z_K - chi * h_ZZ[j] * D_zz_alpha[j]  - n * h_WW[j] * chi * D_ww_alpha[j]
                            + alpha * (chi * h_ZZ[j] * R_zz[j] + n * chi * h_WW[j] * R_ww[j] + K * (K - 2. * theta)
                            - (D - 1.) * c1 * (1. + c2) * theta + 8. * M_PI * (S[j] - (D - 1.) * rho[j]) / (D - 2.));

            rhs.states[j].A_zz += -2. * alpha * theta * A_zz;
            rhs.states[j].A_ww += -2. * alpha * theta * A_ww;

            //note this is now \hat{\Gamma}^z, not \tilde{\Gamma}^z. So we first remove BSSN specific damping then add new terms.
            rhs.states[j].c_chris_Z +=   sigma_BSSN * aux_constraint[j] - 2. * d_z_alpha * theta * h_ZZ[j] + 2. * alpha * h_ZZ[j] * d_z(v_theta, j)
                                        - 2. * c1 * alpha * theta_Z[j] / chi - 4. * alpha * K * theta_Z[j] / ((D - 1.) * chi)
                                        + sigma_BSSN * theta_Z[j]  * ((3. - D) * d_z_beta + 2. * n * beta_z) / chi;

            rhs.states[j].alpha +=  2. * alpha * alpha * f * theta; //4. * alpha * theta; //corrects for K -> K - 2 * theta in 1+log slicing
            //rhs.states[j].beta +=  (evolve_shift)? ( -2. * gamma_fac * theta_Z[j] / chi): 0.; //replace theta hat -> theta tilde in CCZ4 moving puncture gauge; seems to produce instability...

            if (!theta_off)
                rhs.theta[j] = beta * d_z(v_theta, j) + 0.5 * alpha * (chi * h_ZZ[j] * R_zz[j] - h_ZZ[j] * h_ZZ[j] * A_zz * A_zz + n * K * K / (D - 1.)
                            + n * (chi * h_WW[j] * R_ww[j] - A_ww * A_ww / (h_ww * h_ww))  - 2. * K * theta - 2. * theta_Z[j] * d_z_alpha / alpha
                            - c1 * (D + c2 * n) * theta  - 16. * M_PI * rho[j]);
            else
                rhs.theta[j] = 0.;
        }

        //add damping in away from edges for now; may need to add for edges-- we'll see
        if (damping_factor != 0. && j < n_gridpoints - 3 * pow(2, max_ref - 1))
            {
                int res_fac = pow(2, slice_ptr->get_refinement_level(j, refinement_points) - 1);

                vector<int> J = {j - 3 * res_fac, j - 2 * res_fac, j - 1 * res_fac, j, j + 1 * res_fac, j + 2 * res_fac, j + 3 * res_fac}; //indices at which to take stencil

                if (j < 3) //use symmetry across 0 to fill 7-point stencils
                {
                    for (int& index: J)
                        index = std::max(-index, index);
                }

                //pointers to state information at stencil location
                vector<BSSNState*> sts = {&slice_ptr->states[J[0]], &slice_ptr->states[J[1]], &slice_ptr->states[J[2]], &slice_ptr->states[J[3]], &slice_ptr->states[J[4]], &slice_ptr->states[J[5]], &slice_ptr->states[J[6]]};

                //account for odd parity of beta and contracted christoffels-- seems bad though?
                /*if (j < 3)
                {
                    for (int k = 0; k < 3 - j; k++)
                    {
                        BSSNState s_inner = *sts[k];
                        s_inner.beta *= -1.;
                        s_inner.c_chris_Z *= -1.;
                        *(sts[k]) = &s_inner;
                    }
                }*/

                double d_mult = damping_factor * (pow(dr * res_fac, 5.) / 64. );
                BSSNState damping_corr = d_mult * sevenPointDeriv(dr * res_fac, 6, *sts[0], *sts[1], *sts[2], *sts[3],*sts[4], *sts[5], *sts[6]);

                rhs.states[j] = rhs.states[j] + damping_corr;

                if (use_CCZ4)
                {
                    vector<double*> ths = {&slice_ptr->theta[J[0]], &slice_ptr->theta[J[1]], &slice_ptr->theta[J[2]], &slice_ptr->theta[J[3]], &slice_ptr->theta[J[4]], &slice_ptr->theta[J[5]], &slice_ptr->theta[J[6]]};
                    double theta_corr = d_mult * sevenPointDeriv<double>(dr * res_fac, 6, *ths[0], *ths[1], *ths[2], *ths[3],*ths[4], *ths[5], *ths[6]);
                    rhs.theta[j] = rhs.theta[j] + theta_corr;

                }
            }
    }

    int res_fac = pow(2, max_ref - 1); //number of points skipped per slot in stencil

    //asymptotic states and their derivatives where relevant (can ignore for matter values as they decay exponentially)
    BSSNState asymp_state{chi_asymp(R - dr * res_fac), 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.,alpha_asymp(R - dr * res_fac), 0.};
    BSSNState asymp_deriv{d_chi_asymp(R - dr * res_fac), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., d_alpha_asymp(R - dr * res_fac), 0.}; // r-derivative of asymptotic expansion

    //version without spatially varying asymptotics
    if (!spatially_varying_BC)
    {
        asymp_state = (BSSNState) {1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.,1., 0.};
        asymp_deriv = (BSSNState) {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    }

    //maybe adjust to account for purported 1/r^2 decay in K!
    BSSNState N{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    //N = pow(R, (4. - D))  * N;
    N = pow(R, (4. - D))  * N;


    //outermost 5 active gridpoints
    BSSNState& s1 = slice_ptr->states[last_active_j - 4 * res_fac];
    BSSNState& s2 = slice_ptr->states[last_active_j - 3 * res_fac];
    BSSNState& s3 = slice_ptr->states[last_active_j - 2 * res_fac];
    BSSNState& s4 = slice_ptr->states[last_active_j - 1 * res_fac];
    BSSNState& s5 = slice_ptr->states[last_active_j];

    //limiting characteristic speeds at infinity
    //BSSNState char_speeds{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
    BSSNState char_speeds{sqrt(2. / s4.alpha), 1., 1., 1., 1., sqrt(2. / s4.alpha), 1., 1., 1., 1., 1., sqrt(2. / s4.alpha), 1.};

    rhs.states[last_active_j - res_fac] = (-1.) * char_speeds * ( p4_stencil(dr * res_fac, s1, s2, s3, s4, s5) + (0.5 * D - 1.) * N * (s4 - asymp_state) /  (dr * (last_active_j - res_fac))- asymp_deriv);

    //update variable asymptotic states to outermost edge
    if (spatially_varying_BC)
    {asymp_state.chi = chi_asymp(R);  asymp_state.alpha = alpha_asymp(R);
    asymp_deriv.chi = d_chi_asymp(R);  asymp_deriv.alpha = d_alpha_asymp(R);}

    rhs.states[last_active_j] = (-1) * char_speeds * ( p5_stencil(dr * res_fac, s1, s2, s3, s4, s5) + (0.5 * D - 1.) * N * (s5 - asymp_state) / (dr * last_active_j) - asymp_deriv);


    if (use_CCZ4)
    {
        double& t1 = slice_ptr->theta[last_active_j - 4 * res_fac];
        double& t2 = slice_ptr->theta[last_active_j - 3 * res_fac];
        double& t3 = slice_ptr->theta[last_active_j - 2 * res_fac];
        double& t4 = slice_ptr->theta[last_active_j - 1 * res_fac];
        double& t5 = slice_ptr->theta[last_active_j];

        rhs.theta[last_active_j - res_fac] = -(p4_stencil(dr * res_fac, t1, t2, t3, t4, t5) + (0.5 * D - 1.) * t4 / (dr * (last_active_j - res_fac)));
        rhs.theta[last_active_j] = -(p5_stencil(dr * res_fac, t1, t2, t3, t4, t5) + (0.5 * D - 1.) * t5 / (dr * (last_active_j)));
    }

    return rhs;
}

//enforce tracelessness of A
void Spacetime::make_A_traceless(BSSNSlice* slice_ptr)
{
    double n = D - 2.;

    for (int j = 0; j < n_gridpoints; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

        const double h_zz = slice_ptr->states[j].h_zz;
        const double h_ww = slice_ptr->states[j].h_ww;
        const double A_zz = slice_ptr->states[j].A_zz;
        const double A_ww = slice_ptr->states[j].A_ww;

        double A = A_zz / h_zz + n * A_ww / h_ww;
        slice_ptr->states[j].A_zz = A_zz - h_zz * A / (D - 1.);
        slice_ptr->states[j].A_ww = A_ww - h_ww * A / (D - 1.);
    }
}

//TODO: make these work for D =/= 4, and MR if this gets working
double Spacetime::slice_mass(BSSNSlice* slice_ptr)
{
    const double n = D - 2;
    double mass = 0;

    for (int j = 0; j <  0.95 * n_gridpoints - 1; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;
        const double& chi = slice_ptr->states[j].chi;

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

        const double& chi = slice_ptr->states[j].chi;
        const double& phi_re = slice_ptr->states[j].phi_re;
        const double& phi_im = slice_ptr->states[j].phi_im;
        const double& K_phi_re = slice_ptr->states[j].K_phi_re;
        const double& K_phi_im = slice_ptr->states[j].K_phi_im;

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
        double chi = slice_ptr->states[j].chi;
        const double h_zz = slice_ptr->states[j].h_zz;
        const double h_ww = slice_ptr->states[j].h_ww;
        const double A_zz = slice_ptr->states[j].A_zz;
        const double A_ww = slice_ptr->states[j].A_ww;
        const double K = slice_ptr->states[j].K;
        //const double beta = slice_ptr->states[j].beta;

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

    data_file <<  std::setprecision (16) << (j + grid_offset) * dr << "   " << s.states[j].chi << "    " << s.states[j].h_zz << "    " << s.states[j].h_ww  << "    " << s.states[j].A_zz
        << "   " << s.states[j].A_ww << "    " << s.states[j].K << "    " << s.states[j].c_chris_Z  << "    " << s.states[j].phi_re << "    " << s.states[j].phi_im
        << "   " << s.states[j].K_phi_re << "    " << s.states[j].K_phi_im << "    " << s.states[j].alpha << "    " << s.states[j].beta << "    " << theta0 << endl;
    }
}

//writes diagnostic quantities to output file
void Spacetime:: write_diagnostics()
{
    int length = (current_slice_ptr->states).size();

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

        const double& phi_re = current_slice_ptr->states[j].phi_re;
        const double& phi_im = current_slice_ptr->states[j].phi_im;

        //const double& chi = current_slice_ptr->states[j].chi;
        //const double& c_chris_Z = current_slice_ptr->states[j].c_chris_Z;

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

        fill_param_array(current_line, "refinement_points = ", refinement_points, quiet);

    }

    cout << sigma_BSSN << eta << endl;
}

//halves number of gridpoints while keeping R. Should only be called before evolving!
void Spacetime::halve_resolution()
{
    int n_old = n_gridpoints;

    vector<BSSNState> old_states(n_old);

    //cut off one point at end to get odd # of gridpoints if needed
    if (n_gridpoints % 2 == 0)
        {
            n_gridpoints -= 1;
            R *= (n_gridpoints) / (n_old);
        }

    for (int j = 0; j < n_gridpoints; j++) //fill old array with states
        old_states[j] = slices[0].states[j];

    n_gridpoints = (n_gridpoints + 1) / 2;

    for (int k = 0; k < n_gridpoints; k++)
        slices[0].states[k] = old_states[2 * k];

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

        //average at points that seem to produce noise to hopefully kill error. should be only 2nd-order...
        BSSNSlice& s1 = *current_slice_ptr;
        s1.states[j1] = 0.5 * (s1.states[j1 - step] + s1.states[j1 + step]);
        s1.states[j3] = 0.5 * (s1.states[j3 - step] + s1.states[j3 + step]);

        start_point = j4;

        //cout << "Success on level " << level << endl;

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
        //old_momenta[j] = slices[0].states[j].K_phi_im;
        double correction_sq = slices[0].states[j].K_phi_im * slices[0].states[j].K_phi_im + Ham[j] / (32. * M_PI);

        if (correction_sq < 0. && failed == 0)
            {
                failed = 1;
                cout << "WARNING: trick to adjust initial field momentum may have failed starting on gridpoint " << j << endl;
            }

        if (correction_sq >= 0.)
             slices[0].states[j].K_phi_im = - sqrt(correction_sq);
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
            double local_dtK =  (slice_current.states[j].K - slice_last.states[j].K ) / dt;
            dtK_L2 += dr * local_dtK * local_dtK * z * z   * pow(slice_current.states[j].chi, -1.5);
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
        double phase = std::arg( std::complex<double>(s.states[j].phi_re, s.states[j].phi_im));
        //double K_phase = std::arg( std::complex<double>(s.states[j].K_phi_re, s.states[j].K_phi_im));
        double r = j * dr;

        s.states[j].phi_re += cos(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states[j].phi_im += sin(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states[j].K_phi_re += 0.5 * sin(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states[j].alpha);
        s.states[j].K_phi_im -= 0.5 * cos(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states[j].alpha);
    }

    solve_initial_ham();
}

//read data from BosonStar to spacetime and construct initial time slice
void Spacetime::initialize(BosonStar& boson_star)
{
    wave_mode = 0; //make 1 for MR testing purposes only
    //D = SPACEDIM + 1.;
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

    gamma_fac = 0.75; // initialize to 3/4 by default (before params read) for backwards compatibility w/ older params files
    spatially_varying_BC = 1; //as above
    one_log_fac = 2.0; //as above

    refinement_points = {};
    read_parameters();

    if (D != 4.)
        spatially_varying_BC = 0; //only try this for D = 4 for now

    //cout << refinement_points.size() << endl;

    if (refinement_points[0] <= 0) //signal to disable any refinement and use all points
         refinement_points.clear();//refinement_points = {};

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

    //solve BS at higher resolution and read in data to first slice, if not starting in other mode (checkpoint / thinshell read / gaussian start)
    if (BS_resolution_factor > 1 && start_time == 0 && !read_thinshell && !boson_star.gaussian_start)
    {
        boson_star.n_gridpoints = boson_star.n_gridpoints * BS_resolution_factor - BS_resolution_factor + 1;

        cout << " \nRe-solving BS at higher resolution with " << boson_star.n_gridpoints <<  endl;

        boson_star.solve();
        boson_star.write_field();
        boson_star.fill_isotropic_arrays();
        boson_star.write_isotropic();

        //using higher-res omega seems like an improvement
        omega = boson_star.omega;
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

    cout << "About to read" << endl;

    if (start_time == 0)
        slices[0].read_BS_data(boson_star, BS_resolution_factor, isotropic, cell_centered);

    grid_offset = 0.5 * (double)cell_centered;

    cout << "Read BS data" << endl;

    if (run_spacetime_solver && start_time == 0)
        solve_initial_ham();

    //cut off outermost 2 gripoints, where the christoffel symbols will be generally polluted by garbage due to not having data to take derivatives there. Might be better to just extrapolate long term.
    n_gridpoints -= 2 ;

    R *= (n_gridpoints - 1.0) / (n_gridpoints + 1.0); //also need to rescale R to avoid stretching solution
    slices[0].R = R;

    int n_old = n_gridpoints;
    n_gridpoints = round(cutoff_frac * n_gridpoints); //shrink domain by cutoff_frac, ideally to remove detritus in H
    R = (R * (n_gridpoints - 1.0)) / (n_old - 1.0); //TODO: think about this

    slices[0].states.resize(n_gridpoints);
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

    if (make_tangherlini || slices[0].states[0].chi < 10 * min_chi)
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
    cout << "Dynamical N =" << charge0 << ", M = "  << mass0 << " and binding energy is " << mass0 - mu * charge0  << endl;

}


//Time evolution! Must have initialized using a constructed BosonStar first.
void Spacetime::evolve()
{
    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    cout << "dr = " << dr << ", dt = " << dt <<  "   " << endl;

    int num_timesteps = ceil(stop_time / dt);
    int last_checkpoint_time = 0;

    // if(store_A0) A0_values.resize(num_timesteps);

    //write constraint norms at each timestep to file
    std::ofstream constraints_file{"constraint_norms.dat"};
    if (!constraints_file)
    {
        std::cerr << "constraint_norms.dat could not be opened for writing!\n";
        exit(1);
    }

    cout <<" \n Will evolve with " << num_timesteps << " time steps \n" << endl;

    //s_i are returned RHS's, t represents temporary RHS + current_slice quantities that must be stored so derivatives can be accessed
    BSSNSlice s1, s2, s3, s4, t1, t2, t3;

    for (int time_step = 0; time_step < num_timesteps; time_step++)
    {
        double t = start_time + time_step * dt;

        //fill out array until we've reached maximum number of stored slices, then update last element + rotate at end.
        int n = (time_step > max_stored_slices - 2) ? (max_stored_slices - 2) : time_step;

        //if(store_A0)
            //A0_values[time_step] = test_ctr;

        double phase_ctr, phase_last;
        if (time_step > 0) phase_last = phase_ctr;
        phase_ctr = std::arg( std::complex<double>(slices[n].states[0].phi_re, slices[n].states[0].phi_im));

        double omega_approx = (phase_ctr - phase_last) / dt;

        if (omega_approx < -M_PI / dt) omega_approx += 2 * M_PI / dt; //corrects jumps due to branch cuts
        double phase_diff= phase_ctr - std::fmod(omega * dt * time_step, M_PI ); //difference between actual and expected phase at center of BS
        if (phase_diff < 0.) phase_diff += M_PI;

        if (n > 0) compute_dtK(n);
        M = slice_mass(current_slice_ptr); // maybe remove if causes bad bdry oscillations?

        double A0 = sqrt (pow(slices[n].states[0].phi_re,2) + pow(slices[n].states[0].phi_im,2));
        double A1 = sqrt (pow(slices[n].states[1].phi_re,2) + pow(slices[n].states[1].phi_im,2));

        double A_ctr =  (cell_centered) ? (A0- 0.5 * dr * (A1 - A0)) : A0; //just linearly extrapolate to r = 0 for now

        if (time_step % write_CN_interval == 0) //write time-dependent diagnostics to constraint_norms.dat
            constraints_file << std::setprecision (10) << start_time + dt * time_step << "   " << Ham_L2  << "   " << Mom_L2 <<  "   " << slices[n].states[0].chi << "   "
            << A_ctr << "   "  << phase_diff   << "   " <<  M << "   "  << slice_charge(current_slice_ptr)
            << "   " << dtK_L2 << "   " << omega_approx << "   " << slices[n].states[0].alpha<<  endl;

        slices[n + 1].states.resize(n_gridpoints);
        slices[n + 1].R = R;
        slices[n + 1].refinement_points = refinement_points;
        slices[n + 1].active_points = active_points;

        if (make_tangherlini || slices[n].states[0].chi < 10 * min_chi)
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
        kill_refinement_noise(); //running this on every timestep appears to be best approach...
        make_A_traceless(current_slice_ptr);

        //enforce minimum chi
        for (BSSNState& s: slices[n + 1].states)
            {if (s.chi < min_chi) s.chi = min_chi; }

        if (time_step % write_interval == 0)
        {
            //current_slice_ptr->write_slice();
            write_current_slice();
            write_diagnostics();
        }

        //write checkpoint files
        //TODO make this work for pre C++17

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

        if (isnan(slices[n + 1].states[0].chi))
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
    }

}


#endif /*EVOLUTIONVARIABLES*/
