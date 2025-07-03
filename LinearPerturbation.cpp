#ifndef LINEARPERTURBATION_CPP_
#define LINEARPERTURBATION_CPP_

#include "BosonStar.h"
#include "LinearPerturbation.h"
#include "mathutils.h"
#include <sstream>
#include<iomanip>

using namespace std;


//overloads for adding / scalar multiplication of PertState objects
PertState operator+(const PertState& s1, const PertState& s2)
{
    return (PertState){s1.F + s2.F, s1.Fp + s2.Fp, s1.L + s2.L, s1.Lp + s2.Lp};
}

PertState operator-(const PertState& s1, const PertState& s2)
{
    return (PertState){s1.F - s2.F, s1.Fp - s2.Fp, s1.L - s2.L, s1.Lp - s2.Lp};
}

PertState operator*(double c, const PertState& s)
{
    return (PertState){c * s.F, c * s.Fp, c * s.L, c * s.Lp};
}

//returns rhs of the radial ODEs that F, L satisfy
PertState LinearPerturbation::pert_rhs(double r, PertState s, double chi_sq, long double gamma)
{
    if (r == 0.)
    {
        double Fpp0 = ( 2. + (2. * omega * omega - chi_sq) / exp(2. * bg->state[0].phi) + solitonic * (72. * pow(A_central / sigma, 4)  //the small-r expansion leaves L''(0), here gamma, undetermined--
                      - 32. * pow(A_central / sigma, 2))) / 3. - gamma / (8. * M_PI * A_central * A_central); // we keep it arbitrary but express F''(0) in terms of it.

        return (PertState) {0, Fpp0, 0, gamma};
    }

    int j = floor(r / dr); //smallest index below r

    double A0, X0, alpha0, Ap0; //background values of A, X, alpha, and A' (via BS eta) at r

    // detect mid-gridpoint r-values and use linear interpolation of bg quantities where needed.
    if ( abs(j * dr - r) < 0.1 * dr )
        {A0 = bg->state[j].A; X0 = bg->state[j].X; alpha0 = exp(bg->state[j].phi); Ap0 = X0 * (bg->state[j].eta);}
    else if (abs((j + 1) * dr - r) < 0.1 * dr )
        {A0 = bg->state[j + 1].A; X0 = bg->state[j + 1].X; alpha0 = exp(bg->state[j + 1].phi); Ap0 = X0 * (bg->state[j + 1].eta);}
    else
        {
            A0 = 0.5*( bg->state[j].A + bg->state[j + 1].A) ; X0 = 0.5*( bg->state[j].X + bg->state[j + 1].X);
            alpha0 = exp(0.5*( bg->state[j].phi + bg->state[j + 1].phi)); Ap0 = X0 * 0.5 *( bg->state[j].eta + bg->state[j + 1].eta);
        }

    double X0_sq = X0 * X0;
    double A0_sq = A0 * A0;
    double omega_sq = omega * omega;
    double alpha0_sq = alpha0 * alpha0;
    double sigma_sq = sigma * sigma;

    //now compute background derivatives using respective bg ODEs
    double Xp0 = X0 * ( 0.5 * (1 - X0_sq) / r + 2. * M_PI * r * X0_sq * ( (Ap0 * Ap0) / (X0_sq) + pow(omega * A0 / alpha0, 2) + bg->V(A0)  ));
    double alphap0 = alpha0 * ( 0.5 * (-1 + X0_sq) / r + 2. * M_PI * r * X0_sq * ( (Ap0 * Ap0) / (X0_sq) + pow(omega * A0 / alpha0, 2) - bg->V(A0)  ));


    //for now: interpolate 3 surrounding points at r + dr using cubic legendre + use 2nd order approx. to get X''. Can almost certainly do better!!!
    int j0 = bound(j - 1, 0, n_gridpoints - 4);

    vector<double> radii{dr * j0, dr *(j0 + 1), dr * (j0 + 2), dr * (j0 + 3)};
    vector<double> XPP{bg->state[j0].X, bg->state[j0 + 1].X, bg->state[j0 + 2].X, bg->state[j0 + 3].X};

    //double Xpp1 = lagrange_interp(r - dr, radii, XPP); double Xpp2 = lagrange_interp(r, radii, XPP); double Xpp3 = lagrange_interp(r + dr, radii, XPP);
    //double Xpp0 = (Xpp1 + Xpp3 - 2. * Xpp2) / (dr * dr);


    //EXPERIMENTAL: analytic expression for the background value of X''
    double Xpp0 = 3.*Xp0*Xp0 / X0 - 2. * Xp0 / r + 4. * M_PI * (A0_sq + 4.*solitonic *(A0_sq * A0_sq * A0_sq/(sigma_sq*sigma_sq) - A0_sq * A0_sq/sigma_sq )
           + 2. * r * A0 * Ap0 + omega_sq * A0_sq * (1 - r * alphap0 / alpha0) / alpha0_sq  + r * Ap0 * A0 *
           solitonic * (24. * A0_sq*A0_sq / (sigma_sq*sigma_sq) - 16. * A0_sq / sigma_sq ) ) * X0_sq * X0
           - 4. * M_PI * Ap0*Ap0 * (1. + r * alphap0 / alpha0) * X0;

    const double& F = s.F;
    const double& Fp = s.Fp;
    const double& L = s.L;
    const double& Lp = s.Lp;

    //the actual pulsation equations
    double Fpp = ((2. + ( 2. * omega_sq - chi_sq ) / alpha0_sq + 8. * A0 * Ap0 * r * M_PI  + 8. * solitonic * A0_sq
               *(r * A0 * Ap0 * M_PI * (12. * A0_sq - 8. * sigma_sq) + 9. * A0_sq - 4 * sigma_sq) / pow(sigma_sq, 2)  ) * X0_sq + 2. * pow(Ap0 / A0, 2 )) * F
               + ((1 - omega_sq / (alpha0_sq) + solitonic * (12. * pow(A0_sq / sigma_sq, 2) - 8. *  A0_sq / sigma_sq)) * X0_sq
               - 1. / (4. * M_PI * pow(r * A0, 2))  -  pow(Ap0 / A0, 2) - Ap0 / (A0 * r) - Ap0 * alphap0 / (A0 * alpha0) + Xp0 * (Ap0 + 1./(2. * M_PI * r * A0)) / (A0 * X0)) * L
               + ( -alphap0 / alpha0 + Xp0 / X0 - 2. / r) * Fp - Lp / (4. * M_PI * r * A0_sq);


    double Lpp = 32. * ((M_PI * r *(A0 * Ap0 + alphap0 * A0_sq / alpha0) + 4. * pow(A0, 3)* M_PI * r * solitonic
               * (3. * pow(A0,3) * alphap0 + 9. * A0_sq * Ap0 * alpha0 - 2. * A0 * alphap0 * sigma_sq - 4. * alpha0 * Ap0 * sigma_sq) / (alpha0 * pow(sigma_sq, 2 )) ) * X0_sq
               + r * M_PI * (0.5 * Xp0 * A0_sq + solitonic * (6. * pow(A0, 6) * Xp0 / pow(sigma_sq, 2) - 4. * pow(A0_sq,2) * Xp0 / sigma_sq )  ) * X0 - M_PI * Ap0 * Ap0) * F
               + (16. * M_PI * Ap0 * Ap0 - chi_sq * X0_sq / alpha0_sq - 2. * pow(alphap0 / alpha0, 2) + 2. / (r * r) - 4. * alphap0 / (r * alpha0)
               +2.* (Xpp0  + 2. * alphap0 * Xp0 / alpha0 - Xp0 / r) / X0 - 4. * pow(Xp0 / X0, 2)) * L
               + 32. * M_PI * ( - Ap0 * A0 + 0.5 * A0_sq * r * X0_sq + r * X0_sq * solitonic * (6. * pow(A0_sq, 3) / pow(sigma_sq, 2) - 4. * pow(A0_sq, 2) / sigma_sq))  * Fp
               + 3. * (-alphap0 / alpha0 + Xp0 / X0) * Lp;

    //if (r < 0.1) cout<< setprecision(16) << "Results: [" << r << ", " << Fp << ", " << Fpp << ", " << Lp << ", " << Lpp << "]" << endl;

    return (PertState){Fp, Fpp, Lp, Lpp};
}

//construct perturbation solutions for r, gamma
void LinearPerturbation::rk4_solve(double chi_sq, long double gamma)
{
    pert.resize(n_gridpoints);
    pert[0] = (PertState){1, 0, 0, 0};

    //inter-level state values for RK4 evolution
    PertState s1, s2, s3, s4;

    //fill in perturbations on grid using RK4 evolution
    for (int j = 0; j < n_gridpoints - 1; j++)
    {
        double r = j * dr;

        s1 = pert_rhs(r, pert[j], chi_sq, gamma);
        s2 = pert_rhs(r + dr / 2., pert[j] + 0.5 * dr * s1, chi_sq, gamma);
        s3 = pert_rhs(r + dr / 2., pert[j] + 0.5 * dr * s2, chi_sq, gamma);
        s4 = pert_rhs(r + dr, pert[j] + dr * s3, chi_sq, gamma);

        pert[j + 1] = pert[j] + (dr / 6.) * (s1 + 2. * s2 + 2. * s3 + s4);

        if (isnan(pert[j].F) || isnan(pert[j].Fp) || isnan(pert[j].L) || isnan(pert[j].Lp))
        {
            //cerr << " \nWARNING: Perturbation values have become nan on step " << j << endl;
            blowup_point = j;
            return;
        }

        //if (j < 5) cout<< setprecision(16) << "Results: [" << (dr) * s3.F << ", " << (dr) * s3.Fp << ", " << (dr) * s3.L << ", " <<  (dr) * s3.Lp << "] " << dr << endl;

    }
    blowup_point = n_gridpoints - 1; //if blowup does not happen, set to outer boundary.
}

//returns number of zero crossings in L
int LinearPerturbation::count_zero_crossings()
{
    int zero_crossings = 0;
    int max_search_val =  round((n_gridpoints - 1) * cutoff_radius / R); // avoid counting spurious numerical zero crossings

    //check for zero crossings up to blowup point (neglecting one extra point as we sometimes get spurious zero crossings on explosion)
    for (int j = 1; j < max_search_val; j++)
    {
        if ((pert[j].L == 0.) || (pert[j].L > 0. && pert[j - 1].L < 0.) || (pert[j].L < 0. && pert[j - 1].L > 0.) )
            zero_crossings++;
    }

    //cout << "\n" << zero_crossings << " zero crossings found" << endl;

    return zero_crossings;
}

// given chi_sq, attempts to use interval bisection to find the gamma-value minimizing L at bdry.
long double LinearPerturbation::get_best_gamma(double chi_sq, bool quiet)
{
    long double sign = 1.;

    //relevant for chi_sq < 0: account for the fact that our cutoff position in the oscillating L profile may begin with L > 0 for small gamma, not the other way around
    rk4_solve(chi_sq, 0.0);
    if (pert[blowup_point].L < 0)
        sign = -1.;

    //lower bound on gamma
    long double lower_guess = gamma0;
    long double epsilon = bg->freq_epsilon; // just use the epsilon for the frequency; this seems to be OK generally

    rk4_solve(chi_sq, lower_guess);
    //double& L_inf = pert[blowup_point].L;

    int counter = 0;

    while (pert[blowup_point].L * sign < 0. && counter <= 10 ) //primitively looks for suitable lower bound if not immediately found
    {
        lower_guess -= 0.1 * lower_guess + 0.01;
        rk4_solve(chi_sq, lower_guess);
        counter++;

        if (counter > 10) cout << "WARNING: Failed to find lower gamma guess with gamma0 = " << gamma0 << " and chi_sq = " << chi_sq << endl;
    }

    long double upper_guess = abs( lower_guess * 2. + 0.1);

    rk4_solve(chi_sq, upper_guess);
    counter = 0;

    while (pert[blowup_point].L * sign > 0. && counter <= 10 ) //primitively looks for suitable lower bound if not immediately found
    {
        lower_guess += 0.1 * upper_guess + 0.01;
        rk4_solve(chi_sq, upper_guess);
        counter++;

        if (counter > 10) cout << "WARNING: Failed to find upper gamma guess with gamma0 = " << gamma0 << " and chi_sq = " << chi_sq << endl;
    }

    counter = 0;
    long double midpoint = 0.5 * (upper_guess + lower_guess);

    while ( (upper_guess - lower_guess) > epsilon && counter < 50)
    {
        //replace upper/lower bound with midpoint if it has greater/as many or fewer zero crossings than desired, so both ultimately converge on boundary value between eigen and eigen + 1 zero crossings as needed
        midpoint = 0.5 * (upper_guess + lower_guess);
        rk4_solve(chi_sq, midpoint);

        //cout << midpoint << endl;

        if (pert[blowup_point].L * sign < 0.)
            upper_guess = midpoint;
        else
            lower_guess = midpoint;

        counter++;

    }
    if (!quiet)
        cout << " \nBest gamma with chi_sq = " << chi_sq  << " is gamma =" << lower_guess << " with noether_pert = " << get_noether_perturbation() << " and " <<
        count_zero_crossings() << " zero crossings." << endl;

    solved_gamma = lower_guess;
    return lower_guess;
}

//use interval bisection to obtain chi_sq, obtaining the optimal gamma for each chi_sq value tested by the same.
double LinearPerturbation::get_chi_sq()
{
    double epsilon = chi_epsilon;

    double lower_guess = chi_sq0;
    double upper_guess = chi_sq0 + chi_range;

    get_best_gamma(lower_guess, 1);

    //for chi_sq < 0 this can happen if chi_sq is too small; try to find a guess with no zero crossings if this happens
    if (count_zero_crossings() > 0)
    {
        if (lower_guess >= 0.)
            cout << "WARNING: lower guess of chi_sq = "  <<  lower_guess << " too high!" << endl;
        else
        {

            cout << "WARNING: zero crossings found with negative lower guess chi_sq = "  <<  lower_guess << ", attempting to find good lower guess" << endl;
            double lower_search = lower_guess;

            int j = 1, k = 1;

            while (count_zero_crossings() > 0 && k < 4)
            {
                lower_search = lower_guess * (1. - j * pow(10., -k));
                j++;

                get_best_gamma(lower_search, 1);

                if (j >= 10)
                {
                    j = 0;
                    k++;
                }
            }

            if (count_zero_crossings() > 0)
                cout << "Search failed..." << endl;
            else
            {
                cout << "Search succeeded with chi_sq = " << lower_search << endl;
                lower_guess = lower_search;
            }

        }
    }

    get_best_gamma(upper_guess, 1);
    if (count_zero_crossings() == 0)
        cout << "WARNING: upper guess of chi_sq = " << upper_guess << " too low!" << endl;

    int counter = 0;

    long double midpoint = 0.5 * (upper_guess + lower_guess);

    while (upper_guess - lower_guess > epsilon && counter < 20 )
    {
        midpoint = 0.5 * (upper_guess + lower_guess);
        //cout << midpoint << endl;
        get_best_gamma(midpoint, 1);

        if (count_zero_crossings() > 0)
            upper_guess = midpoint;
        else
            lower_guess = midpoint;

        counter++;
    }

    cout << "Obtained chi_sq = " << lower_guess << " with noether_pert = " << get_noether_perturbation() << endl;

    solved_chi_sq = lower_guess;
    return lower_guess;

}


//alternative: try using secant method to get chi, searching for root of noether perturbation
double LinearPerturbation::get_chi_sq_newton()
{
    double epsilon = chi_epsilon;

    double chi_sq = chi_sq0;
    double chi_sq_prev = chi_sq0 * 0.99; //0.99 is ad hoc

    get_best_gamma(chi_sq, 1);
    double noether = get_noether_perturbation();

    get_best_gamma(chi_sq_prev, 1);
    double noether_prev = get_noether_perturbation();

    int counter = 0;

    while(abs(chi_sq - chi_sq_prev) > epsilon && counter < 50 && abs(noether) > 1e-11)
    {
        double chi_sq_new = chi_sq - noether * (chi_sq - chi_sq_prev) / (noether - noether_prev);

        get_best_gamma(chi_sq_new, 1);
        double noether_new = get_noether_perturbation();

        chi_sq_prev = chi_sq;
        noether_prev = noether;
        chi_sq = chi_sq_new;
        noether = noether_new;

        counter++;

        cout << "chi_sq = " << chi_sq << ", noether = " << noether << endl;
    }

    cout << "Obtained chi_sq = " << chi_sq << " with noether_pert = " << get_noether_perturbation() << endl;

    solved_chi_sq = chi_sq;
    return chi_sq;


}

double LinearPerturbation::get_noether_perturbation()
{
    if (blowup_point < n_gridpoints - 1) cout << "\n WARNING: perturbation solution nan'd; Noether perturbation likely to be inaccurate" << endl;

    int edge_point = n_gridpoints - 1;
    double R_edge = (R * edge_point) / (n_gridpoints - 1);

    FieldState& s = bg->state[edge_point];

    //both terms included as F generically grows exponentially (indeed 1st term usually dominates)
    noether_perturbation = 8. * M_PI * exp(s.phi) * s.A * s.eta * R_edge * R_edge * pert[edge_point].F / (omega)
                                - exp(s.phi) * R_edge * pert[edge_point].L / (omega * s.X );


    //cout << 8. * M_PI * exp(s.phi) * s.A * s.eta * R * R * pert[edge_point].F / omega << endl;
    return noether_perturbation;
}

void LinearPerturbation::write_pert(string filename)
{

    std::ofstream data_file{filename};

    // If we couldn't open the output file stream for writing
    if (!data_file)
    {
        // Print an error and exit
        std::cerr << "pert.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 0; j < n_gridpoints; j++)
    {
        data_file << std::setprecision (10) <<j * dr << "   " << pert[j].F << "    " << pert[j].Fp << "    " << pert[j].L << "    " << pert[j].Lp << "    " << "toby" << endl;
    }
}

void LinearPerturbation::write_chi_results()
{
     ofstream data_file{"chi.dat"};
     data_file << std::setprecision (10) <<  bg->A_central << "     " << solved_chi_sq << "     " << solved_gamma << "     "    << get_noether_perturbation() << "     " << bg->M << "     " << bg->omega << endl;
}

//cycles through perturbations as function of A0 in single-valued region
void LinearPerturbation::pert_cycle(double A0, double dA, int n_stars)
{
    bg->A_central = A0;
    ofstream data_file{"chi.dat"};

    double chi_prev, gamma_prev;

    //use if we want to switch behaviour when chi^2 switches sign
    bool chi_sq_positive = 1;
    if (chi_sq0 < 0.)
        chi_sq_positive = 0;

    for (int j = 0; j < n_stars; j++)
    {
        bg->solve(1);
        omega = bg->omega;
        A_central = bg-> A_central;

        if (j > 1)
        {
            gamma0 = 2.0 * solved_gamma - gamma_prev; //just linearly extrapolate the gamma guess, if this is a little off in either direction it should be OK

            if (chi_sq_positive)
                chi_sq0 = 0.95 * solved_chi_sq - 2.0 * abs(solved_chi_sq - chi_prev); //must be more careful here, as chi_sq0 must be < the true chi_sq: factors are a bit ad hoc
            else
                chi_sq0 = 2.0 * solved_chi_sq - chi_prev; //using Newton's method, want to be close, but no hard buffer needed
        }

        if (j > 0)
        {
            chi_prev = solved_chi_sq;
            gamma_prev = solved_gamma;
        }

        //use bisection if we expect chi_sq > 0 and Newton's method otherwise
        if (chi_sq_positive)
            get_chi_sq();
        else
            get_chi_sq_newton();

        data_file << std::setprecision (10) <<  bg->A_central << "     " << solved_chi_sq << "     " << solved_gamma << "     "    << get_noether_perturbation() << "     " << bg->M <<  "     " << bg->omega << endl;
        bg->A_central += dA;
        chi_epsilon = 0.0001 * solved_chi_sq + 10e-15;

        if (solved_chi_sq <= 0.)
            chi_sq_positive = 0;
        else
            chi_sq_positive = 1;
    }

}

//find and return value of first extremum in L
double LinearPerturbation::last_extremum()
{
    double ext_val = 0.;

    for (int j = n_gridpoints - 2; j > 0; j--)
    {
        if (abs(pert[j].L) > abs(pert[j - 1].L) &&  abs(pert[j].L) > abs(pert[j + 1].L))
        {
            ext_val = pert[j].L;
            break;
        }
    }

    return ext_val;
}

//testing RHS only (WARNING: BREAKS BACKGROUND BOSON STAR!!!)
PertState LinearPerturbation::test_rhs (double r, FieldState f, PertState s, double chi_sq, long double gamma)
{
    for (int j = 0; j < n_gridpoints; j++)
    {
        bg->state[j] = f;
    }

    PertState p = pert_rhs(r, s, chi_sq, gamma);

    cout<< "Results: [" << p.F << ", " << p.Fp << ", " << p.L << ", " << p.Lp << "]" << endl;
    return p;
}

//read params specific to computing linear perturbations
void LinearPerturbation::read_parameters(bool quiet)
{
    ifstream params{ "BSParams.par" };

    // Print an error and exit if file cannot open
    if (!params)
    {
        std::cerr << "Could not open BSParams.par\n";
        abort();
    }

    string current_line{};

    while (getline(params, current_line))
    {
        fill_parameter(current_line, "chi_sq0 = ", chi_sq0, quiet);
        fill_parameter(current_line, "gamma0 = ", gamma0, quiet);
        fill_parameter(current_line, "chi_range = ", chi_range, quiet);
        fill_parameter(current_line, "chi_epsilon = ", chi_epsilon, quiet);
        fill_parameter(current_line, "cutoff_radius = ", cutoff_radius, quiet);
    }
}
#endif /* LINEARPERTURBATION_CPP_ */
