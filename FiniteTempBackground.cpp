// Finite-temperature background read-in for SBSE: runs or loads the relaxation solver (finite_temp/), converts its
// compactified profile (rbar = r / (1 + r)) to SBSE's uniform physical grid, and fills the BosonStar arrays that
// BSSNSlice::read_BS_data uses, plus the background thermodynamics in ft_state. See FiniteTempBackground.h.

#include "BosonStar.h"
#include "FiniteTempBackground.h"
#include "finite_temp/FiniteTempRelaxer.h"
#include "spline.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

namespace
{

// Interpolates the relaxer's compactified-grid profile onto physical radii with cubic splines in r through the
// knots r_i = rbar_i / (1 - rbar_i). Splining in r (not rbar) lets the origin parity be imposed exactly:
// sigma, alpha, Xs are even (zero slope at r = 0) and eta = sigma' / Xs is odd (zero curvature at r = 0).
class ProfileInterpolator
{
public:
    struct Fields { double sigma, eta, alpha, Xs; };

    // r_max: largest physical radius that will be requested. Returns false if the relaxer grid does not reach it.
    bool build(const FiniteTempRelaxer& relaxer, double r_max)
    {
        const RelaxerState& st = relaxer.state();
        const size_t n = st.rbar.size();
        std::vector<double> sig, eta, alp, xs;

        // Stop a few knots past r_max (the rest cluster towards rbar = 1 and are never needed); never use rbar = 1.
        int beyond = 0;
        for (size_t i = 0; i < n && beyond < 4; ++i)
        {
            if (st.rbar[i] >= 1. - 1e-14)
                break;
            const double r = relaxer.areal_radius(st.rbar[i]);
            r_knots.push_back(r);
            sig.push_back(st.sigma[i]);
            eta.push_back(st.eta[i]);
            alp.push_back(st.alpha[i]);
            xs.push_back(st.Xs[i]);
            if (r > r_max)
                beyond++;
        }

        if (r_knots.size() < 4 || r_knots.back() < r_max)
        {
            cerr << "ERROR: finite-T relaxer grid only reaches r = " << (r_knots.empty() ? 0. : r_knots.back())
                 << " but SBSE needs r = " << r_max << "; raise rbar_max towards 1 in the relaxer" << endl;
            return false;
        }

        s_sigma.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_alpha.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_Xs.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_eta.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0);

        s_sigma.set_points(r_knots, sig, tk::spline::cspline);
        s_eta.set_points(r_knots, eta, tk::spline::cspline);
        s_alpha.set_points(r_knots, alp, tk::spline::cspline);
        s_Xs.set_points(r_knots, xs, tk::spline::cspline);
        return true;
    }

    Fields operator()(double r) const
    {
        return Fields{s_sigma(r), s_eta(r), s_alpha(r), s_Xs(r)};
    }

    // knot spacing of the relaxer grid around radius r (for reporting the effective input resolution)
    double knot_spacing(double r) const
    {
        for (size_t i = 1; i < r_knots.size(); ++i)
            if (r_knots[i] >= r)
                return r_knots[i] - r_knots[i - 1];
        return r_knots.size() > 1 ? r_knots.back() - r_knots[r_knots.size() - 2] : 0.;
    }

private:
    std::vector<double> r_knots;
    tk::spline s_sigma, s_eta, s_alpha, s_Xs;
};

template <typename T>
void warn_if_differs(const char* name, T relaxer_val, T sbse_val)
{
    if (relaxer_val != sbse_val)
        cout << "WARNING: finite-T relaxer parameter " << name << " = " << relaxer_val
             << " overridden by BSParams.par value " << sbse_val << endl;
}

} // namespace

void BosonStar::read_finite_temp()
{
    if (D != 4.)
    {
        cerr << "ERROR: finite_temp = 1 requires D = 4 (the relaxer equations are four-dimensional)" << endl;
        exit(1);
    }
    if (isotropic)
    {
        cerr << "ERROR: finite_temp = 1 requires isotropic = 0 (backgrounds are in polar-areal coordinates)" << endl;
        exit(1);
    }

    FiniteTempRelaxer relaxer;
    relaxer.read_parameters(ft_par_file, true);

    // The scalar potential must be the one SBSE evolves, so BSParams.par wins. Both codes use
    // V = mu^2 A^2 + lambda A^4 (times (1 - 2 A^2/sigma^2)^2 on the mass term if solitonic), so the values map 1:1.
    RelaxerParams& rp = relaxer.params();
    warn_if_differs("mu", rp.mu, mu);
    warn_if_differs("lambda", rp.lambda, lambda);
    warn_if_differs("solitonic", rp.solitonic, solitonic);
    if (solitonic) warn_if_differs("sigma_scale", rp.sigma_scale, sigma);
    warn_if_differs("G", rp.G, G);
    ft_cph_epsilon = rp.cph_epsilon;
    ft_q = rp.pressure_interp_power;
    ft_suppression = rp.use_suppression_factor;
    rp.mu = mu;
    rp.lambda = lambda;
    rp.solitonic = solitonic;
    rp.sigma_scale = sigma;
    rp.G = G;
    if (G != 1.)
        cout << "WARNING: the finite-T relaxer equations assume G = 1" << endl;

    if (ft_relax)
    {
        // A single thermalised model is wanted, not the omega family the standalone relaxer can trace afterwards.
        if (!rp.first_thermalized_only)
            cout << "NOTE: setting first_thermalized_only = 1 for the in-SBSE relaxer run" << endl;
        rp.first_thermalized_only = true;
        cout << "\nRunning finite-temperature relaxer (" << ft_par_file << ")..." << endl;
        if (!relaxer.run(false))
        {
            cerr << "ERROR: finite-temperature relaxation failed" << endl;
            exit(1);
        }
        //polish the final model (a no-op after a converged temperature stage; relaxes the bare seed if T_target = 0).
        //With force_continue the relaxer can return unconverged models, so refuse anything clearly not a solution.
        const bool polished = relaxer.relax_current_state();
        cout.unsetf(std::ios::floatfield); //the relaxer leaves cout in scientific format
        cout << std::setprecision(6);
        if (!polished)
        {
            const double res = relaxer.last_residual_norm();
            if (!(res <= ft_max_residual))
            {
                cerr << "ERROR: finite-T relaxer did not converge (residual " << res << " > ft_max_residual = "
                     << ft_max_residual << ")" << endl;
                exit(1);
            }
            cout << "WARNING: finite-T model residual " << res << " is above newton_tolerance but within ft_max_residual" << endl;
        }
        relaxer.write_profile(ft_profile_file);
    }
    else if (!relaxer.load_profile(ft_profile_file))
    {
        cerr << "ERROR: could not read finite-T profile " << ft_profile_file << endl;
        exit(1);
    }
    cout.unsetf(std::ios::floatfield);
    cout << std::setprecision(6);

    // Physical uniform grid, as for the T = 0 solver (vertex-centred; read_BS_data handles cell-centred offsets).
    const double dr = R / (n_gridpoints - 1.);
    ProfileInterpolator interp;
    if (!interp.build(relaxer, R + 4. * dr))
        exit(1);

    state.resize(n_gridpoints);
    radius_array.resize(n_gridpoints);
    ft_state.resize(n_gridpoints);

    omega = relaxer.omega_infinity();
    omega_pre_rescale = omega; // the relaxer already has alpha(infinity) = 1
    T_inf = relaxer.temperature();

    for (int k = 0; k < n_gridpoints; k++)
    {
        const double r = k * dr;
        const ProfileInterpolator::Fields f = interp(r);
        radius_array[k] = r;
        state[k].A = f.sigma;
        state[k].X = f.Xs;
        state[k].phi = std::log(std::max(1e-300, f.alpha));
        state[k].eta = f.eta;

        const FiniteTempRelaxer::PointThermo th = relaxer.thermo_at(f.sigma, f.eta, f.alpha, (double)omega);
        FiniteTempPoint& p = ft_state[k];
        p.T = th.T; p.mu = th.mu; p.P = th.P; p.s = th.s; p.n = th.n;
        p.Upsilon = th.Upsilon; p.F = th.F; p.w = th.w; p.rho = th.rho;
    }
    state[0].eta = 0.;
    state[0].X = 1.;
    blowup_point = n_gridpoints;

    A_central = state[0].A;
    alpha_central = std::exp(state[0].phi);

    //global quantities on the SBSE grid, as for the T = 0 solver: M = m(R), r_99 from the mass function, and the
    //total (superfluid + normal) charge N = int 4 pi r^2 X (Upsilon mu + n) dr (draft sign, > 0)
    M = m(n_gridpoints - 1);
    r_99 = R;
    for (int k = 0; k < n_gridpoints; k++)
        if (m(k) >= 0.99 * M) { r_99 = radius_array[k]; break; }
    noether_charge = 0.;
    for (int k = 1; k < n_gridpoints; k++)
    {
        auto f = [&](int j) { return 4. * M_PI * radius_array[j] * radius_array[j] * state[j].X
                                     * (ft_state[j].Upsilon * ft_state[j].mu + ft_state[j].n); };
        noether_charge += 0.5 * (f(k) + f(k - 1)) * dr;
    }
    binding_energy = M - mu * noether_charge;
    compactness = M / r_99;

    // Consistency check of the conversion: Hamiltonian constraint in polar-areal gauge, m' = 4 pi r^2 rho with
    // m = r (1 - 1/X^2) / 2, using the interpolated X and the EoS rho evaluated on the SBSE grid.
    double max_res = 0., max_src = 0.;
    for (int k = 2; k < n_gridpoints - 2; k++)
    {
        const double r = radius_array[k];
        if (r > r_99 * 1.5) break;
        auto m_of = [&](int j) { return 0.5 * radius_array[j] * (1. - 1. / (state[j].X * state[j].X)); };
        const double dm = (-m_of(k + 2) + 8. * m_of(k + 1) - 8. * m_of(k - 1) + m_of(k - 2)) / (12. * dr);
        const double src = 4. * M_PI * r * r * ft_state[k].rho;
        max_res = std::max(max_res, std::abs(dm - src));
        max_src = std::max(max_src, std::abs(src));
    }

    cout << "\nRead finite-temperature background: T_inf = " << T_inf << ", omega = " << (double)omega
         << ", sigma0 = " << A_central << ", alpha0 = " << alpha_central << "\n  M(R) = " << M << " (relaxer M_inf = "
         << relaxer.asymptotic_mass() << "), N(R) = " << noether_charge << " (relaxer " << relaxer.particle_number() << ")"
         << ", r_99 = " << r_99 << "; relaxer knot spacing at r_99 = " << interp.knot_spacing(r_99)
         << ", at R = " << interp.knot_spacing(R) << " (SBSE dr = " << dr << ")"
         << "\n  Hamiltonian-constraint check on SBSE grid (r < 1.5 r_99): max |m' - 4 pi r^2 rho| = " << max_res
         << " (max |4 pi r^2 rho| = " << max_src << ")" << endl;

    write_finite_temp();
}

void BosonStar::write_finite_temp(std::string filename)
{
    std::ofstream out{filename};
    if (!out)
    {
        cerr << filename << " could not be opened for writing" << endl;
        return;
    }
    out << "# T_inf = " << T_inf << "  omega = " << (double)omega << "  (draft sign convention, phi ~ e^{-i omega t})\n";
    out << "# r A X alpha eta T mu P s n Upsilon F w rho N_S N_N Sigma\n";
    out << std::setprecision(12);
    for (int k = 0; k < (int)ft_state.size(); k++)
    {
        const FiniteTempPoint& p = ft_state[k];
        const FiniteTempConserved c = ft_static_conserved(p);
        out << radius_array[k] << ' ' << state[k].A << ' ' << state[k].X << ' ' << std::exp(state[k].phi) << ' '
            << state[k].eta << ' ' << p.T << ' ' << p.mu << ' ' << p.P << ' ' << p.s << ' ' << p.n << ' '
            << p.Upsilon << ' ' << p.F << ' ' << p.w << ' ' << p.rho << ' ' << c.N_S << ' ' << c.N_N << ' '
            << c.Sigma << '\n';
    }
}
