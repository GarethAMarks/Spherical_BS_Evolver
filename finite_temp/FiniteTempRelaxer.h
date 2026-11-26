#ifndef FINITE_TEMP_RELAXER_SOLVER_H_
#define FINITE_TEMP_RELAXER_SOLVER_H_

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "compactification.h"

struct RelaxerParams
{
    int n_gridpoints = 257;
    double rbar_min = 0.0;
    double rbar_max = 1.0;

    double G = 1.0;
    double mu = 1.0;
    double lambda = 0.0;
    bool solitonic = false;
    double sigma_scale = 1.0;

    double sigma_central = 0.08;
    double alpha_central = 1.0;
    double mu_chem_infty = 0.85;

    bool use_shooting_seed = true;
    double seed_R = 40.0;
    int seed_n_gridpoints = 2001;
    int seed_eigen = 0;
    double seed_frequency_guess = 0.85;
    double seed_freq_epsilon = 1e-10;
    double seed_r_match_fac = 0.75;

    double temperature_target = 0.0;
    double temperature_step = 0.0;
    bool force_continue = false;
    bool first_thermalized_only = false;
    bool use_suppression_factor = true;
    int decay_exponent = 10;
    double cph_epsilon = 1e-8;
    double pressure_interp_power = 1.0;

    double omega_target = 0.885;
    double omega_step = 0.001;
    double charge_step = 1e-4;
    double sigma_jump_on_fail = 5.0;
    bool increment_omega = false;
    int num_models = 100;

    std::vector<double> hi_sigma_vals;
    std::vector<double> lo_sigma_vals;
    double resolution_factor = 1.0;

    int max_newton_iterations = 16;
    double newton_tolerance = 1e-10;
    double jacobian_fd_step = 1e-6;
    double line_search_shrink = 0.5;
    double line_search_min_step = 1e-4;

    std::string output_profile = "profile.dat";
    bool save_model_profiles = false; // write model_NNNN.dat for every converged family member
};

struct RelaxerState
{
    std::vector<double> rbar;
    std::vector<double> sigma;
    std::vector<double> eta;
    std::vector<double> alpha;
    std::vector<double> Xs;
    std::vector<double> omega;
    std::vector<double> rho;
};

class FiniteTempRelaxer
{
public:
    FiniteTempRelaxer(std::shared_ptr<Compactification> comp = nullptr)
        : compactification_(comp ? comp : std::make_shared<StandardCompactification>()) {}

    struct ThermoConjugates
    {
        double s = 0.0;
        double n = 0.0;
        double F = 0.0;
        double Upsilon = 0.0;
    };

    void read_parameters(bool quiet = false);
    void read_parameters(const std::string& filename, bool quiet = false);
    bool run(bool quiet = false);
    void write_profile(const std::string& filename = "finite_temp_profile.dat") const;
    void refine_grid(double res_factor);

    // ---- Interface used by SBSE (Spherical_BS_Evolver) to read relaxed backgrounds ----

    // Load a profile previously written by write_profile(). The rbar grid is taken from the file and the
    // omitted outer point (rbar = rbar_max, i.e. r = infinity) is re-appended with its boundary values.
    // The global (Tolman) temperature T_inf = T_local * alpha is recovered from the file. The EoS and
    // potential parameters are NOT stored in the file; they come from params() as currently set.
    bool load_profile(const std::string& filename, bool quiet = false);

    // One Newton relaxation of the current state at the current temperature (boundary values re-imposed first).
    // Returns true if the residual reached newton_tolerance. Used by SBSE to polish the final model, which also
    // relaxes a bare shooting seed when run() stopped at T = 0 without a Newton stage.
    bool relax_current_state();
    double last_residual_norm() const { return last_residual_norm_; }

    RelaxerParams& params() { return p_; }
    const RelaxerParams& params() const { return p_; }
    const RelaxerState& state() const { return state_; }
    double temperature() const { return current_temperature_; } // T_inf (local T = T_inf / alpha)
    double omega_infinity() const { return state_.omega.empty() ? p_.mu_chem_infty : state_.omega.back(); }

    // Thermodynamics of a static background point: local T = T_inf / alpha, mu = omega / alpha,
    // X = eta^2 (= g^rr (d_r sigma)^2), Xi = -mu^2 (draft sign convention, phi ~ sigma e^{-i omega t}).
    struct PointThermo
    {
        double T = 0.0;       // local temperature
        double mu = 0.0;      // local chemical potential (> 0 in the draft convention)
        double X = 0.0;
        double Xi = 0.0;
        double P = 0.0;       // total pressure P = L + P_thermal
        double s = 0.0;       // entropy density
        double n = 0.0;       // normal number density
        double F = 0.0;
        double Upsilon = 0.0;
        double w = 0.0;       // normal enthalpy density T s + mu n
        double rho = 0.0;     // energy density n^a n^b T_ab = mu^2 Upsilon + w - P (static, v = 0)
    };
    PointThermo thermo_at(double sigma, double eta, double alpha, double omega) const;

    // EoS P(T, mu, sigma^2, X, Xi) and its conjugates. Note: with use_suppression_factor the thermal
    // part depends on the central amplitude state().sigma.front(), so a profile must be loaded first.
    double pressure_Pthermal(double T, double mu_chem, double sigma_sq, double X, double Xi) const;
    ThermoConjugates conjugates_from_pressure(double T, double mu_chem, double sigma_sq, double X, double Xi) const;

    double areal_radius(double rbar) const;
    double asymptotic_mass() const;
    double radius_at_mass_fraction(double mass_fraction) const;
    double particle_number() const;

private:
    std::shared_ptr<Compactification> compactification_;
    struct Residuals
    {
        std::vector<double> sigma;
        std::vector<double> eta;
        std::vector<double> alpha;
        std::vector<double> Xs;
        std::vector<double> omega;
        std::vector<double> rho;
    };

    RelaxerParams p_{};
    RelaxerState state_{};
    double current_temperature_ = 0.0;
    int last_newton_iterations_ = 0;
    double last_residual_norm_ = 0.0;
    double target_particle_number_ = 0.0;
    double last_saved_sigma0_ = 0.0;
    double last_sigma_step_ = 0.0;
    std::set<double> crossed_hi_sigma_thresholds_;
    std::set<double> crossed_lo_sigma_thresholds_;

    void build_grid();
    void initialize_seed();
    void enforce_boundary_values();
    void append_model_summary(const std::string& filename) const;

    double metric_mass(double r, double Xs) const;
    double particle_number_from_rho() const;
    size_t active_field_count() const;

    double potential(double sigma) const;
    double potential_of_sigma_sq(double sigma_sq) const;
    double d_potential_d_sigma2(double sigma) const;
    double ddV(double sigma_sq) const;
    double cph_sq(double sigma_sq, double Xi) const;
    double scalar_pressure_term(double sigma_sq, double X, double Xi) const;
    double thermal_pressure_term(double T, double mu_chem, double sigma_sq, double X, double Xi) const;
    double compute_suppression_factor(double sigma_sq) const;

    double d1_rbar(const std::vector<double>& u, int i, double h) const;
    double d2_rbar(const std::vector<double>& u, int i, double h) const;
    double d1_r_from_rbar(double u_b, double rbar) const;
    double d2_r_from_rbar(double u_b, double u_bb, double rbar) const;

    double local_temperature(double alpha) const;

    std::vector<double> pack_state() const;
    void unpack_state(const std::vector<double>& x);

    Residuals compute_residuals() const;
    std::vector<double> pack_residual_vector(const Residuals& R) const;

    double residual_norm(const Residuals& R) const;
    double vector_norm_inf(const std::vector<double>& v) const;

    void build_jacobian_fd(const std::vector<double>& x_base,
                           const std::vector<double>& r_base,
                           std::vector<double>& J) const;
    bool solve_linear_system(std::vector<double> A,
                             std::vector<double> b,
                             std::vector<double>& x) const;
    bool newton_solve_stage();

    std::vector<double> stage_targets(bool charge_mode, double start_value) const;
    bool shoot_for_sigma(double sigma0);

    // Diagnostic methods
    void print_field_diagnostics(const std::string& label) const;
    void print_residual_components(int index) const;
};

#endif // FINITE_TEMP_RELAXER_SOLVER_H_