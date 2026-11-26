#ifndef FINITE_TEMP_RELAXER_SEED_SHOOTER_H_
#define FINITE_TEMP_RELAXER_SEED_SHOOTER_H_

#include <vector>
#include <string>

namespace seedbs
{

struct SeedState
{
    double A {};
    double X {};
    double phi {};
    double eta {};

    SeedState() = default;
    SeedState(double A_, double X_, double phi_, double eta_) : A(A_), X(X_), phi(phi_), eta(eta_) {}
};

SeedState operator+(const SeedState& s1, const SeedState& s2);
SeedState operator-(const SeedState& s1, const SeedState& s2);
SeedState operator*(double c, const SeedState& s);

class BosonStarSeedShooter
{
public:
    double D = 4.0;
    double mu = 1.0;
    double lambda = 0.0;
    bool solitonic = false;
    double G = 1.0;
    double sigma_scale = 1.0;
    int eigen = 0;
    double alpha_central = 1.0;
    double A_central = 0.08;
    double frequency_guess = 0.85;
    double r_match_fac = 0.75;
    double R = 40.0;
    int n_gridpoints = 2001;
    double freq_epsilon = 1e-10;

    std::vector<SeedState> state;
    std::vector<double> radius_array;
    int blowup_point = 0;
    long double omega = 0.0;

    bool solve(bool quiet = false);
     void write_field(const std::string& filename = "seed_bs.dat") const;

private:
    double V(const double A) const;
    double dV(const double A) const;
    double ddV(const double A) const;

    SeedState state_RHS(const double radius, long double frequency, SeedState s, bool asymptotic_region) const;
    SeedState state_expansion(const double radius, long double frequency) const;

    void rk4_solve(long double freq);
    int count_zero_crossings() const;
    long double find_frequency(bool quiet = false);

    int find_last_minimum() const;
    bool fill_asymptotic(bool quiet = false);

    void rescale_lapse(double phi_shift);
    double m(int j) const;
};

} // namespace seedbs

#endif // FINITE_TEMP_RELAXER_SEED_SHOOTER_H_
