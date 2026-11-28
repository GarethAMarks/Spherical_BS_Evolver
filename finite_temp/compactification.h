#ifndef COMPACTIFICATION_H
#define COMPACTIFICATION_H

#include <cmath>
#include <stdexcept>
#include <iostream>

/**
 * @class Compactification
 * @brief Abstract interface for coordinate compactification schemes.
 *
 * Defines the transformation between physical radius r ∈ [0, ∞) and 
 * compactified radius r_bar ∈ [0, 1).
 *
 * To implement a custom scheme, derive from this class and override all methods
 * with closed-form expressions. Ensure derivatives are analytically correct.
 *
 * Example: Standard scheme r_bar = r / (1 + r) ⟺ r = r_bar / (1 - r_bar)
 */
class Compactification {
public:
    virtual ~Compactification() = default;

    /**
     * @brief Forward transformation: physical radius → compactified radius
     * @param r Physical radius (r ≥ 0)
     * @return Compactified radius r_bar ∈ [0, 1)
     */
    virtual double r_of_rbar(double rbar) const = 0;

    /**
     * @brief Inverse transformation: compactified radius → physical radius
     * @param rbar Compactified radius (0 ≤ r_bar < 1)
     * @return Physical radius r ≥ 0
     */
    virtual double rbar_of_r(double r) const = 0;

    /**
     * @brief First derivative: dr_bar / dr
     * @param r Physical radius
     * @return Derivative value
     *
     * This is the Jacobian of the inverse transformation.
     * For standard scheme: dr_bar/dr = (1 - r_bar)^2 = 1 / (1 + r)^2
     */
    virtual double drbar_dr(double r) const = 0;

    /**
     * @brief Second derivative: d²r_bar / dr²
     * @param r Physical radius
     * @return Derivative value
     *
     * For standard scheme: d²r_bar/dr² = -2(1 - r_bar)^3 = -2 / (1 + r)^3
     */
    virtual double d2rbar_dr2(double r) const = 0;

    /**
     * @brief Jacobian factor for first derivatives in r-space.
     * Encodes the chain rule: df/dr = (dr_bar/dr)^(-1) * df/dr_bar
     *
     * @param rbar Compactified radius
     * @return Factor to multiply df/dr_bar by to get df/dr
     *
     * For standard scheme: (1 - r_bar)^2
     */
    virtual double d1_r_factor(double rbar) const = 0;

    /**
     * @brief Jacobian factor for second-derivative term in r-space.
     * Encodes: d²f/dr² = (dr_bar/dr)^(-2) * d²f/dr_bar² + (d²r_bar/dr²) * (dr_bar/dr)^(-3) * df/dr_bar
     *
     * @param rbar Compactified radius
     * @return Factor to multiply d²f/dr_bar² by
     *
     * For standard scheme: (1 - r_bar)^4
     */
    virtual double d2_r_factor_first(double rbar) const = 0;

    /**
     * @brief Jacobian factor for first-derivative term in second derivative.
     * See d2_r_factor_first() for context.
     *
     * @param rbar Compactified radius
     * @return Factor to multiply df/dr_bar by (includes d²r_bar/dr² term)
     *
     * For standard scheme: -2(1 - r_bar)^3
     */
    virtual double d2_r_factor_second(double rbar) const = 0;

    /**
     * @brief Validate that the inverse transformation is consistent.
     * Tests: r_bar = r_bar_of_r(r_of_rbar(r_bar)) for a sample of points.
     *
     * @param tolerance Maximum allowed error
     * @return true if all tested points pass; false if any exceed tolerance
     *
     * Issues a warning to stdout if validation fails.
     */
    virtual bool validate_inverse(double tolerance = 1e-12) const = 0;
};

/**
 * @class StandardCompactification
 * @brief Standard boson star compactification: r_bar = r / (1 + r)
 *
 * Maps [0, ∞) ↔ [0, 1).
 * r_bar = 0 corresponds to the origin (r = 0).
 * r_bar → 1 corresponds to spatial infinity (r → ∞).
 *
 * Derivatives:
 * - dr_bar/dr = 1 / (1 + r)² = (1 - r_bar)²
 * - d²r_bar/dr² = -2 / (1 + r)³ = -2(1 - r_bar)³
 */
class StandardCompactification : public Compactification {
public:
    // optional length scale L: r = L r_bar / (1 - r_bar) (L = 1: the standard map). Larger L moves grid points outward.
    explicit StandardCompactification(double L = 1.0) : L_(L > 0.0 ? L : 1.0) {}
    double scale() const { return L_; }

    double r_of_rbar(double rbar) const override {
        const double one_minus = std::max(1e-12, 1.0 - rbar);
        return L_ * rbar / one_minus;
    }

    double rbar_of_r(double r) const override {
        return r / (L_ + r);
    }

    double drbar_dr(double r) const override {
        const double denom = L_ + r;
        return L_ / (denom * denom);
    }

    double d2rbar_dr2(double r) const override {
        const double denom = L_ + r;
        return -2.0 * L_ / (denom * denom * denom);
    }

    double d1_r_factor(double rbar) const override {
        const double one_minus = 1.0 - rbar;
        return one_minus * one_minus / L_;
    }

    double d2_r_factor_first(double rbar) const override {
        const double one_minus = 1.0 - rbar;
        return one_minus * one_minus * one_minus * one_minus / (L_ * L_);
    }

    double d2_r_factor_second(double rbar) const override {
        const double one_minus = 1.0 - rbar;
        return -2.0 * one_minus * one_minus * one_minus / (L_ * L_);
    }

    bool validate_inverse(double tolerance) const override {
        // Test a sample of points across the domain
        const int n_test = 20;
        bool all_pass = true;

        for (int i = 0; i < n_test; ++i) {
            double rbar_test = i / static_cast<double>(n_test - 1) * 0.999;  // Avoid exact 1.0
            double r = r_of_rbar(rbar_test);
            double rbar_recovered = rbar_of_r(r);
            double error = std::abs(rbar_test - rbar_recovered);

            if (error > tolerance) {
                std::cerr << "StandardCompactification::validate_inverse() WARNING:\n"
                          << "  At r_bar = " << rbar_test << ", r = " << r << ":\n"
                          << "  Recovered r_bar = " << rbar_recovered << "\n"
                          << "  Error = " << error << " (tolerance = " << tolerance << ")\n";
                all_pass = false;
            }
        }

        if (all_pass) {
            std::cout << "StandardCompactification::validate_inverse() passed "
                      << "(" << n_test << " test points, tolerance = " << tolerance << ")\n";
        }
        return all_pass;
    }

private:
    double L_ = 1.0;
};

#endif // COMPACTIFICATION_H
