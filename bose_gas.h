#ifndef BOSE_GAS_H
#define BOSE_GAS_H
// Free non-relativistic gas for the finite-temperature EoS, P_gas = (m/2pi)^{3/2} T^{5/2} F(x), x = (mu - m)/T.
//   Boltzmann statistics:     F(x) = e^x
//   Bose-Einstein statistics: F(x) = Li_{5/2}(e^x) for x < 0, dF/dx = Li_{3/2}(e^x).
// A free Bose gas does not exist for mu > m (it condenses), so for x >= 0 the Bose branch is continued C^1 by
// F = zeta(5/2) e^{kappa0 x}, kappa0 = zeta(3/2)/zeta(5/2), which only matters where the interpolation hands over to phonons.
// Thermodynamics: n_gas = P_gas kappa / T, s_gas = P_gas (5/2 - x kappa) / T, with kappa = F'/F.
#include <cmath>

namespace bose_gas
{
// zeta(5/2 - k), k = 0..15
static const double kZeta52[16] = {1.341487257250917, 2.612375348685488, -1.4603545088095868, -0.20788622497735457,
    -0.025485201889833036, 0.008516928777850331, 0.004441011335479432, -0.0030916692472158338, -0.0026714580198992244,
    0.0027467679395368687, 0.00326903957260022, -0.00441603287300489, -0.006672172296466641, 0.011146122473942813,
    0.02039697871594279, -0.04057496748119458};
static const double kGammaM32 = 2.363271801207355;  // Gamma(-3/2)
static const double kGammaM12 = -3.544907701811032; // Gamma(-1/2)

// ln F(x) and kappa(x) = F'(x)/F(x).  stats = 0: Boltzmann, 1: Bose-Einstein.
inline void lnF_kappa(double x, int stats, double& lnF, double& kappa)
{
    if (stats == 0) { lnF = x; kappa = 1.; return; }
    if (x >= 0.)
    {
        const double k0 = kZeta52[1] / kZeta52[0];
        lnF = std::log(kZeta52[0]) + k0 * x; kappa = k0; return;
    }
    if (x < -0.7)
    {
        // series in z = e^x <= 0.497: Li_s(z) = z sum_k z^{k-1} / k^s
        const double z = std::exp(x);
        double S5 = 0., S3 = 0., zk = 1.;
        for (int k = 1; k <= 60; k++)
        {
            const double kk = (double)k;
            S5 += zk / (kk * kk * std::sqrt(kk));
            S3 += zk / (kk * std::sqrt(kk));
            zk *= z;
            if (zk < 1e-18) break;
        }
        lnF = x + std::log(S5); kappa = S3 / S5; return;
    }
    // -0.7 <= x < 0: Li_s(e^x) = Gamma(1-s) (-x)^{s-1} + sum_k zeta(s-k) x^k / k!
    const double mx = -x;
    double L5 = kGammaM32 * mx * std::sqrt(mx), L3 = kGammaM12 * std::sqrt(mx);
    double term = 1.;
    for (int k = 0; k < 15; k++)
    {
        L5 += kZeta52[k] * term;
        L3 += kZeta52[k + 1] * term;
        term *= x / (double)(k + 1);
    }
    lnF = std::log(L5); kappa = L3 / L5;
}
} // namespace bose_gas
#endif
