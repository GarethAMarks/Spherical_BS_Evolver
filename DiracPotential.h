#ifndef DIRACPOTENTIAL_H_
#define DIRACPOTENTIAL_H_

#include <cmath>

/*
 Self-interaction of the two-spinor Dirac field (DiracStar, ComplexSpinorField), as a function of the scalar density
 S = sum_m psibar_m psi_m:
     L = sum_m -psibar_m (gamma^a nabla_a + mu) psi_m - U(S),
     U(S) = (lambda/2) S^2 + [solitonic] mu S ( (1 - S/S0)^2 - 1 ),     S0 = sigma^2 / (4 pi),
 so that with solitonic = 1 (and lambda = 0) the mass term plus self-interaction is the degenerate-vacuum potential
     W(S) = mu S + U(S) = mu S (1 - S/S0)^2,
 vanishing at S = 0 and at the second vacuum S = S0. At the centre of a star S = A_central^2 / (2 pi), so the vacuum is at
 A_central = sigma / sqrt(2), as for the solitonic boson stars (V = mu^2 A^2 (1 - 2 A^2/sigma^2)^2). In terms of
 Herdeiro & Radu (2020), mu S - lambda_HR S^2 + nu S^3 with beta = nu mu / lambda_HR^2 = 1/4 (degenerate case).
 The spinors obey the Dirac equation with mu_eff = mu + U'(S); with solitonic = 1, mu_eff = mu (1 - S/S0)(1 - 3 S/S0)
 (+ lambda S), which vanishes in the second vacuum. The stress tensor gains g_ab P(S), P = S U'(S) - U(S).
*/
struct DiracPotential
{
    double mu = 1.;
    double lambda = 0.;
    bool solitonic = false;
    double sigma = 1.;

    double S0() const { return sigma * sigma / (4. * M_PI); }

    double U(double S) const
    {
        double u = 0.5 * lambda * S * S;
        if (solitonic)
        {
            double y = S / S0();
            u += mu * S * (y * y - 2. * y);
        }
        return u;
    }

    double dU(double S) const
    {
        double du = lambda * S;
        if (solitonic)
        {
            double y = S / S0();
            du += mu * (3. * y * y - 4. * y);
        }
        return du;
    }

    double mu_eff(double S) const { return mu + dU(S); }
    double P(double S) const { return S * dU(S) - U(S); }
    bool free_field() const { return lambda == 0. && !solitonic; }
};

#endif /* DIRACPOTENTIAL_H_ */
