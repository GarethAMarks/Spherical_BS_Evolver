"""
Self-interactions of the Dirac field in the 3+1 evolution (ComplexSpinorField), signature (-,+,+,+):

    L = sum_m -psibar_m (gamma^a nabla_a + mu) psi_m - U(S),   S = sum_m psibar_m psi_m,
    U(S) = (lambda/2) S^2 + kappa S^3     (quartic: kappa = 0; degenerate vacuum mu S + U = mu S (1 - S/S0)^2:
                                            lambda = -4 mu / S0, kappa = mu / S0^2).

Each spinor obeys (gamma^a nabla_a + mu_eff) psi_m = 0 with mu_eff = mu + U'(S). In the variables of ComplexSpinorField,
S = (|u|^2 - |v|^2) / (2 pi), so mu_eff is an algebraic function of the fields at each point. The mass enters the
derivation of check_dirac_3p1.py only algebraically (no derivatives of mu appear anywhere), so the evolution equations
are those of the free field with mu -> mu_eff:
    d_t u = beta u' + (alpha K/2 - i alpha mu_eff) u + i alpha [ ... ],   d_t v = beta v' + (alpha K/2 + i alpha mu_eff) v - i alpha [ ... ],
and the kinetic stress tensor obeys the same identities with mu -> mu_eff (in particular rho_kin = mu_eff S + S_rr^ + 2 S_thth^).
The full stress tensor is T_ab = T^kin_ab + g_ab P(S), P = S U' - U (see check_dirac_star.py), so
    rho     = mu S + U(S) + S_rr^ + 2 S_thth^       (= rho_kin - P),
    S_rr^   = Re( vbar u' - u vbar' ) / (2 pi sqrt(g_rr)) + P(S),
    S_thth^ = -Re( u vbar ) / (2 pi r sqrt(g_T)) + P(S),
    j_r^    = Im( ubar u' + vbar v' ) / (2 pi sqrt(g_rr))    (unchanged),
with S_rr^, S_thth^ on the right of the rho line the kinetic (lambda-independent) parts.

This script checks these sources against the evolution equations through local energy and momentum conservation,
    d_t rho + (r^2 j_r)' / r^2 = 0,     d_t j_r + S_rr' + 2 (S_rr - S_thth) / r = 0,
identically in flat space (alpha = 1, beta = 0, K = 0, g = flat), for arbitrary complex u(r), v(r). Conservation fixes the
normalisation of the g_ab P term (c P with c != 1 fails), complementing the curved, stationary check in check_dirac_star.py.
Run:  python3 derivations/check_dirac_selfint_3p1.py   (a few seconds)
"""
import sympy as sp

r = sp.symbols('r', positive=True)
mu, lam, kap, k, eps = sp.symbols('mu lambda kappa k epsilon', real=True)
I = sp.I
uR, uI, vR, vI = [sp.Function(n, real=True)(r) for n in ('u_R', 'u_I', 'v_R', 'v_I')]
D = lambda f: sp.diff(f, r)

S = (uR**2 + uI**2 - vR**2 - vI**2) / (2*sp.pi)
U = lam/2*S**2 + kap*S**3
P = S*(lam*S + 3*kap*S**2) - U
mu_eff = mu + lam*S + 3*kap*S**2

# flat-space evolution equations (BSSN form with alpha = chi = h = 1, beta = K = 0: L = 0, cz = cw = 1), real and imaginary
# parts of  d_t u = -i mu_eff u + i (v' + 2 v / r),   d_t v = i mu_eff v - i u'
dt_fields = {uR: mu_eff*uI - (D(vI) + 2*vI/r), uI: -mu_eff*uR + D(vR) + 2*vR/r,
             vR: -mu_eff*vI + D(uI), vI: mu_eff*vR - D(uR)}

def d_t(e):
    """time derivative of an expression in u, v and their r-derivatives, along the evolution equations"""
    shifted = e.subs({f: f + eps*df for f, df in dt_fields.items()}, simultaneous=True)
    return sp.diff(shifted, eps).subs(eps, 0).doit()

def sources(kcoef):
    Srr_kin = (vR*D(uR) + vI*D(uI) - uR*D(vR) - uI*D(vI)) / (2*sp.pi)     # Re(vbar u' - u vbar') / (2 pi)
    Stt_kin = -(uR*vR + uI*vI) / (2*sp.pi*r)                             # -Re(u vbar) / (2 pi r)
    rho = mu*S + U + P - kcoef*P + Srr_kin + 2*Stt_kin                    # rho_kin - c P  (rho_kin = mu_eff S + kinetic)
    j = (uR*D(uI) - uI*D(uR) + vR*D(vI) - vI*D(vR)) / (2*sp.pi)           # Im(ubar u' + vbar v') / (2 pi)
    return rho, j, Srr_kin + kcoef*P, Stt_kin + kcoef*P

def residuals(kcoef):
    rho, j, Srr, Stt = sources(kcoef)
    e_res = d_t(rho) + sp.diff(r**2*j, r)/r**2
    p_res = d_t(j) + sp.diff(Srr, r) + 2*(Srr - Stt)/r
    return sp.simplify(sp.expand(e_res)), sp.simplify(sp.expand(p_res))

e0, p0 = residuals(1)
assert e0 == 0, e0
assert p0 == 0, p0
print('energy and momentum conservation with rho = mu S + U + S_rr^kin + 2 S_thth^kin and S_ij += g_ij P, P = S U\' - U: OK')
ek, pk = residuals(k)
for res in (ek, pk):
    # the residual is (c - 1) times a nonzero expression
    assert res != 0 and sp.simplify(res.subs(k, 1)) == 0 and sp.simplify(sp.diff(res, k, 2)) == 0
print('conservation fixes the normalisation of the g_ab P term: OK')
