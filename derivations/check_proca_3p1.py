"""
Symbolic check that the 3+1 Proca evolution equations used in ComplexVectorField follow from the covariant
field equations  nabla_a F^{ab} = mu^2 A^b  (with the Lorenz condition nabla_a A^a = 0 they imply), for a general
spherically symmetric D-dimensional spacetime with arbitrary lapse, radial shift and time-dependent metric:

    ds^2 = -alpha^2 dt^2 + g_rr (dr + beta dt)^2 + g_T r^2 dOmega_{n}^2 ,   n = D - 2 (symbolic).

Variables (one real Proca field; the complex field is two copies):
    X_r = gamma_r^a A_a,   X_phi = -n^a A_a,   E^r = gamma^{rr} F_{r a} n^a.
Equations checked (Z = 0):
    d_t X_r   = beta d_r X_r + X_r d_r beta - alpha (gamma_rr E^r + d_r X_phi) - X_phi d_r alpha
    d_t X_phi = beta d_r X_phi + alpha K X_phi - alpha D_i X^i - X^i d_i alpha
    d_t E^r   = beta d_r E^r - E^r d_r beta + alpha K E^r + alpha mu^2 X^r
where the Gauss constraint is C = D_i E^i + mu^2 X_phi = 0. The E equation equals the exact field equation plus
beta * C (exact - ours = -beta C), i.e. they agree on the constraint surface; this is the standard choice (as in
Zilhao, Witek & Cardoso 2015) and the script verifies the constraint term explicitly.
Run:  python3 derivations/check_proca_3p1.py
"""
import sympy as sp

t, r = sp.symbols('t r', positive=True)
n, mu = sp.symbols('n mu', positive=True)
F = lambda name: sp.Function(name)(t, r)
alpha, beta, grr, gT, At, Ar = F('alpha'), F('beta'), F('g_rr'), F('g_T'), F('A_t'), F('A_r')

# (t, r) block of the inverse 4-metric
gtt_up, gtr_up, grr_up = -1/alpha**2, beta/alpha**2, 1/grr - beta**2/alpha**2
sqrt_gamma = sp.sqrt(grr) * (gT * r**2)**(n/2)
sqrt_g = alpha * sqrt_gamma

F_tr = sp.diff(Ar, t) - sp.diff(At, r)
F_up_tr = (gtt_up*grr_up - gtr_up**2) * F_tr          # F^{tr}
A_up_t = gtt_up*At + gtr_up*Ar
A_up_r = gtr_up*At + grr_up*Ar

lorenz = (sp.diff(sqrt_g*A_up_t, t) + sp.diff(sqrt_g*A_up_r, r)) / sqrt_g
eq_r = sp.diff(sqrt_g*F_up_tr, t)/sqrt_g - mu**2*A_up_r     # b = r component of nabla_a F^{ab} - mu^2 A^b
eq_t = sp.diff(sqrt_g*(-F_up_tr), r)/sqrt_g - mu**2*A_up_t  # b = t component

# 3+1 variables
X = Ar
P = -(At - beta*Ar)/alpha                     # X_phi = -n^a A_a, n_a = (-alpha, 0)
E_dn = -F_tr/alpha                            # E_r = F_{r a} n^a = F_{rt}/alpha
E = E_dn/grr                                  # E^r
X_up = X/grr

# extrinsic curvature trace from the actual metric time derivative: K_ij = -(d_t gamma_ij - Lie_beta gamma_ij)/(2 alpha)
K_rr = -(sp.diff(grr, t) - (beta*sp.diff(grr, r) + 2*grr*sp.diff(beta, r)))/(2*alpha)
gTT = gT*r**2
K_TT = -(sp.diff(gTT, t) - beta*sp.diff(gTT, r))/(2*alpha)
K = K_rr/grr + n*K_TT/gTT

div = lambda Vr: sp.diff(sqrt_gamma*Vr, r)/sqrt_gamma
gauss = div(E) + mu**2*P

# eliminate second time derivatives / d_t A_t using the field equations
dtAt = sp.solve(lorenz, sp.diff(At, t))[0]
dt2Ar = sp.solve(eq_r, sp.diff(Ar, t, 2))[0]
def on_shell(expr):
    expr = sp.sympify(expr)
    expr = expr.subs(sp.diff(Ar, t, 2), dt2Ar)
    expr = expr.subs(sp.diff(At, t, r), sp.diff(dtAt, r))
    return expr.subs(sp.diff(At, t), dtAt)

checks = {
    'd_t X_r':   (sp.diff(X, t), beta*sp.diff(X, r) + X*sp.diff(beta, r) - alpha*(grr*E + sp.diff(P, r)) - P*sp.diff(alpha, r), 0),
    'd_t X_phi': (sp.diff(P, t), beta*sp.diff(P, r) + alpha*K*P - alpha*div(X_up) - X_up*sp.diff(alpha, r), 0),
    'd_t E^r':   (sp.diff(E, t), beta*sp.diff(E, r) - E*sp.diff(beta, r) + alpha*K*E + alpha*mu**2*X_up, -beta*gauss),
}
ok = True
for name, (exact, ours, allowed) in checks.items():
    diff = sp.simplify(on_shell(exact - ours) - on_shell(allowed))
    ok &= (diff == 0)
    print("%-10s  exact - ours - (%s) = %s" % (name, '-beta * C' if allowed != 0 else '0', diff))

# Gauss constraint is the t component of the field equation:  eq_t = -C / alpha
g_diff = sp.simplify(eq_t + gauss/alpha)
ok &= (g_diff == 0)
print("Gauss law: eq_t + C/alpha = %s" % g_diff)
print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
