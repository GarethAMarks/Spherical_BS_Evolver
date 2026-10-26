"""
Symbolic derivation/check of the equations solved by DiracStar (Einstein-Dirac stars), in the code's metric convention
and the (-,+,+,+) signature:

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,     alpha = exp(phi),  r areal,  D = 4.

Matter: two Dirac spinors of mass mu, occupying the j = 1/2 (kappa = -1, s-wave) states with m_j = +1/2 and -1/2,
which gives a spherically symmetric stress tensor. Conventions:
    Clifford algebra {gamma^a, gamma^b} = 2 eta^{ab}, eta = diag(-1, 1, 1, 1)  (gamma^a = -i gamma_Dirac^a),
    Dirac equation   (gamma^a nabla_a + mu) psi = 0,   L = -psibar (gamma^a nabla_a + mu) psi,  psibar = psi^dag i gamma^0,
    stress tensor    T_ab = (1/4) (X_ab + X_ba),  X_ab = psibar gamma_a nabla_b psi - (nabla_a psibar) gamma_b psi
                     (the overall factor is fixed so that rho = T_00 = (omega/alpha) psi^dag psi > 0, as it must be for
                      positive-frequency stationary states),
    number density   n = sum_m psi^dag psi (orthonormal frame), Noether charge Q = int n sqrt(gamma) d^3x.
Spinor ansatz (orthonormal tetrad aligned with Cartesian axes, standard Dirac representation):
    psi_m = e^{-i omega t} alpha^{-1/2} r^{-1} ( G(r) Omega_{-1,m},  i F(r) Omega_{+1,m} ),
    Omega_{-1,m} = Y_00 chi_m,   Omega_{+1,m} = -(sigma . n) Omega_{-1,m}.

The script
  1. computes the spin connection and checks it vanishes in flat space (alpha = X = 1);
  2. inserts the ansatz into the Dirac equation and checks that it separates into the radial equations
         G' = X [  G/r + (mu + omega/alpha) F ]
         F' = X [ -F/r + (mu - omega/alpha) G ]
     (the alpha^{-1/2} and 1/r factors absorb the spin-connection terms; in flat space these are the textbook
      kappa = -1 radial equations);
  3. computes T_ab summed over both spinors, checks it is static (T_0i = 0) and spherically symmetric, and reads off
         rho = (omega/alpha) n,
         p_r = [ (omega/alpha)(F^2 + G^2) - mu (G^2 - F^2) + 2 F G / r ] / (2 pi alpha r^2),
         p_t = -F G / (2 pi alpha r^3),
     with n = (F^2 + G^2) / (2 pi alpha r^2);
  4. checks conservation  p_r' + (rho + p_r) alpha'/alpha + 2 (p_r - p_t)/r = 0  identically, for arbitrary alpha, X;
  5. prints the Einstein equations (G = 1) used by the solver:
         X'   = X (1 - X^2) / (2 r) + 4 pi r X^3 rho
         phi' = (X^2 - 1) / (2 r) + 4 pi r X^2 p_r.
Run:  python3 derivations/check_dirac_star.py   (takes about a minute)
"""
import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi', positive=True)
w, mu = sp.symbols('omega mu', positive=True)
al = sp.Function('alpha')(r)
X = sp.Function('X')(r)
F = sp.Function('F')(r)
G = sp.Function('G')(r)
coords = [t, r, th, ph]
I = sp.I

# ---------------- metric, tetrad, spin connection ----------------
g = sp.diag(-al**2, X**2, r**2, r**2 * sp.sin(th)**2)
ginv = g.inv()
eta = sp.diag(-1, 1, 1, 1)

def christoffel(g, ginv):
    Gam = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                Gam[a][b][c] = (sum(ginv[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b])
                                                           - sp.diff(g[b, c], coords[d])) for d in range(4)) / 2)
    return Gam

Gam = christoffel(g, ginv)

nvec = [sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)]
thvec = [sp.cos(th)*sp.cos(ph), sp.cos(th)*sp.sin(ph), -sp.sin(th)]
phvec = [-sp.sin(ph), sp.cos(ph), 0]
# e_a^mu (rows a = 0..3, columns mu): Cartesian-aligned orthonormal tetrad
E = sp.zeros(4, 4)
E[0, 0] = 1/al
for i in range(3):
    E[i+1, 1] = nvec[i] / X
    E[i+1, 2] = thvec[i] / r
    E[i+1, 3] = phvec[i] / (r * sp.sin(th))
# co-tetrad e^a_mu = eta^{ab} g_{mu nu} e_b^nu
Ecov = sp.zeros(4, 4)  # e^a_mu
for a in range(4):
    for m in range(4):
        Ecov[a, m] = sum(eta[a, b] * g[m, n] * E[b, n] for b in range(4) for n in range(4))

# omega_mu^a_b = e^a_nu (d_mu e_b^nu + Gamma^nu_{mu lam} e_b^lam);  omega_{mu ab} = eta_ac omega_mu^c_b
omega = [[[0]*4 for _ in range(4)] for _ in range(4)]
for m in range(4):
    for a in range(4):
        for b in range(4):
            omega[m][a][b] = (sum(Ecov[a, n] * (sp.diff(E[b, n], coords[m])
                                         + sum(Gam[n][m][l] * E[b, l] for l in range(4))) for n in range(4)))
omega_low = [[[sum(eta[a, c] * omega[m][c][b] for c in range(4)) for b in range(4)] for a in range(4)] for m in range(4)]

# After differentiating, expressions are evaluated at generic points with rational trigonometric values and with the radial
# functions and their derivatives replaced by symbols, so that every check reduces to rational-function arithmetic.
a_, x_ = sp.symbols('a x', positive=True)
da_, dx_, F_, dF_, G_, dG_ = sp.symbols('da dx F_ dF_ G_ dG_', real=True)
def at_point(e, pt):
    e = sp.sympify(e)
    e = e.subs({sp.Derivative(al, r): da_, sp.Derivative(X, r): dx_, sp.Derivative(F, r): dF_, sp.Derivative(G, r): dG_})
    e = e.subs({al: a_, X: x_, F: F_, G: G_})
    e = e.subs({sp.sin(th): pt[0], sp.cos(th): pt[1], sp.sin(ph): pt[2], sp.cos(ph): pt[3]})
    return e
PTS = [(sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(5, 13), sp.Rational(12, 13)),
       (sp.Rational(8, 17), sp.Rational(-15, 17), sp.Rational(-7, 25), sp.Rational(24, 25))]
def is_zero(e):
    return sp.cancel(sp.together(sp.expand(e))) == 0

flat_ok = all(is_zero(at_point(omega_low[m][a][b], pt).subs({a_: 1, da_: 0, x_: 1, dx_: 0}))
              for pt in PTS for m in range(4) for a in range(4) for b in range(4))
assert flat_ok, 'spin connection should vanish in flat space for a Cartesian tetrad'
print('1. spin connection vanishes in flat space: OK')

# ---------------- gamma matrices ----------------
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]])
sig = [s1, s2, s3]
Z2, I2 = sp.zeros(2), sp.eye(2)
def block(a, b, c, d):
    M = sp.zeros(4)
    M[:2, :2] = a; M[:2, 2:] = b; M[2:, :2] = c; M[2:, 2:] = d
    return M
gD = [block(I2, Z2, Z2, -I2)] + [block(Z2, s, -s, Z2) for s in sig]   # Dirac rep, eta_D = diag(1,-1,-1,-1)
gam = [-I * gm for gm in gD]                                           # mostly-plus Clifford algebra
for a in range(4):
    for b in range(4):
        assert sp.expand(gam[a]*gam[b] + gam[b]*gam[a] - 2*eta[a, b]*sp.eye(4)) == sp.zeros(4)
gam_low = [sum((eta[a, b] * gam[b] for b in range(4)), sp.zeros(4)) for a in range(4)]
beta_mat = I * gam[0]   # psibar = psi^dag beta

# Gamma_mu = (1/8) omega_{mu ab} [gamma^a, gamma^b]
Gspin = [sum(((omega_low[m][a][b] / 8) * (gam[a]*gam[b] - gam[b]*gam[a]) for a in range(4) for b in range(4)), sp.zeros(4))
         for m in range(4)]

# ---------------- spinor ansatz ----------------
Y00 = 1 / sp.sqrt(4*sp.pi)
sig_n = sum((nvec[i] * sig[i] for i in range(3)), sp.zeros(2))
chis = [sp.Matrix([1, 0]), sp.Matrix([0, 1])]

def psi_m(chi):
    up = G * Y00 * chi
    lo = I * F * (-(sig_n * Y00 * chi))
    v = sp.Matrix([up[0], up[1], lo[0], lo[1]])
    return sp.exp(-I*w*t) * v / (sp.sqrt(al) * r)

def nabla_coord(psi):
    return [sp.diff(psi, coords[m]) + Gspin[m] * psi for m in range(4)]

def nabla_frame(psi):
    nc = nabla_coord(psi)
    return [sum((E[a, m] * nc[m] for m in range(4)), sp.zeros(4, 1)) for a in range(4)]

# radial equations to be verified
dG = X * (G/r + (mu + w/al) * F)
dF = X * (-F/r + (mu - w/al) * G)
on_shell = {sp.Derivative(G, r): dG, sp.Derivative(F, r): dF}

def simp(e):
    return sp.cancel(sp.together(sp.expand(e)))
dG_val, dF_val = at_point(dG, PTS[0]), at_point(dF, PTS[0])   # no angular dependence

dirac_ok = True
nab = []
for chi in chis:
    psi = psi_m(chi)
    nf = nabla_frame(psi)
    nab.append((psi, nf))
    res = sum((gam[a] * nf[a] for a in range(4)), sp.zeros(4, 1)) + mu * psi
    for pt in PTS:
        for comp in res:
            c = at_point(comp * sp.exp(I*w*t), pt).subs({dG_: dG_val, dF_: dF_val})
            if not is_zero(c):
                dirac_ok = False
print('2. Dirac equation separates into  G\' = X[G/r + (mu + w/alpha) F],  F\' = X[-F/r + (mu - w/alpha) G]:', 'OK' if dirac_ok else 'FAILED')
assert dirac_ok

# ---------------- stress tensor ----------------
def dagger(v):
    return v.conjugate().T

T = sp.zeros(4)
n_dens = 0
for psi, nf in nab:
    psibar = dagger(psi) * beta_mat
    nfbar = [dagger(v) * beta_mat for v in nf]
    Xab = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            Xab[a, b] = (psibar * gam_low[a] * nf[b])[0] - (nfbar[a] * gam_low[b] * psi)[0]
    T += (Xab + Xab.T) / 4
    n_dens += (dagger(psi) * psi)[0]

def ev(e, pt):
    return simp(at_point(e, pt).subs({dG_: dG_val, dF_: dF_val}))
pt0 = PTS[0]
vec_at = lambda v: [at_point(c, pt0) for c in v]
nn, tt, pp = vec_at(nvec), vec_at(thvec), vec_at(phvec)
Tn = sp.Matrix(4, 4, lambda a, b: ev(T[a, b], pt0))
rho = simp(Tn[0, 0])
p_r = simp(sum(nn[i]*nn[j]*Tn[i+1, j+1] for i in range(3) for j in range(3)))
p_t = simp(sum(tt[i]*tt[j]*Tn[i+1, j+1] for i in range(3) for j in range(3)))
p_f = simp(sum(pp[i]*pp[j]*Tn[i+1, j+1] for i in range(3) for j in range(3)))
off = [simp(Tn[0, i]) for i in range(1, 4)] + [simp(sum(nn[i]*tt[j]*Tn[i+1, j+1] for i in range(3) for j in range(3))),
                                               simp(sum(nn[i]*pp[j]*Tn[i+1, j+1] for i in range(3) for j in range(3))),
                                               simp(sum(tt[i]*pp[j]*Tn[i+1, j+1] for i in range(3) for j in range(3)))]
n_dens = ev(n_dens, pt0)

assert all(o == 0 for o in off), off
assert simp(p_t - p_f) == 0
# angular independence: same values at the second point
assert simp(ev(T[0, 0], PTS[1]) - rho) == 0
assert simp(ev(sum(at_point(nvec[i], PTS[1]) * at_point(nvec[j], PTS[1]) * T[i+1, j+1] for i in range(3) for j in range(3)), PTS[1]) - p_r) == 0
print('3. T_ab static and spherically symmetric (T_0i = 0, off-diagonal spatial parts vanish, p_theta = p_phi): OK')

n_expected = (F_**2 + G_**2) / (2*sp.pi*a_*r**2)
assert simp(n_dens - n_expected) == 0
assert simp(rho - (w/a_) * n_expected) == 0
pr_expected = ((w/a_) * (F_**2 + G_**2) - mu * (G_**2 - F_**2) + 2*F_*G_/r) / (2*sp.pi*a_*r**2)
pt_expected = -(F_*G_ / r) / (2*sp.pi*a_*r**2)
assert simp(p_r - pr_expected) == 0, sp.factor(p_r)
assert simp(p_t - pt_expected) == 0, ('p_t', sp.factor(p_t))
print('   n   = (F^2 + G^2) / (2 pi alpha r^2)')
print('   rho = (omega/alpha) n')
print('   p_r = [ (omega/alpha)(F^2 + G^2) - mu (G^2 - F^2) + 2 F G / r ] / (2 pi alpha r^2)')
print('   p_t = -F G / (2 pi alpha r^3)')

# ---------------- conservation ----------------
# d/dr along a solution: chain rule through a, F_, G_ (alpha' = da arbitrary, X arbitrary)
dpr = sp.diff(pr_expected, r) + sp.diff(pr_expected, a_) * da_ + sp.diff(pr_expected, F_) * dF_val + sp.diff(pr_expected, G_) * dG_val
cons = simp(dpr + (rho + pr_expected) * da_ / a_ + 2 * (pr_expected - pt_expected) / r)
assert cons == 0, cons
print('4. p_r\' + (rho + p_r) alpha\'/alpha + 2 (p_r - p_t)/r = 0 identically (any alpha, X): OK')

print('5. Einstein equations (G = 1):  X\' = X (1 - X^2)/(2r) + 4 pi r X^3 rho,   phi\' = (X^2 - 1)/(2r) + 4 pi r X^2 p_r')
print('   Noether charge Q = int n X r^2 dOmega dr = int 2 (F^2 + G^2) X / alpha dr')
