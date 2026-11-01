"""
Linearised NON-RADIAL (polar, axisymmetric M = 0, total angular momentum L) perturbation equations of free spherical Dirac
stars (two spinors, m_j = +-1/2, no self-interaction), signature (-,+,+,+), G = 1. Conventions as in check_dirac_star.py
and dirac_pulsation.py:

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,  alpha = exp(phi0),
    psi_m = e^{-i w t} ( g(r) Y00 chi_m, -i f(r) (sigma.n) Y00 chi_m )      (f, g the physical amplitudes of dirac_relax.py),
    Cartesian-aligned orthonormal tetrad, gamma^a = -i gamma_Dirac^a, psibar = psi^dag i gamma^0.

Symmetry. The background is invariant under a rotation combined with the U(2) flavour rotation of the doublet (psi_+, psi_-),
so perturbations carry an integer total angular momentum L. For M = 0 and polar parity (-1)^L the + spinor perturbation is
    dpsi_+ = e^{-i w t} ( a1 Y chi_up + a2 e^{i phi} Y' chi_dn,  i (sigma.n)(b1 Y chi_up + b2 e^{i phi} Y' chi_dn) ),
    Y = P_L(cos theta), Y' = dY/dtheta,
i.e. the upper component lies in the (l = L, m_j = 1/2) space spanned by the kappa = -(L+1) and kappa = L channels, and the
lower one in its image under sigma.n. a_i, b_i are complex: a_i = aR_i + i aI_i etc., real functions of (t, r). The -
spinor is fixed by the symmetry R_x(pi) x (flavour swap), under which the polar metric perturbation is even/odd as (-1)^L:
    dpsi_- = e^{-i w t} ( a1 Y chi_dn + s2 a2 e^{-i phi} Y' chi_up,  i (sigma.n)(b1 Y chi_dn + s2 b2 e^{-i phi} Y' chi_up) ),
with s2 = -1 (checked below: s2 = +1 leaves axial parts in T_ab, i.e. mixes in the L +- 1 axial sectors).

Metric perturbation (general polar, used for the gauge check; Regge-Wheeler gauge h0 = h1 = Gs = 0 for the mode solver):
    h_tt = alpha^2 H0 Y,  h_tr = H1 Y,  h_rr = X^2 H2 Y,  h_tA = h0 D_A Y,  h_rA = h1 D_A Y,
    h_AB = r^2 (K Y Omega_AB + Gs D_A D_B Y),
with a symmetric tetrad perturbation  de_a^mu = -(1/2) h^mu_nu e_a^nu.

Derivation: everything is linearised in epsilon; derivatives in theta are taken symbolically with Y'' = -cot Y' - Lam Y
(Lam = L(L+1)), and the results evaluated at generic rational points on the sphere, with Y(theta0) -> y0, Y'(theta0) -> y1
treated as independent (valid for general L). The projections are
    Dirac (+ spinor):   upper = A1 Y chi_up + A2 e^{i phi} Y' chi_dn,   (sigma.n) lower = C1 Y chi_up + C2 e^{i phi} Y' chi_dn,
    Einstein E_mn = dG_mn - 8 pi dT_mn:  E_tt = Ett Y, E_tr = Etr Y, E_rr = Err Y, E_tth = Etth Y', E_rth = Erth Y',
    E_thth = Ea Y + Eb Y'',  E_phph / sin^2 = Ea Y + Eb cot Y'.
Checks:
  1. background: zeroth-order Dirac and Einstein equations vanish with the dirac_relax.py background equations;
  2. closure: the projections hold (no leftover angular structure) at two different points, the axial components
     E_tph, E_rph, E_thph vanish, and the - spinor's Dirac residual is the symmetry image of the + one;
  3. gauge: a pure-gauge perturbation h = L_zeta g0, dpsi = Kosmann Lie derivative of psi0 along a polar zeta,
     zeta = (T Y, R Y, S Y' / r^2, 0) with arbitrary T(t,r), R(t,r), S(t,r), satisfies all projected equations identically.
Output: relaxation/dirac_nonradial_coeffs.py with, for each projected equation in the RW gauge, the coefficients of the
fields y = (aR1, aI1, aR2, aI2, bR1, bI1, bR2, bI2, H0, H1, H2, K) and their derivatives d_t, d_r, d_tt, d_tr, d_rr,
as functions of (r, w, mu, Lam, g, f, X, alpha).
Run: python3 derivations/dirac_nonradial.py   (several minutes)
"""
import os, sys, time
import sympy as sp

T0 = time.time()
def log(*a):
    print('[%5.0fs]' % (time.time() - T0), *a, flush=True)

t, r, th, ph = sp.symbols('t r theta phi_', real=True)
coords = [t, r, th, ph]
w, mu, Lam = sp.symbols('omega mu Lambda', real=True)
I = sp.I
phi0 = sp.Function('phi0', real=True)(r)
al = sp.exp(phi0)
X = sp.Function('X', real=True)(r)
gb = sp.Function('g', real=True)(r)
fb = sp.Function('f', real=True)(r)
Y = sp.Function('Y', real=True)(th)
Yp = sp.diff(Y, th)
s_, c_ = sp.sin(th), sp.cos(th)

def Fn(name):
    return sp.Function(name, real=True)(t, r)

MAT = ['aR1', 'aI1', 'aR2', 'aI2', 'bR1', 'bI1', 'bR2', 'bI2']
MET = ['H0', 'H1', 'H2', 'K']
GAU = ['h0', 'h1', 'Gs']
fld = {n: Fn(n) for n in MAT + MET + GAU}
a1 = fld['aR1'] + I * fld['aI1']; a2 = fld['aR2'] + I * fld['aI2']
b1 = fld['bR1'] + I * fld['bI1']; b2 = fld['bR2'] + I * fld['bI2']

# ---------------- background ----------------
g0 = sp.diag(-al**2, X**2, r**2, r**2 * s_**2)
g0i = sp.diag(-1 / al**2, 1 / X**2, 1 / r**2, 1 / (r**2 * s_**2))
eta = sp.diag(-1, 1, 1, 1)

def christoffel0():
    G = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                G[a][b][c] = sp.expand(sum(g0i[a, d] * (sp.diff(g0[d, b], coords[c]) + sp.diff(g0[d, c], coords[b])
                                                     - sp.diff(g0[b, c], coords[d])) for d in range(4)) / 2)
    return G
Gam0 = christoffel0()

nvec = [s_ * sp.cos(ph), s_ * sp.sin(ph), c_]
thvec = [c_ * sp.cos(ph), c_ * sp.sin(ph), -s_]
phvec = [-sp.sin(ph), sp.cos(ph), 0]
E0 = sp.zeros(4, 4)
E0[0, 0] = 1 / al
for i in range(3):
    E0[i + 1, 1] = nvec[i] / X
    E0[i + 1, 2] = thvec[i] / r
    E0[i + 1, 3] = phvec[i] / (r * s_)

s1 = sp.Matrix([[0, 1], [1, 0]]); s2m = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]])
sig = [s1, s2m, s3]
Z2, I2 = sp.zeros(2), sp.eye(2)
def block(a, b, c, d):
    M = sp.zeros(4)
    M[:2, :2] = a; M[:2, 2:] = b; M[2:, :2] = c; M[2:, 2:] = d
    return M
gam = [-I * gm for gm in [block(I2, Z2, Z2, -I2)] + [block(Z2, s, -s, Z2) for s in sig]]
gam_low = [sum((eta[a, b] * gam[b] for b in range(4)), sp.zeros(4)) for a in range(4)]
beta = I * gam[0]
sig_n = sum((nvec[i] * sig[i] for i in range(3)), sp.zeros(2))
up, dn = sp.Matrix([1, 0]), sp.Matrix([0, 1])
Y00 = 1 / sp.sqrt(4 * sp.pi)

def spinor(upper, lower):
    return sp.Matrix([upper[0], upper[1], lower[0], lower[1]])

ph_t = sp.exp(-I * w * t)
psi0 = {+1: ph_t * spinor(gb * Y00 * up, -I * fb * sig_n * Y00 * up),
        -1: ph_t * spinor(gb * Y00 * dn, -I * fb * sig_n * Y00 * dn)}

def dpsi(s2):
    Bp = (a1 * Y * up + a2 * sp.exp(I * ph) * Yp * dn, b1 * Y * up + b2 * sp.exp(I * ph) * Yp * dn)
    Bm = (a1 * Y * dn + s2 * a2 * sp.exp(-I * ph) * Yp * up, b1 * Y * dn + s2 * b2 * sp.exp(-I * ph) * Yp * up)
    return {+1: ph_t * spinor(Bp[0], I * sig_n * Bp[1]), -1: ph_t * spinor(Bm[0], I * sig_n * Bm[1])}

# ---------------- perturbed metric, tetrad, connection (first order) ----------------
def metric_pert(H0, H1, H2, K, h0, h1, Gs):
    h = sp.zeros(4)
    h[0, 0] = al**2 * H0 * Y
    h[0, 1] = h[1, 0] = H1 * Y
    h[1, 1] = X**2 * H2 * Y
    h[0, 2] = h[2, 0] = h0 * Yp
    h[1, 2] = h[2, 1] = h1 * Yp
    h[2, 2] = r**2 * (K * Y + Gs * sp.diff(Y, th, 2))
    h[3, 3] = r**2 * (K * Y * s_**2 + Gs * s_ * c_ * Yp)
    return h

def first_order(h):
    hup = g0i * h                     # h^mu_nu
    Gam1 = [[[sp.expand(-sum(hup[a, e] * Gam0[e][b][c] for e in range(4))
                        + sum(g0i[a, d] * (sp.diff(h[d, c], coords[b]) + sp.diff(h[d, b], coords[c]) - sp.diff(h[b, c], coords[d]))
                              for d in range(4)) / 2) for c in range(4)] for b in range(4)] for a in range(4)]
    E1 = sp.zeros(4, 4)
    for a in range(4):
        for m in range(4):
            E1[a, m] = -sum(hup[m, n] * E0[a, n] for n in range(4)) / 2
    return hup, Gam1, E1

Ecov0 = sp.Matrix(4, 4, lambda a, m: sum(eta[a, b] * g0[m, n] * E0[b, n] for b in range(4) for n in range(4)))

def spin_connection(E, Ecov, Gam):
    """omega_{mu ab} lowered; returns Gamma_mu = (1/8) omega_{mu ab} [gamma^a, gamma^b]"""
    out = []
    for m in range(4):
        M = sp.zeros(4)
        for a in range(4):
            for b in range(a + 1, 4):
                w_ab = sum(eta[a, a] * Ecov[a, n] * (sp.diff(E[b, n], coords[m]) + sum(Gam[n][m][l] * E[b, l] for l in range(4)))
                           for n in range(4))
                if w_ab != 0:
                    M += (w_ab / 4) * (gam[a] * gam[b] - gam[b] * gam[a])   # (1/8)(w_ab [a,b] + w_ba [b,a]) with w_ba = -w_ab
        out.append(M)
    return out

log('background spin connection')
Gs0 = spin_connection(E0, Ecov0, Gam0)

def linear_connection(h, hup, Gam1, E1):
    """first-order spin connection: d/deps of Gamma_mu(E0 + eps E1, Gam0 + eps Gam1)"""
    Ecov1 = sp.Matrix(4, 4, lambda a, m: sum(eta[a, b] * (g0[m, n] * E1[b, n] + h[m, n] * E0[b, n]) for b in range(4) for n in range(4)))
    out = []
    for m in range(4):
        M = sp.zeros(4)
        for a in range(4):
            for b in range(a + 1, 4):
                w1 = 0
                for n in range(4):
                    d0 = sp.diff(E0[b, n], coords[m]) + sum(Gam0[n][m][l] * E0[b, l] for l in range(4))
                    d1 = sp.diff(E1[b, n], coords[m]) + sum(Gam0[n][m][l] * E1[b, l] + Gam1[n][m][l] * E0[b, l] for l in range(4))
                    w1 += eta[a, a] * (Ecov1[a, n] * d0 + Ecov0[a, n] * d1)
                if w1 != 0:
                    M += (w1 / 4) * (gam[a] * gam[b] - gam[b] * gam[a])
        out.append(M)
    return Ecov1, out

# ---------------- evaluation at points ----------------
PTS = [(sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(5, 13), sp.Rational(12, 13)),
       (sp.Rational(8, 17), sp.Rational(-15, 17), sp.Rational(-7, 25), sp.Rational(24, 25))]
y0, y1 = sp.symbols('y0 y1', real=True)
# background derivative rules (free field, dirac_relax.py)
gS, fS, XS, aS = sp.symbols('g_ f_ X_ al_', real=True)
def bg_rules():
    n = (gb**2 + fb**2) / (2 * sp.pi); S = (gb**2 - fb**2) / (2 * sp.pi)
    rho = (w / al) * n
    pr = (w / al) * n - mu * S + fb * gb / (sp.pi * r)
    dphi = (X**2 - 1) / (2 * r) + 4 * sp.pi * r * X**2 * pr
    dX = X * (1 - X**2) / (2 * r) + 4 * sp.pi * r * X**3 * rho
    dg = (X - 1) * gb / r - dphi * gb / 2 + X * (mu + w / al) * fb
    df = -(1 + X) * fb / r - dphi * fb / 2 + X * (mu - w / al) * gb
    first = {sp.Derivative(phi0, r): dphi, sp.Derivative(X, r): dX, sp.Derivative(gb, r): dg, sp.Derivative(fb, r): df}
    rules = [dict(first)]
    cur = {F: e for F, e in [(phi0, dphi), (X, dX), (gb, dg), (fb, df)]}
    for n in range(2, 5):
        nxt = {F: sp.diff(e, r).subs(first) for F, e in cur.items()}
        rules.append({sp.Derivative(F, (r, n)): e for F, e in nxt.items()})
        cur = nxt
    return rules[::-1]          # highest order first
BG_RULES = bg_rules()
TO_SYM = {gb: gS, fb: fS, X: XS}

FS = {}
def field_symbols(names):
    rep = {}
    for n in names:
        F = fld[n]
        for k, dv in [('_tt', (t, 2)), ('_tr', (t, r)), ('_rr', (r, 2)), ('_t', (t,)), ('_r', (r,))]:
            sym = sp.Symbol(n + k, real=True)
            FS[n + k] = sym
            rep[sp.Derivative(F, *dv) if len(dv) > 1 and dv[0] != dv[-1] or len(dv) == 1 else sp.Derivative(F, dv)] = sym
        FS[n] = sp.Symbol(n, real=True)
    return rep

def reduce_Y(e):
    Y2 = -sp.cot(th) * Yp - Lam * Y
    Y3 = sp.diff(Y2, th).subs(sp.Derivative(Y, (th, 2)), Y2)
    Y4 = sp.diff(Y3, th).subs(sp.Derivative(Y, (th, 2)), Y2)
    e = e.subs(sp.Derivative(Y, (th, 4)), Y4).subs(sp.Derivative(Y, (th, 3)), Y3).subs(sp.Derivative(Y, (th, 2)), Y2)
    return e

def at_point(e, pt, extra_fields=()):
    e = sp.sympify(e)
    e = reduce_Y(e)
    for rule in BG_RULES:
        e = e.subs(rule)
    e = e.replace(sp.tan, lambda x: sp.sin(x) / sp.cos(x)).replace(sp.cot, lambda x: sp.cos(x) / sp.sin(x))
    def _ex(x):
        a = sp.expand(x.args[0]); k = sp.simplify(a.coeff(ph) / I)
        return sp.exp(a - I * k * ph) * (pt[3] + I * pt[2])**k
    e = e.replace(lambda x: isinstance(x, sp.exp) and x.has(ph), _ex)
    e = e.subs({sp.sin(th): pt[0], sp.cos(th): pt[1], sp.sin(ph): pt[2], sp.cos(ph): pt[3]})
    e = e.subs({sp.Derivative(Y, th): y1}).subs(Y, y0)
    assert not e.has(th) and not e.has(ph), e.atoms(sp.Function)
    # field derivatives -> symbols
    rep = {}
    for n in list(fld):
        F = fld[n]
        rep[sp.Derivative(F, (t, 2))] = sp.Symbol(n + '_tt', real=True)
        rep[sp.Derivative(F, t, r)] = sp.Symbol(n + '_tr', real=True)
        rep[sp.Derivative(F, (r, 2))] = sp.Symbol(n + '_rr', real=True)
        rep[sp.Derivative(F, (r, 3))] = sp.Symbol(n + '_rrr', real=True)
        rep[sp.Derivative(F, (t, 2), r)] = sp.Symbol(n + '_ttr', real=True)
        rep[sp.Derivative(F, t, (r, 2))] = sp.Symbol(n + '_trr', real=True)
        rep[sp.Derivative(F, (t, 3))] = sp.Symbol(n + '_ttt', real=True)
    e = e.subs(rep)
    rep = {}
    for n in list(fld):
        F = fld[n]
        rep[sp.Derivative(F, t)] = sp.Symbol(n + '_t', real=True)
        rep[sp.Derivative(F, r)] = sp.Symbol(n + '_r', real=True)
    e = e.subs(rep)
    e = e.subs({F: sp.Symbol(n, real=True) for n, F in fld.items()})
    e = e.subs(sp.exp(phi0), aS).subs(sp.exp(-phi0), 1 / aS).subs(phi0, sp.log(aS))
    e = e.subs(TO_SYM)
    assert not e.has(sp.Derivative), [d for d in e.atoms(sp.Derivative)][:3]
    return e

def simp(e):
    return sp.cancel(sp.together(sp.expand(e)))

# ---------------- linearised equations for a given perturbation ----------------
def linearised(h, dps, label):
    """returns at each point: Dirac residual (4-vectors) for both spinors and the Einstein components E_mn"""
    log(label, ': first-order geometry')
    hup, Gam1, E1 = first_order(h)
    Ecov1, Gs1 = linear_connection(h, hup, Gam1, E1)
    log(label, ': spinor derivatives')
    res = {}
    # coordinate derivatives of spinors (symbolic in angles), then frame derivatives
    dcache = {}
    for m in (+1, -1):
        dcache[m] = ([sp.diff(psi0[m], x) for x in coords], [sp.diff(dps[m], x) for x in coords])
    # Einstein tensor perturbation (symbolic), needs derivatives of Gam1
    log(label, ': linearised Ricci')
    dGam1 = [[[[sp.diff(Gam1[a][b][c], coords[e]) for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
    out = []
    for ip, pt in enumerate(PTS):
        P = lambda e: at_point(e, pt)
        log(label, ': point', ip)
        E0p, E1p = E0.applyfunc(P), E1.applyfunc(P)
        Ec0p, Ec1p = Ecov0.applyfunc(P), Ecov1.applyfunc(P)
        G0p = [M.applyfunc(P) for M in Gs0]
        G1p = [M.applyfunc(P) for M in Gs1]
        dirac, nf0, nf1 = {}, {}, {}
        for m in (+1, -1):
            p0 = psi0[m].applyfunc(P); dp = dps[m].applyfunc(P)
            d0 = [v.applyfunc(P) for v in dcache[m][0]]; d1 = [v.applyfunc(P) for v in dcache[m][1]]
            cov0 = [d0[k] + G0p[k] * p0 for k in range(4)]
            cov1 = [G1p[k] * p0 + d1[k] + G0p[k] * dp for k in range(4)]
            nf0[m] = [sum((E0p[a, k] * cov0[k] for k in range(4)), sp.zeros(4, 1)) for a in range(4)]
            nf1[m] = [sum((E1p[a, k] * cov0[k] + E0p[a, k] * cov1[k] for k in range(4)), sp.zeros(4, 1)) for a in range(4)]
            D0 = sum((gam[a] * nf0[m][a] for a in range(4)), sp.zeros(4, 1)) + mu * p0
            D1 = sum((gam[a] * nf1[m][a] for a in range(4)), sp.zeros(4, 1)) + mu * dp
            ephase = sp.exp(I * w * t)
            dirac[m] = ([simp(x * ephase) for x in D0], [simp(x * ephase) for x in D1], (p0, dp))
        # stress tensor (frame), first order
        T0f, T1f = sp.zeros(4), sp.zeros(4)
        for m in (+1, -1):
            p0, dp = dirac[m][2]
            pb0 = p0.conjugate().T * beta; pb1 = dp.conjugate().T * beta
            for a in range(4):
                for b in range(4):
                    nb0 = nf0[m][b]; nb1 = nf1[m][b]
                    na0b = nf0[m][a].conjugate().T * beta; na1b = nf1[m][a].conjugate().T * beta
                    x0 = (pb0 * gam_low[a] * nb0)[0] - (na0b * gam_low[b] * p0)[0]
                    x1 = ((pb1 * gam_low[a] * nb0)[0] + (pb0 * gam_low[a] * nb1)[0]
                          - (na1b * gam_low[b] * p0)[0] - (na0b * gam_low[b] * dp)[0])
                    T0f[a, b] += x0 / 4; T0f[b, a] += x0 / 4
                    T1f[a, b] += x1 / 4; T1f[b, a] += x1 / 4
        T1c = sp.zeros(4)
        for mm in range(4):
            for nn in range(4):
                T1c[mm, nn] = simp(sum(Ec1p[a, mm] * Ec0p[b, nn] * T0f[a, b] + Ec0p[a, mm] * Ec1p[b, nn] * T0f[a, b]
                                       + Ec0p[a, mm] * Ec0p[b, nn] * T1f[a, b] for a in range(4) for b in range(4)))
        T0c = sp.Matrix(4, 4, lambda mm, nn: simp(sum(Ec0p[a, mm] * Ec0p[b, nn] * T0f[a, b] for a in range(4) for b in range(4))))
        # Ricci: background and first order
        Gm0 = [[[P(Gam0[a][b][c]) for c in range(4)] for b in range(4)] for a in range(4)]
        Gm1 = [[[P(Gam1[a][b][c]) for c in range(4)] for b in range(4)] for a in range(4)]
        dG0 = [[[[P(sp.diff(Gam0[a][b][c], coords[e])) for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
        dG1 = [[[[P(dGam1[a][b][c][e]) for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
        def ricci(dG, GA, GB, GC=None):
            pass
        R0 = sp.zeros(4); R1 = sp.zeros(4)
        for b in range(4):
            for c in range(4):
                R0[b, c] = sum(dG0[a][b][c][a] - dG0[a][b][a][c] for a in range(4)) + sum(
                    Gm0[a][a][d] * Gm0[d][b][c] - Gm0[a][c][d] * Gm0[d][b][a] for a in range(4) for d in range(4))
                R1[b, c] = sum(dG1[a][b][c][a] - dG1[a][b][a][c] for a in range(4)) + sum(
                    Gm1[a][a][d] * Gm0[d][b][c] + Gm0[a][a][d] * Gm1[d][b][c] - Gm1[a][c][d] * Gm0[d][b][a]
                    - Gm0[a][c][d] * Gm1[d][b][a] for a in range(4) for d in range(4))
        g0p, g0ip, hp = g0.applyfunc(P), g0i.applyfunc(P), h.applyfunc(P)
        Rs0 = sum(g0ip[a, b] * R0[a, b] for a in range(4) for b in range(4))
        hupup = g0ip * hp * g0ip
        Rs1 = sum(g0ip[a, b] * R1[a, b] - hupup[a, b] * R0[a, b] for a in range(4) for b in range(4))
        Ein0 = sp.Matrix(4, 4, lambda b, c: simp(R0[b, c] - g0p[b, c] * Rs0 / 2 - 8 * sp.pi * T0c[b, c]))
        Ein1 = sp.Matrix(4, 4, lambda b, c: simp(R1[b, c] - hp[b, c] * Rs0 / 2 - g0p[b, c] * Rs1 / 2 - 8 * sp.pi * T1c[b, c]))
        out.append(dict(dirac=dirac, Ein0=Ein0, Ein1=Ein1, pt=pt))
    return out

def project(o):
    """projected coefficients at one point (see docstring); returns dict and the list of closure residuals"""
    pt = o['pt']; s, c, sp_, cp = pt
    eip = cp + I * sp_
    cot = c / s
    sn = sp.Matrix([[c, s * (cp - I * sp_)], [s * (cp + I * sp_), -c]])   # sigma.n at the point
    bad = []
    D1 = o['dirac'][+1][1]
    up_ = D1[:2]; lo_ = sn * sp.Matrix(D1[2:])
    A1 = simp(up_[0]); A2 = simp(up_[1] / eip); C1 = simp(lo_[0]); C2 = simp(lo_[1] / eip)
    # closure: A1, C1 proportional to y0; A2, C2 to y1
    bad += [sp.diff(A1, y1), sp.diff(C1, y1), sp.diff(A2, y0), sp.diff(C2, y0)]
    co = dict(A1=simp(A1 / y0), A2=simp(A2 / y1), C1=simp(C1 / y0), C2=simp(C2 / y1))
    E = o['Ein1']
    co['Ett'] = simp(E[0, 0] / y0); bad.append(sp.diff(E[0, 0], y1))
    co['Etr'] = simp(E[0, 1] / y0); bad.append(sp.diff(E[0, 1], y1))
    co['Err'] = simp(E[1, 1] / y0); bad.append(sp.diff(E[1, 1], y1))
    co['Etth'] = simp(E[0, 2] / y1); bad.append(sp.diff(E[0, 2], y0))
    co['Erth'] = simp(E[1, 2] / y1); bad.append(sp.diff(E[1, 2], y0))
    Q = simp(E[3, 3] / s**2); Pq = simp(E[2, 2])
    Eb = simp(sp.diff(Q, y1) / cot); Ea = simp(sp.diff(Q, y0))
    bad += [simp(sp.diff(Pq, y1) + Eb * cot), simp(sp.diff(Pq, y0) - Ea + Lam * Eb)]
    co['Ea'] = Ea; co['Eb'] = Eb
    axial = [E[0, 3], E[1, 3], E[2, 3]]
    return co, bad, axial

FIELDS = MAT + MET
DERIVS = ['', '_t', '_r', '_tt', '_tr', '_rr']

def write_coeffs(CO):
    """write the projected RW-gauge equations as coefficient tables: row = sum_{field, deriv} coeff * d(field)"""
    rows = {}
    for k in ('A1', 'A2', 'C1', 'C2'):
        e = sp.expand(CO[k])
        rows[k + 're'] = sp.expand((e + sp.conjugate(e)) / 2)
        rows[k + 'im'] = sp.expand((e - sp.conjugate(e)) / (2 * I))
    for k in ('Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb'):
        rows[k] = sp.expand(CO[k])
    syms = [sp.Symbol(n + d, real=True) for n in FIELDS for d in DERIVS]
    out = ['"""generated by derivations/dirac_nonradial.py: projected linearised non-radial (polar, M = 0) equations of free',
           '   Dirac stars in the Regge-Wheeler gauge.  row(name) returns a dict {(field, deriv): coefficient array} with',
           '   deriv in "", "_t", "_r", "_tt", "_tr", "_rr";  fields aR1 aI1 aR2 aI2 bR1 bI1 bR2 bI2 H0 H1 H2 K.',
           '   Arguments: r, omega, mu, Lam = L(L+1), g, f, X, alpha (background, dirac_relax conventions)."""',
           'import numpy as np', 'from numpy import pi, sqrt', '', 'ROWS = {}', '']
    summary = []
    for name, e in rows.items():
        terms = []
        rest = e
        for n in FIELDS:
            for d in DERIVS:
                sym = sp.Symbol(n + d, real=True)
                cf = sp.factor(sp.cancel(e.coeff(sym)))
                rest = rest - cf * sym
                if cf != 0:
                    terms.append(((n, d), cf))
        assert sp.simplify(sp.expand(rest)) == 0, (name, rest)
        out.append('def _%s(r, omega, mu, Lam, g_, f_, X_, al_):' % name)
        out.append('    z = np.zeros_like(r)')
        out.append('    return {')
        for (n, d), cf in terms:
            code = sp.pycode(cf).replace('math.', '')
            code = code.replace('omega', 'omega').replace('Lambda', 'Lam')
            out.append('        (%r, %r): (%s) + z,' % (n, d, code))
        out.append('    }')
        out.append("ROWS[%r] = _%s" % (name, name))
        out.append('')
        summary.append('%-6s: %s' % (name, ' '.join(n + d for (n, d), _ in terms)))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'relaxation', 'dirac_nonradial_coeffs.py')
    open(path, 'w').write('\n'.join(out))
    log('coefficients written to relaxation/dirac_nonradial_coeffs.py')
    for sline in summary:
        print('   ' + sline)

if __name__ == '__main__':
    quick = '--quick' in sys.argv
    zero_gauge = {'h0': 0, 'h1': 0, 'Gs': 0}
    # ---------- 1/2. RW-gauge equations, background and closure checks ----------
    h_rw = metric_pert(fld['H0'], fld['H1'], fld['H2'], fld['K'], 0, 0, 0)
    for s2 in (-1, +1):
        o = linearised(h_rw, dpsi(s2), 'RW s2=%+d' % s2)
        co_pts, ok_close, ok_axial = [], True, True
        for k, ok in enumerate(o):
            if s2 == -1 and k == 0:
                bg = [simp(x) for x in ok['dirac'][+1][0]] + [simp(x) for x in ok['dirac'][-1][0]] + list(ok['Ein0'])
                assert all(b == 0 for b in bg), [b for b in bg if b != 0][:2]
                log('1. background Dirac and Einstein equations satisfied: OK')
            co, bad, axial = project(ok)
            co_pts.append(co)
            ok_close &= all(simp(b) == 0 for b in bad)
            ok_axial &= all(simp(a) == 0 for a in axial)
        same = all(simp(co_pts[0][k] - co_pts[1][k]) == 0 for k in co_pts[0])
        log('   s2 = %+d: projection closes: %s, same at both points: %s, axial E_ti/E_thph vanish: %s'
            % (s2, ok_close, same, ok_axial))
        if s2 == -1:
            assert ok_close and same and ok_axial
            CO = co_pts[0]
            # - spinor: its Dirac residual must vanish whenever the + one does (symmetry image)
            D1m = o[0]['dirac'][-1][1]
            pt = o[0]['pt']; s, c, sp_, cp = pt
            emip = cp - I * sp_
            sn = sp.Matrix([[c, s * (cp - I * sp_)], [s * (cp + I * sp_), -c]])
            lo_m = sn * sp.Matrix(D1m[2:])
            # expected: upper = A1 Y chi_dn + s2 A2 e^{-i phi} Y' chi_up, (sigma.n) lower = C1 Y chi_dn + s2 C2 e^{-i phi} Y' chi_up
            exp_res = [simp(D1m[1] - CO['A1'] * y0), simp(D1m[0] - s2 * CO['A2'] * emip * y1),
                       simp(lo_m[1] - CO['C1'] * y0), simp(lo_m[0] - s2 * CO['C2'] * emip * y1)]
            assert all(e == 0 for e in exp_res), exp_res
            log('2. RW gauge, s2 = -1: projections close, are theta/phi independent, no axial parts; - spinor is the '
                'symmetry image of the + spinor: OK')
            write_coeffs(CO)
            if '--coeffs-only' in sys.argv:
                sys.exit(0)
        if quick:
            break
    # ---------- 3. gauge check: a pure-gauge perturbation solves the general-polar equations ----------
    Tz, Rz, Sz = Fn('Tz'), Fn('Rz'), Fn('Sz')
    fld.update({'Tz': Tz, 'Rz': Rz, 'Sz': Sz})
    zeta_up = [Tz * Y, Rz * Y, Sz * Yp / r**2, 0]
    zeta_lo = [sum(g0[m, n] * zeta_up[n] for n in range(4)) for m in range(4)]
    # h = L_zeta g0 = nabla_m zeta_n + nabla_n zeta_m
    def cov_d(m, n):
        return sp.diff(zeta_lo[n], coords[m]) - sum(Gam0[l][m][n] * zeta_lo[l] for l in range(4))
    hg = sp.Matrix(4, 4, lambda m, n: sp.expand(cov_d(m, n) + cov_d(n, m)))
    # read off the polar functions: h_tt = al^2 H0 Y etc.
    H0g = sp.simplify(hg[0, 0] / (al**2 * Y)); H1g = sp.simplify(hg[0, 1] / Y); H2g = sp.simplify(hg[1, 1] / (X**2 * Y))
    h0g = sp.simplify(hg[0, 2] / Yp); h1g = sp.simplify(hg[1, 2] / Yp)
    Gsg = sp.simplify(2 * Sz / r**2)
    Kg = sp.simplify(sp.expand((hg[2, 2] - r**2 * Gsg * sp.diff(Y, th, 2)) / (r**2 * Y)))
    hcheck = metric_pert(H0g, H1g, H2g, Kg, h0g, h1g, Gsg) - hg
    assert all(sp.simplify(reduce_Y(e)) == 0 for e in hcheck), 'L_zeta g is not of the polar form'
    log('3. L_zeta g0 read off as polar (H0, H1, H2, K, h0, h1, G) functions')
    # Kosmann derivative: L_zeta psi = zeta^m nabla_m psi - (1/8) nabla_[a zeta_b] [gamma^a, gamma^b] psi, in frame indices
    found = False
    for sgn in (-1, +1):
        dps = {}
        for m in (+1, -1):
            p = psi0[m]
            covp = [sp.diff(p, coords[k]) + Gs0[k] * p for k in range(4)]
            term1 = sum((zeta_up[k] * covp[k] for k in range(4)), sp.zeros(4, 1))
            # nabla_m zeta_n antisymmetrised, converted to frame: E_a^m E_b^n
            term2 = sp.zeros(4, 1)
            for a in range(4):
                for b in range(a + 1, 4):
                    wab = sum(E0[a, mm] * E0[b, nn] * (cov_d(mm, nn) - cov_d(nn, mm)) for mm in range(4) for nn in range(4)) / 2
                    wab_low = eta[a, a] * eta[b, b] * wab
                    if wab_low != 0:
                        term2 += wab_low * (gam[a] * gam[b] - gam[b] * gam[a]) * p / 4
            dps[m] = term1 + sgn * term2
        o = linearised(metric_pert(H0g, H1g, H2g, Kg, h0g, h1g, Gsg), dps, 'gauge sgn=%+d' % sgn)
        # the gauge perturbation is not in the a_i, b_i parametrisation; check the raw residuals and Einstein components
        ok = True
        for oo in o:
            for m in (+1, -1):
                ok &= all(simp(x) == 0 for x in oo['dirac'][m][1])
            ok &= all(simp(x) == 0 for x in oo['Ein1'])
        log('   Kosmann sign %+d: pure gauge perturbation solves the linearised Dirac + Einstein equations: %s' % (sgn, ok))
        if ok:
            found = True
            break
    assert found, 'gauge check failed'
    log('3. gauge check: OK')
