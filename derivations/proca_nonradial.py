"""
Linearised NON-RADIAL (polar, axisymmetric M = 0, multipole L) perturbation equations of free spherical Proca stars,
signature (-,+,+,+), G = 1, D = 4, in the conventions of ProcaStar:

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,  alpha = exp(phi0),
    A = e^{-i w t} ( f dt + i g dr ),  E = f' - w g,  g = -w E / (alpha^2 mu^2)   (free field),
    field equation  nabla_a F^{ab} = mu^2 A^b,
    T_ab = (1/2)(F_ac Fbar_b^c + Fbar_ac F_b^c) - (1/4) g_ab F_cd Fbar^cd
           + (mu^2/2)(A_a Abar_b + Abar_a A_b) - (mu^2/2) g_ab A_c Abar^c        (normalisation checked against ProcaStar).

Perturbation (polar parity (-1)^L, M = 0; Y = P_L(cos theta), Y' = dY/dtheta):
    dA_t = e^{-i w t} a0 Y,   dA_r = e^{-i w t} a1 Y,   dA_theta = e^{-i w t} a2 Y',   dA_phi = 0,
with complex a_k = aRk + i aIk, real functions of (t, r). Metric perturbation as in dirac_nonradial.py:
    h_tt = alpha^2 H0 Y, h_tr = H1 Y, h_rr = X^2 H2 Y, h_tA = h0 D_A Y, h_rA = h1 D_A Y, h_AB = r^2 (K Y Omega_AB + Gs D_A D_B Y),
with the Regge-Wheeler gauge h0 = h1 = Gs = 0 for the mode solver.

Projections (general L through Y'' = -cot Y' - Lam Y, evaluated at generic rational points with Y -> y0, Y' -> y1):
    Proca  P^b = nabla_a F^{ab} - mu^2 A^b:  P^t = Pt Y, P^r = Pr Y, P^theta = Pth Y' / r^2 (Pth defined as r^2 x coefficient),
           times e^{i w t}, split into real and imaginary parts;
    Einstein E_mn = dG_mn - 8 pi dT_mn: as in dirac_nonradial.py (Ett, Etr, Err, Etth, Erth, Ea, Eb).
Checks: 1. background (ProcaStar's equations and stress-tensor normalisation); 2. closure at two points, P^phi and the
axial Einstein components vanish; 3. pure gauge: dA = L_zeta A0, h = L_zeta g0 for a general polar zeta (T Y, R Y, S Y'/r^2)
solves all projected equations (in the general polar form, including h0, h1, Gs).
Output: relaxation/proca_nonradial_coeffs.py (RW gauge), same layout as dirac_nonradial_coeffs.py, arguments
(r, omega, mu, Lam, f, E, X, alpha).
Run: python3 derivations/proca_nonradial.py
"""
import os, sys, time
import sympy as sp

T0 = time.time()
def log(*a):
    print('[%5.0fs]' % (time.time() - T0), *a, flush=True)

t, r, th, ph = sp.symbols('t r theta phi_', real=True)
coords = [t, r, th, ph]
w, mu, Lam = sp.symbols('omega mu Lambda', real=True)
eps = sp.Symbol('epsilon', real=True)
I = sp.I
phi0 = sp.Function('phi0', real=True)(r)
al = sp.exp(phi0)
X = sp.Function('X', real=True)(r)
fb = sp.Function('f', real=True)(r)
Eb_ = sp.Function('E', real=True)(r)
gb = -w * Eb_ / (al**2 * mu**2)
Y = sp.Function('Y', real=True)(th)
Yp = sp.diff(Y, th)
s_, c_ = sp.sin(th), sp.cos(th)

def Fn(name):
    return sp.Function(name, real=True)(t, r)
MAT = ['aR0', 'aI0', 'aR1', 'aI1', 'aR2', 'aI2']
MET = ['H0', 'H1', 'H2', 'K']
GAU = ['h0', 'h1', 'Gs', 'Tz', 'Rz', 'Sz']
fld = {n: Fn(n) for n in MAT + MET + GAU}
a0 = fld['aR0'] + I * fld['aI0']; a1 = fld['aR1'] + I * fld['aI1']; a2 = fld['aR2'] + I * fld['aI2']

g0 = sp.diag(-al**2, X**2, r**2, r**2 * s_**2)
g0i = sp.diag(-1 / al**2, 1 / X**2, 1 / r**2, 1 / (r**2 * s_**2))
Gam0 = [[[sp.expand(sum(g0i[a, d] * (sp.diff(g0[d, b], coords[c]) + sp.diff(g0[d, c], coords[b]) - sp.diff(g0[b, c], coords[d]))
                       for d in range(4)) / 2) for c in range(4)] for b in range(4)] for a in range(4)]
ph_t = sp.exp(-I * w * t)
A0 = [ph_t * fb, ph_t * I * gb, 0, 0]

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

# ---------------- background rules (ProcaStar, free field) ----------------
def bg_rules():
    gg = gb
    E_term = Eb_**2 / (al**2 * X**2)
    two_rho = E_term + mu**2 * (fb**2 / al**2 + gg**2 / X**2)
    two_pr = -E_term + mu**2 * (gg**2 / X**2 + fb**2 / al**2)
    dX = X * (2 * sp.pi * r * X**2 * two_rho - (X**2 - 1) / (2 * r))
    dphi = (X**2 - 1) / (2 * r) + 2 * sp.pi * r * X**2 * two_pr
    df = Eb_ + w * gg
    dE = mu**2 * fb * X**2 - 2 * Eb_ / r + Eb_ * (dphi + dX / X)
    first = {sp.Derivative(phi0, r): dphi, sp.Derivative(X, r): dX, sp.Derivative(fb, r): df, sp.Derivative(Eb_, r): dE}
    rules = [dict(first)]
    cur = {phi0: dphi, X: dX, fb: df, Eb_: dE}
    for n in range(2, 5):
        nxt = {F: sp.diff(e, r).subs(first) for F, e in cur.items()}
        rules.append({sp.Derivative(F, (r, n)): e for F, e in nxt.items()})
        cur = nxt
    return rules[::-1], two_rho / 2
BG_RULES, RHO_PS = bg_rules()
fS, ES, XS, aS = sp.symbols('f_ E_ X_ al_', real=True)

PTS = [(sp.Rational(3, 5), sp.Rational(4, 5)), (sp.Rational(8, 17), sp.Rational(-15, 17))]
y0, y1 = sp.symbols('y0 y1', real=True)

def reduce_Y(e):
    Y2 = -sp.cot(th) * Yp - Lam * Y
    Y3 = sp.diff(Y2, th).subs(sp.Derivative(Y, (th, 2)), Y2)
    Y4 = sp.diff(Y3, th).subs(sp.Derivative(Y, (th, 2)), Y2)
    return e.subs(sp.Derivative(Y, (th, 4)), Y4).subs(sp.Derivative(Y, (th, 3)), Y3).subs(sp.Derivative(Y, (th, 2)), Y2)

def at_point(e, pt):
    e = reduce_Y(sp.sympify(e))
    for rule in BG_RULES:
        e = e.subs(rule)
    e = e.replace(sp.tan, lambda x: sp.sin(x) / sp.cos(x)).replace(sp.cot, lambda x: sp.cos(x) / sp.sin(x))
    e = e.subs({sp.sin(th): pt[0], sp.cos(th): pt[1]})
    e = e.subs({sp.Derivative(Y, th): y1}).subs(Y, y0)
    assert not e.has(th), e.atoms(sp.Function)
    rep = {}
    for n, F in fld.items():
        for k, dv in [('_ttt', ((t, 3),)), ('_ttr', ((t, 2), r)), ('_trr', (t, (r, 2))), ('_rrr', ((r, 3),)),
                      ('_tt', ((t, 2),)), ('_tr', (t, r)), ('_rr', ((r, 2),))]:
            rep[sp.Derivative(F, *dv)] = sp.Symbol(n + k, real=True)
    e = e.subs(rep)
    e = e.subs({sp.Derivative(F, t): sp.Symbol(n + '_t', real=True) for n, F in fld.items()})
    e = e.subs({sp.Derivative(F, r): sp.Symbol(n + '_r', real=True) for n, F in fld.items()})
    e = e.subs({F: sp.Symbol(n, real=True) for n, F in fld.items()})
    e = e.subs(sp.exp(phi0), aS).subs(sp.exp(-phi0), 1 / aS).subs(phi0, sp.log(aS))
    e = e.subs({fb: fS, Eb_: ES, X: XS})
    assert not e.has(sp.Derivative), list(e.atoms(sp.Derivative))[:3]
    return e

def simp(e):
    return sp.cancel(sp.together(sp.expand(e)))

def linearised(h, dA, label):
    """Proca residual P^b and Einstein components (zeroth and first order) at the points"""
    log(label, ': geometry')
    g = g0 + eps * h
    gi = g0i - eps * g0i * h * g0i                      # first order inverse
    sqrtg0 = al * X * r**2 * s_
    sqrtg = sqrtg0 * (1 + eps * sum(g0i[a, a] * h[a, a] for a in range(4)) / 2)
    A = [A0[a] + eps * dA[a] for a in range(4)]
    Ab = [sp.conjugate(x) for x in A]
    F = sp.Matrix(4, 4, lambda a, b: sp.diff(A[b], coords[a]) - sp.diff(A[a], coords[b]))
    Fb = sp.Matrix(4, 4, lambda a, b: sp.diff(Ab[b], coords[a]) - sp.diff(Ab[a], coords[b]))
    Fup = gi * F * gi.T
    Fbup = gi * Fb * gi.T
    Aup = gi * sp.Matrix(A); Abup = gi * sp.Matrix(Ab)
    def lin0(e):
        return e.subs(eps, 0)
    def lin1(e):
        return sp.diff(e, eps).subs(eps, 0)
    # Proca: P^b = (1/sqrtg) d_a (sqrtg F^{ab}) - mu^2 A^b
    log(label, ': Proca')
    P = []
    for b in range(4):
        e = sum(sp.diff(sqrtg * Fup[a, b], coords[a]) for a in range(4)) / sqrtg - mu**2 * Aup[b]
        P.append((lin0(e), lin1(e)))
    # stress tensor (lower indices)
    log(label, ': stress tensor')
    FF = sum(F[a, b] * Fbup[a, b] for a in range(4) for b in range(4))
    AA = sum(A[a] * Abup[a] for a in range(4))
    T = sp.zeros(4)
    for a in range(4):
        for b in range(a, 4):
            Fa_Fb = sum(F[a, c] * Fb[b, d] * gi[c, d] for c in range(4) for d in range(4))
            Fba_F = sum(Fb[a, c] * F[b, d] * gi[c, d] for c in range(4) for d in range(4))
            e = (Fa_Fb + Fba_F) / 2 - g[a, b] * FF / 4 + mu**2 * (A[a] * Ab[b] + Ab[a] * A[b]) / 2 - mu**2 * g[a, b] * AA / 2
            T[a, b] = T[b, a] = e
    T0 = T.applyfunc(lin0); T1 = T.applyfunc(lin1)
    # linearised Einstein tensor
    log(label, ': Einstein')
    hup = g0i * h
    Gam1 = [[[sp.expand(-sum(hup[a, e] * Gam0[e][b][c] for e in range(4))
                        + sum(g0i[a, d] * (sp.diff(h[d, c], coords[b]) + sp.diff(h[d, b], coords[c]) - sp.diff(h[b, c], coords[d]))
                              for d in range(4)) / 2) for c in range(4)] for b in range(4)] for a in range(4)]
    out = []
    for ip, pt in enumerate(PTS):
        Pp = lambda e: at_point(e, pt)
        log(label, ': point', ip)
        Gm0 = [[[Pp(Gam0[a][b][c]) for c in range(4)] for b in range(4)] for a in range(4)]
        Gm1 = [[[Pp(Gam1[a][b][c]) for c in range(4)] for b in range(4)] for a in range(4)]
        dG0 = [[[[Pp(sp.diff(Gam0[a][b][c], coords[e])) for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
        dG1 = [[[[Pp(sp.diff(Gam1[a][b][c], coords[e])) for e in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
        R0 = sp.zeros(4); R1 = sp.zeros(4)
        for b in range(4):
            for c in range(4):
                R0[b, c] = sum(dG0[a][b][c][a] - dG0[a][b][a][c] for a in range(4)) + sum(
                    Gm0[a][a][d] * Gm0[d][b][c] - Gm0[a][c][d] * Gm0[d][b][a] for a in range(4) for d in range(4))
                R1[b, c] = sum(dG1[a][b][c][a] - dG1[a][b][a][c] for a in range(4)) + sum(
                    Gm1[a][a][d] * Gm0[d][b][c] + Gm0[a][a][d] * Gm1[d][b][c] - Gm1[a][c][d] * Gm0[d][b][a]
                    - Gm0[a][c][d] * Gm1[d][b][a] for a in range(4) for d in range(4))
        g0p, g0ip, hp = g0.applyfunc(Pp), g0i.applyfunc(Pp), h.applyfunc(Pp)
        Rs0 = sum(g0ip[a, b] * R0[a, b] for a in range(4) for b in range(4))
        hupup = g0ip * hp * g0ip
        Rs1 = sum(g0ip[a, b] * R1[a, b] - hupup[a, b] * R0[a, b] for a in range(4) for b in range(4))
        T0p = T0.applyfunc(Pp); T1p = T1.applyfunc(Pp)
        Ein0 = sp.Matrix(4, 4, lambda b, c: simp(R0[b, c] - g0p[b, c] * Rs0 / 2 - 8 * sp.pi * T0p[b, c]))
        Ein1 = sp.Matrix(4, 4, lambda b, c: simp(R1[b, c] - hp[b, c] * Rs0 / 2 - g0p[b, c] * Rs1 / 2 - 8 * sp.pi * T1p[b, c]))
        ephase = sp.exp(I * w * t)
        Pp0 = [simp(Pp(p[0]) * ephase) for p in P]
        Pp1 = [simp(Pp(p[1]) * ephase) for p in P]
        out.append(dict(P0=Pp0, P1=Pp1, Ein0=Ein0, Ein1=Ein1, T0=T0p, pt=pt))
    return out

def project(o):
    s, c = o['pt']; cot = c / s
    bad = []
    P1 = o['P1']
    co = {}
    co['Pt'] = simp(P1[0] / y0); bad.append(sp.diff(P1[0], y1))
    co['Pr'] = simp(P1[1] / y0); bad.append(sp.diff(P1[1], y1))
    co['Pth'] = simp(r**2 * P1[2] / y1); bad.append(sp.diff(P1[2], y0))
    bad.append(P1[3])
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
    bad += [E[0, 3], E[1, 3], E[2, 3]]
    return co, bad

FIELDS = MAT + MET
DERIVS = ['', '_t', '_r', '_tt', '_tr', '_rr']

def write_coeffs(CO):
    rows = {}
    for k in ('Pt', 'Pr', 'Pth'):
        e = sp.expand(CO[k])
        rows[k + 're'] = sp.expand((e + sp.conjugate(e)) / 2)
        rows[k + 'im'] = sp.expand((e - sp.conjugate(e)) / (2 * I))
    for k in ('Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb'):
        rows[k] = sp.expand(CO[k])
    out = ['"""generated by derivations/proca_nonradial.py: projected linearised non-radial (polar, M = 0) equations of free',
           '   Proca stars in the Regge-Wheeler gauge.  ROWS[name](r, omega, mu, Lam, f, E, X, alpha) returns a dict',
           '   {(field, deriv): coefficient array}, deriv in "", "_t", "_r", "_tt", "_tr", "_rr";',
           '   fields aR0 aI0 aR1 aI1 aR2 aI2 H0 H1 H2 K."""',
           'import numpy as np', 'from numpy import pi, sqrt', '', 'ROWS = {}', '']
    summary = []
    for name, e in rows.items():
        terms = []; rest = e
        for n in FIELDS:
            for d in DERIVS:
                sym = sp.Symbol(n + d, real=True)
                cf = sp.factor(sp.cancel(e.coeff(sym)))
                rest = rest - cf * sym
                if cf != 0:
                    terms.append(((n, d), cf))
        assert sp.simplify(sp.expand(rest)) == 0, (name, rest)
        out.append('def _%s(r, omega, mu, Lam, f_, E_, X_, al_):' % name)
        out.append('    z = np.zeros_like(r)')
        out.append('    return {')
        for (n, d), cf in terms:
            out.append('        (%r, %r): (%s) + z,' % (n, d, sp.pycode(cf).replace('math.', '').replace('Lambda', 'Lam')))
        out.append('    }')
        out.append('ROWS[%r] = _%s' % (name, name))
        out.append('')
        summary.append('%-6s: %s' % (name, ' '.join(n + d for (n, d), _ in terms)))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'relaxation', 'proca_nonradial_coeffs.py')
    open(path, 'w').write('\n'.join(out))
    log('coefficients written to relaxation/proca_nonradial_coeffs.py')
    for sl in summary:
        print('   ' + sl)

if __name__ == '__main__':
    # ---------- 1/2. RW gauge ----------
    h = metric_pert(fld['H0'], fld['H1'], fld['H2'], fld['K'], 0, 0, 0)
    dA = [ph_t * a0 * Y, ph_t * a1 * Y, ph_t * a2 * Yp, 0]
    o = linearised(h, dA, 'RW')
    # background: stress-tensor normalisation and equations
    rho_pt = at_point(RHO_PS, PTS[0])
    T00 = o[0]['T0'][0, 0]
    ratio = simp(T00 / (aS**2 * rho_pt))
    assert ratio == 1, ratio
    for oo in o:
        assert all(simp(x) == 0 for x in oo['P0']), oo['P0']
        assert all(x == 0 for x in oo['Ein0']), [x for x in oo['Ein0'] if x != 0][:1]
    log('1. stress tensor matches ProcaStar (rho), background Proca and Einstein equations satisfied: OK')
    cos_ = []
    for oo in o:
        co, bad = project(oo)
        assert all(simp(b) == 0 for b in bad), [simp(b) for b in bad if simp(b) != 0][:2]
        cos_.append(co)
    assert all(simp(cos_[0][k] - cos_[1][k]) == 0 for k in cos_[0])
    log('2. projections close at both points, P^phi and axial Einstein components vanish, theta-independent: OK')
    write_coeffs(cos_[0])
    if '--no-gauge' in sys.argv:
        sys.exit(0)
    # ---------- 3. gauge check ----------
    Tz, Rz, Sz = fld['Tz'], fld['Rz'], fld['Sz']
    zeta_up = [Tz * Y, Rz * Y, Sz * Yp / r**2, 0]
    zeta_lo = [sum(g0[m, n] * zeta_up[n] for n in range(4)) for m in range(4)]
    def cov_d(m, n):
        return sp.diff(zeta_lo[n], coords[m]) - sum(Gam0[l][m][n] * zeta_lo[l] for l in range(4))
    hg = sp.Matrix(4, 4, lambda m, n: sp.expand(cov_d(m, n) + cov_d(n, m)))
    H0g = sp.simplify(hg[0, 0] / (al**2 * Y)); H1g = sp.simplify(hg[0, 1] / Y); H2g = sp.simplify(hg[1, 1] / (X**2 * Y))
    h0g = sp.simplify(hg[0, 2] / Yp); h1g = sp.simplify(hg[1, 2] / Yp); Gsg = sp.simplify(2 * Sz / r**2)
    Kg = sp.simplify(sp.expand((hg[2, 2] - r**2 * Gsg * sp.diff(Y, th, 2)) / (r**2 * Y)))
    assert all(sp.simplify(reduce_Y(e)) == 0 for e in metric_pert(H0g, H1g, H2g, Kg, h0g, h1g, Gsg) - hg)
    dAg = [sum(zeta_up[c] * sp.diff(A0[a], coords[c]) + A0[c] * sp.diff(zeta_up[c], coords[a]) for c in range(4)) for a in range(4)]
    og = linearised(metric_pert(H0g, H1g, H2g, Kg, h0g, h1g, Gsg), dAg, 'gauge')
    ok = all(simp(x) == 0 for oo in og for x in list(oo['P1']) + list(oo['Ein1']))
    log('3. pure gauge (h = L_zeta g0, dA = L_zeta A0) solves the linearised Proca + Einstein equations:', ok)
    assert ok
