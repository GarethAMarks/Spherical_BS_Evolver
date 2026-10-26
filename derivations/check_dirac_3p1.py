"""
Derivation of the 3+1 (modified-cartoon) evolution equations for the two-spinor Dirac field of DiracStar, for a general
time-dependent spherically symmetric metric with lapse and radial shift, in the (-,+,+,+) signature:

    ds^2 = -alpha^2 dt^2 + g_rr (dr + beta dt)^2 + g_T r^2 dOmega^2,    beta = beta^r,
    (in BSSN variables on the cartoon z-axis: r -> z, g_rr = h_zz / chi, g_T = h_ww / chi).

Tetrad (fixes the local Lorentz freedom): time gauge e_0 = n = (d_t - beta d_r)/alpha, and the symmetric spatial triad,
which in spherical symmetry is the Cartesian-aligned one
    e_i = n_i g_rr^{-1/2} d_r + theta_i (r g_T^{1/2})^{-1} d_theta + phi_i (r g_T^{1/2} sin theta)^{-1} d_phi
(on the z-axis: e_z = g_zz^{-1/2} d_z, e_x = g_ww^{-1/2} d_x, ...). It is an algebraic function of the metric, so no extra
variables are evolved; its time derivative enters through the spin connection.

Matter (as in DiracStar.h): two spinors, (gamma^a nabla_a + mu) psi = 0, {gamma^a, gamma^b} = 2 eta^{ab}, with
    psi_m = ( u(t,r) Omega_{-1,m},  i v(t,r) Omega_{+1,m} ),    m = +-1/2,
u, v complex; u = G / (sqrt(alpha) r), v = F / (sqrt(alpha) r) for the stationary stars (the physical spinor amplitudes
written by DiracStar). On the z-axis u is even and v is odd in z.

The script
  1. solves the Dirac equation for d_t u, d_t v (checking the ansatz is consistent for both spinors, i.e. that every
     component of the Dirac equation is satisfied);
  2. checks current conservation  d_mu (sqrt(-g) J^mu) = 0  identically (arbitrary alpha, beta, g_rr, g_T, including their
     time dependence), with J^a = i psibar gamma^a psi summed over both spinors;
  3. checks the stationary limit reproduces DiracStar's radial equations;
  4. derives the 3+1 sources rho = n^a n^b T_ab, j_r = -gamma_r^a n^b T_ab, S_rr, S_T (= S_thth / r^2), with the time
     derivatives in T eliminated using 1;
  5. checks the compact form, in which the shift and the metric time derivatives (via the ADM equations) combine into
     beta d_r + alpha K / 2 (the triad's evolution enters only through K: no rotation terms in spherical symmetry);
  6. checks the same equations written in the code's BSSN variables on the z-axis:
         d_t u = beta u' + (alpha K/2 - i alpha mu) u + i alpha [ sqrt(chi/h_zz) (v' + L v) + (sqrt(chi/h_zz) + sqrt(chi/h_ww)) v / z ]
         d_t v = beta v' + (alpha K/2 + i alpha mu) v - i alpha [ sqrt(chi/h_zz) (u' + L u) + (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z ]
         L = (h_ww'/h_ww - chi'/chi + alpha'/alpha) / 2,     ' = d_z,   u even, v odd in z.
     At z = 0: v/z -> v', and (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z -> 0 since h_zz - h_ww = O(z^2).
Run:  python3 derivations/check_dirac_3p1.py   (takes about 8 minutes)
"""
import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi', positive=True)
mu = sp.symbols('mu', positive=True)
I = sp.I
coords = [t, r, th, ph]
fn = lambda name: sp.Function(name)(t, r)
al, be, grr, gT = fn('alpha'), fn('beta'), fn('g_rr'), fn('g_T')
uR, uI, vR, vI = fn('u_R'), fn('u_I'), fn('v_R'), fn('v_I')
u, v = uR + I*uI, vR + I*vI

# ---------------- metric, tetrad, spin connection ----------------
g = sp.zeros(4)
g[0, 0] = -al**2 + grr*be**2
g[0, 1] = g[1, 0] = grr*be
g[1, 1] = grr
g[2, 2] = gT * r**2
g[3, 3] = gT * r**2 * sp.sin(th)**2
ginv = sp.zeros(4)
ginv[0, 0] = -1/al**2
ginv[0, 1] = ginv[1, 0] = be/al**2
ginv[1, 1] = 1/grr - be**2/al**2
ginv[2, 2] = 1/(gT*r**2)
ginv[3, 3] = 1/(gT*r**2*sp.sin(th)**2)
eta = sp.diag(-1, 1, 1, 1)

Gam = [[[sum(ginv[a, d]*(sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d]))
             for d in range(4))/2 for c in range(4)] for b in range(4)] for a in range(4)]

nvec = [sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)]
thvec = [sp.cos(th)*sp.cos(ph), sp.cos(th)*sp.sin(ph), -sp.sin(th)]
phvec = [-sp.sin(ph), sp.cos(ph), 0]
E = sp.zeros(4, 4)  # e_a^mu
E[0, 0] = 1/al
E[0, 1] = -be/al
for i in range(3):
    E[i+1, 1] = nvec[i]/sp.sqrt(grr)
    E[i+1, 2] = thvec[i]/(r*sp.sqrt(gT))
    E[i+1, 3] = phvec[i]/(r*sp.sqrt(gT)*sp.sin(th))
Ecov = sp.Matrix(4, 4, lambda a, m: sum(eta[a, b]*g[m, n]*E[b, n] for b in range(4) for n in range(4)))  # e^a_mu
omega = [[[sum(Ecov[a, n]*(sp.diff(E[b, n], coords[m]) + sum(Gam[n][m][l]*E[b, l] for l in range(4))) for n in range(4))
           for b in range(4)] for a in range(4)] for m in range(4)]
omega_low = [[[sum(eta[a, c]*omega[m][c][b] for c in range(4)) for b in range(4)] for a in range(4)] for m in range(4)]

s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]])
sig = [s1, s2, s3]
Z2, I2 = sp.zeros(2), sp.eye(2)
def block(a_, b_, c_, d_):
    Mx = sp.zeros(4); Mx[:2, :2] = a_; Mx[:2, 2:] = b_; Mx[2:, :2] = c_; Mx[2:, 2:] = d_
    return Mx
gD = [block(I2, Z2, Z2, -I2)] + [block(Z2, s, -s, Z2) for s in sig]
gam = [-I*gm for gm in gD]
gam_low = [sum((eta[a, b]*gam[b] for b in range(4)), sp.zeros(4)) for a in range(4)]
beta_mat = I*gam[0]
Gspin = [sum(((omega_low[m][a][b]/8)*(gam[a]*gam[b] - gam[b]*gam[a]) for a in range(4) for b in range(4)), sp.zeros(4))
         for m in range(4)]

Y00 = 1/sp.sqrt(4*sp.pi)
sig_n = sum((nvec[i]*sig[i] for i in range(3)), sp.zeros(2))
chis = [sp.Matrix([1, 0]), sp.Matrix([0, 1])]
def psi_m(chi):
    up = u*Y00*chi
    lo = I*v*(-(sig_n*Y00*chi))
    return sp.Matrix([up[0], up[1], lo[0], lo[1]])
def nabla_frame(psi):
    nc = [sp.diff(psi, coords[m]) + Gspin[m]*psi for m in range(4)]
    return [sum((E[a, m]*nc[m] for m in range(4)), sp.zeros(4, 1)) for a in range(4)]

# ---------------- evaluation at points: functions and derivatives -> real symbols ----------------
base = {al: 'a', be: 'b', grr: 'grr', gT: 'gT', uR: 'uR', uI: 'uI', vR: 'vR', vI: 'vI'}
positive = {'a', 'grr', 'gT'}
SYM = {}
def S(name):
    if name not in SYM:
        SYM[name] = sp.Symbol(name, positive=True) if name in positive else sp.Symbol(name, real=True)
    return SYM[name]
def to_symbols(e):
    e = sp.sympify(e)
    derivs = sorted(e.atoms(sp.Derivative), key=lambda d: -sum(c for _, c in d.variable_count))
    rep = {}
    for d in derivs:
        nm = base[d.expr] + '_' + ''.join(str(vv)*c for vv, c in d.variable_count)
        rep[d] = S(nm)
    e = e.subs(rep)
    return e.subs({f: S(nm) for f, nm in base.items()})
PTS = [(sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(5, 13), sp.Rational(12, 13)),
       (sp.Rational(8, 17), sp.Rational(-15, 17), sp.Rational(-7, 25), sp.Rational(24, 25))]
def at(e, pt):
    e = to_symbols(e)
    return e.subs({sp.sin(th): pt[0], sp.cos(th): pt[1], sp.sin(ph): pt[2], sp.cos(ph): pt[3]})
def simp(e):
    return sp.cancel(sp.together(sp.expand(e)))

# ---------------- 1. Dirac equation -> d_t u, d_t v ----------------
unk = [S('uR_t'), S('uI_t'), S('vR_t'), S('vI_t')]
nab = []
eqs_all = []
for chi in chis:
    psi = psi_m(chi)
    nf = nabla_frame(psi)
    nab.append((psi, nf))
    res = sum((gam[a]*nf[a] for a in range(4)), sp.zeros(4, 1)) + mu*psi
    for pt in PTS:
        for comp in res:
            c = sp.expand(at(comp, pt))
            eqs_all += [sp.re(c), sp.im(c)]
eqs_all = [simp(e) for e in eqs_all]
eqs_nz = [e for e in eqs_all if e != 0]
sol = sp.solve(eqs_nz[:8], unk, dict=True)
assert len(sol) == 1, sol
sol = {k: simp(vv) for k, vv in sol[0].items()}
assert all(simp(e.subs(sol)) == 0 for e in eqs_nz), 'ansatz inconsistent: Dirac equation not satisfied by all components'
print('1. Dirac equation solved for d_t u, d_t v; all components (both spinors, two angular points) consistent: OK')
dtu = simp(sol[S('uR_t')] + I*sol[S('uI_t')])
dtv = simp(sol[S('vR_t')] + I*sol[S('vI_t')])
assert not any(s in (dtu + dtv).free_symbols for s in (S('a_t'), S('b_t'))), 'lapse/shift time derivatives should drop out'

# ---------------- 2. current conservation ----------------
# J^mu = e_a^mu J^a, J^a = i psibar gamma^a psi; check d_mu(sqrt(-g) J^mu) = 0 after using the evolution equations
def dagger(x):
    return x.conjugate().T
sqrt_mg = al*sp.sqrt(grr)*gT*r**2*sp.sin(th)
Jmu = [0]*4
for psi, nf in nab:
    psibar = dagger(psi)*beta_mat
    Ja = [I*(psibar*gam[a]*psi)[0] for a in range(4)]
    for m in range(4):
        Jmu[m] += sum(E[a, m]*Ja[a] for a in range(4))
div = sum(sp.diff(sqrt_mg*Jmu[m], coords[m]) for m in range(4)) / sqrt_mg
cons_ok = True
for pt in PTS:
    c = at(div, pt).subs(sol)
    if simp(sp.expand(c)) != 0:
        cons_ok = False
print('2. current conservation d_mu(sqrt(-g) J^mu) = 0 identically for any time-dependent alpha, beta, g_rr, g_T:', 'OK' if cons_ok else 'FAILED')
assert cons_ok

# ---------------- 3. stationary limit ----------------
w, Fs, Gs, dF, dG = sp.symbols('omega F G dF dG', real=True)
X, al0 = sp.symbols('X alpha0', positive=True)
# u = e^{-i w t} G/(sqrt(alpha) r), v = e^{-i w t} F/(sqrt(alpha) r), static metric, beta = 0, g_rr = X^2, g_T = 1
dal = sp.Symbol('dalpha', real=True); dX = sp.Symbol('dX', real=True)
ghat = Gs/(sp.sqrt(al0)*r); fhat = Fs/(sp.sqrt(al0)*r)
d_ghat = sp.diff(ghat, r) + sp.diff(ghat, Gs)*dG + sp.diff(ghat, al0)*dal
d_fhat = sp.diff(fhat, r) + sp.diff(fhat, Fs)*dF + sp.diff(fhat, al0)*dal
stat = {S('a'): al0, S('a_r'): dal, S('b'): 0, S('b_r'): 0, S('grr'): X**2, S('grr_r'): 2*X*dX, S('gT'): 1, S('gT_r'): 0,
        S('grr_t'): 0, S('gT_t'): 0,
        S('uR'): ghat, S('uI'): 0, S('uR_r'): d_ghat, S('uI_r'): 0, S('vR'): fhat, S('vI'): 0, S('vR_r'): d_fhat, S('vI_r'): 0}
radial = {dG: X*(Gs/r + (mu + w/al0)*Fs), dF: X*(-Fs/r + (mu - w/al0)*Gs)}
res_u = simp((dtu.subs(stat) - (-I*w*ghat)).subs(radial)); res_v = simp((dtv.subs(stat) - (-I*w*fhat)).subs(radial))
ok3 = res_u == 0 and res_v == 0
if not ok3:
    print('   residual u:', sp.factor(res_u)); print('   residual v:', sp.factor(res_v))
    print('   d_t u (stationary):', sp.factor(simp(dtu.subs(stat)))); print('   d_t v (stationary):', sp.factor(simp(dtv.subs(stat))))
print('3. stationary limit reproduces DiracStar\'s radial equations:', 'OK' if ok3 else 'FAILED')
assert ok3

# ---------------- 4. 3+1 sources ----------------
T = sp.zeros(4)
for psi, nf in nab:
    psibar = dagger(psi)*beta_mat
    nfbar = [dagger(x)*beta_mat for x in nf]
    Xab = sp.Matrix(4, 4, lambda a, b: (psibar*gam_low[a]*nf[b])[0] - (nfbar[a]*gam_low[b]*psi)[0])
    T += (Xab + Xab.T)/4
pt0 = PTS[0]
vec = lambda vv: [at(c, pt0) for c in vv]
nn, tt = vec(nvec), vec(thvec)
def ev(e):
    return simp(sp.expand(at(e, pt0).subs(sol)))
rho = ev(T[0, 0])
J = ev(sum(nn[i]*T[0, i+1] for i in range(3)))        # T_{0 rhat}
Sr = ev(sum(nn[i]*nn[j]*T[i+1, j+1] for i in range(3) for j in range(3)))
St = ev(sum(tt[i]*tt[j]*T[i+1, j+1] for i in range(3) for j in range(3)))
assert simp(ev(T[0, 0]) - simp(sp.expand(at(T[0, 0], PTS[1]).subs(sol)))) == 0
uR_, uI_, vR_, vI_ = S('uR'), S('uI'), S('vR'), S('vI')
uR_r, uI_r, vR_r, vI_r = S('uR_r'), S('uI_r'), S('vR_r'), S('vI_r')
a_, grr_, gT_, mu_ = S('a'), S('grr'), S('gT'), mu
nden = simp(ev(sum((dagger(psi)*psi)[0] for psi, _ in nab)))
# compact complex forms (overbar = complex conjugate, ' = d_r)
uc_, vc_ = uR_ + I*uI_, vR_ + I*vI_
ucr_, vcr_ = uR_r + I*uI_r, vR_r + I*vI_r
cj = sp.conjugate
Re = lambda e: sp.expand((e + cj(e))/2)
Im = lambda e: sp.expand((e - cj(e))/(2*I))
sg, sT = sp.sqrt(grr_), sp.sqrt(gT_)
n_c = (uR_**2 + uI_**2 + vR_**2 + vI_**2)/(2*sp.pi)
Srr_c = Re(cj(vc_)*ucr_ - uc_*cj(vcr_))/(2*sp.pi*sg)
Stt_c = -Re(uc_*cj(vc_))/(2*sp.pi*r*sT)
rho_c = mu*(uR_**2 + uI_**2 - vR_**2 - vI_**2)/(2*sp.pi) + Srr_c + 2*Stt_c
J_c = -Im(cj(uc_)*ucr_ + cj(vc_)*vcr_)/(2*sp.pi*sg)
assert all(simp(x - y) == 0 for x, y in [(nden, n_c), (rho, rho_c), (J, J_c), (Sr, Srr_c), (St, Stt_c)])
print('4. 3+1 sources (orthonormal frame; overbar = complex conjugate):')
print('     n       = (|u|^2 + |v|^2) / (2 pi)')
print('     S_rr^   = Re( vbar u\' - u vbar\' ) / (2 pi sqrt(g_rr))')
print('     S_thth^ = -Re( u vbar ) / (2 pi r sqrt(g_T))')
print('     rho     = mu (|u|^2 - |v|^2) / (2 pi) + S_rr^ + 2 S_thth^      (i.e. T^a_a = -mu psibar psi)')
print('     j_r^    = -T_{0 r^} = Im( ubar u\' + vbar v\' ) / (2 pi sqrt(g_rr))')
print('   in code variables: j_z = sqrt(h_zz/chi) j_r^, S_zz = (h_zz/chi) S_rr^, S_ww = (h_ww/chi) S_thth^')

# ---------------- 5. compact form and BSSN variables ----------------
# Claimed form (complex u, v), with K_rr, K_T = K_thth / r^2 and the ADM equations for the metric time derivatives:
#   d_t u = beta u' + (alpha K / 2) u - i alpha mu u + i [ (alpha/sqrt(g_rr)) (v' + (1/r + g_T'/(2 g_T)) v) + alpha v/(r sqrt(g_T)) + alpha' v/(2 sqrt(g_rr)) ]
#   d_t v = beta v' + (alpha K / 2) v + i alpha mu v - i [ (alpha/sqrt(g_rr)) (u' + (1/r + g_T'/(2 g_T)) u) - alpha u/(r sqrt(g_T)) + alpha' u/(2 sqrt(g_rr)) ]
a, b, grr_s, gT_s = S('a'), S('b'), S('grr'), S('gT')
a_r, b_r, grr_r, gT_r = S('a_r'), S('b_r'), S('grr_r'), S('gT_r')
Krr, KT = sp.symbols('K_rr K_T', real=True)
uc, vc = uR_ + I*uI_, vR_ + I*vI_
ucr, vcr = uR_r + I*uI_r, vR_r + I*vI_r
Ktr = Krr/grr_s + 2*KT/gT_s
adm = {S('grr_t'): -2*a*Krr + b*grr_r + 2*grr_s*b_r, S('gT_t'): -2*a*KT + b*gT_r + 2*b*gT_s/r}
dtu_c = b*ucr + a*Ktr/2*uc - I*a*mu*uc + I*((a/sp.sqrt(grr_s))*(vcr + (1/r + gT_r/(2*gT_s))*vc) + a*vc/(r*sp.sqrt(gT_s)) + a_r*vc/(2*sp.sqrt(grr_s)))
dtv_c = b*vcr + a*Ktr/2*vc + I*a*mu*vc - I*((a/sp.sqrt(grr_s))*(ucr + (1/r + gT_r/(2*gT_s))*uc) - a*uc/(r*sp.sqrt(gT_s)) + a_r*uc/(2*sp.sqrt(grr_s)))
ok5 = simp(sp.expand(dtu.subs(adm) - dtu_c)) == 0 and simp(sp.expand(dtv.subs(adm) - dtv_c)) == 0
print('5. compact form (shift and metric time derivatives combine into beta d_r + alpha K / 2):', 'OK' if ok5 else 'FAILED')
assert ok5

# BSSN: g_rr = h_zz/chi, g_T = h_ww/chi, K_ij = (A_ij + h_ij K/3)/chi with A traceless (A_zz/h_zz + 2 A_ww/h_ww = 0)
chi, hzz, hww, Azz, Kb = sp.symbols('chi h_zz h_ww A_zz K', positive=True)
dchi, dhzz, dhww = sp.symbols('dchi dh_zz dh_ww', real=True)
Aww = -hww*Azz/(2*hzz)
bssn = {grr_s: hzz/chi, gT_s: hww/chi, grr_r: dhzz/chi - hzz*dchi/chi**2, gT_r: dhww/chi - hww*dchi/chi**2,
        Krr: (Azz + hzz*Kb/3)/chi, KT: (Aww + hww*Kb/3)/chi}
assert simp(Ktr.subs(bssn) - Kb) == 0
# claimed BSSN form, written with the cartoon variable z = r
z = r
cz, cw = sp.sqrt(chi/hzz), sp.sqrt(chi/hww)
L = (dhww/hww - dchi/chi + a_r/a)/2
dtu_b = b*ucr + (a*Kb/2 - I*a*mu)*uc + I*a*(cz*(vcr + L*vc) + (cz + cw)*vc/z)
dtv_b = b*vcr + (a*Kb/2 + I*a*mu)*vc - I*a*(cz*(ucr + L*uc) + (cz - cw)*uc/z)
ok6 = simp(sp.expand((dtu_c.subs(bssn) - dtu_b))) == 0 and simp(sp.expand((dtv_c.subs(bssn) - dtv_b))) == 0
print('6. BSSN cartoon form verified:')
print('     d_t u = beta d_z u + (alpha K/2 - i alpha mu) u + i alpha [ sqrt(chi/h_zz) (d_z v + L v) + (sqrt(chi/h_zz) + sqrt(chi/h_ww)) v / z ]')
print('     d_t v = beta d_z v + (alpha K/2 + i alpha mu) v - i alpha [ sqrt(chi/h_zz) (d_z u + L u) + (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z ]')
print('     L = (d_z h_ww / h_ww - d_z chi / chi + d_z alpha / alpha) / 2;   u even, v odd in z')
assert ok6
