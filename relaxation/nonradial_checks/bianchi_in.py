"""Projected Bianchi identities inside the star: matter fields from integrating the 8 Dirac rows (driven by arbitrary
metric fields), so the matter is on-shell; then I_t, I_r, I_theta must vanish."""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from scipy.integrate import solve_ivp
import dirac_nonradial_coeffs as nc
from bg import get
A = float(sys.argv[1]) if len(sys.argv) > 1 else 0.11
L = int(sys.argv[2]) if len(sys.argv) > 2 else 2
zero_metric = len(sys.argv) > 3 and sys.argv[3] == 'nometric'
lam = 0.13 + 0.07j; Lam = L * (L + 1)
sol, d = get(A); w0 = d['omega']
MAT = ['aR1', 'aI1', 'aR2', 'aI2', 'bR1', 'bI1', 'bR2', 'bI2']
DR = ['A1re', 'A1im', 'A2re', 'A2im', 'C1re', 'C1im', 'C2re', 'C2im']
spec = {'H0': (1.0, 0.11 + 0.05j), 'H1': (0.7, -0.08 + 0.1j), 'H2': (-0.4, 0.07 - 0.03j), 'K': (0.9, -0.05 + 0.06j)}
if zero_metric:
    spec = {k: (0., 0.) for k in spec}
def mvals(name, d, r):
    c, k = spec[name]; F = c * np.exp(k * r)
    return {'': F, '_r': k * F, '_rr': k * k * F, '_t': lam * F, '_tr': lam * k * F, '_tt': lam**2 * F}[d]
def bg(r):
    r = np.atleast_1d(r)
    f, g, m, ph, w = sol.sol(r)
    return (r, w0, 1.0, float(Lam), g, f, 1 / np.sqrt(1 - 2 * m / r), np.exp(ph))
def rhs(r, y):
    args = bg(r)
    Mr = np.zeros((8, 8), complex); rest = np.zeros(8, complex)
    for i, name in enumerate(DR):
        for (n, dd), c in nc.ROWS[name](*args).items():
            c = c[0]
            if n in MAT:
                j = MAT.index(n)
                if dd == '_r': Mr[i, j] += c
                elif dd == '': rest[i] += c * y[j]
                elif dd == '_t': rest[i] += c * lam * y[j]
                else: raise ValueError(dd)
            else:
                rest[i] += c * mvals(n, dd, r)
    return -np.linalg.solve(Mr, rest)
r0, r1 = 2.0, 9.0
y0 = np.array([0.3, -0.2, 0.15, 0.1, -0.25, 0.05, 0.2, -0.1], complex)
s = solve_ivp(rhs, (r0, r1), y0, rtol=1e-11, atol=1e-13, dense_output=True)
r = np.linspace(r0, r1, 14001)
Y = s.sol(r); dY = np.array([rhs(ri, Y[:, k]) for k, ri in enumerate(r)]).T
args = bg(r)
def vals(n, dd):
    if n in MAT:
        j = MAT.index(n)
        return {'': Y[j], '_r': dY[j], '_t': lam * Y[j], '_tr': lam * dY[j]}[dd]
    return mvals(n, dd, r)
E, S = {}, {}
for name in ('Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb') + tuple(DR):
    out = 0 * r + 0j; sc = 0 * r
    for (n, dd), c in nc.ROWS[name](*args).items():
        v = vals(n, dd); out = out + c * v; sc = sc + np.abs(c * v)
    E[name], S[name] = out, sc
print('Dirac rows residual (should be ~0):', '%.1e' % max(np.max(np.abs(E[n]) / (S[n] + 1e-300)) for n in DR))
X = args[6]; al = args[7]
dal = np.gradient(al, r); dX = np.gradient(X, r)
G_rtt = al * dal / X**2; G_ttr = dal / al; G_rrr = dX / X
D = lambda f: np.gradient(f, r)
It = -(lam * E['Ett'] - 2 * G_rtt * E['Etr']) / al**2 + (D(E['Etr']) - (G_rrr + G_ttr) * E['Etr']) / X**2 \
     + 2 * E['Etr'] / (r * X**2) - Lam * E['Etth'] / r**2
Ir = -(lam * E['Etr'] - G_rtt * E['Err'] - G_ttr * E['Ett']) / al**2 + (D(E['Err']) - 2 * G_rrr * E['Err']) / X**2 \
     - Lam * E['Erth'] / r**2 + 2 * E['Err'] / (r * X**2) - (2 * E['Ea'] - Lam * E['Eb']) / r**3
Ith = -(lam * E['Etth'] - G_rtt * E['Erth']) / al**2 + (D(E['Erth']) - G_rrr * E['Erth'] - E['Erth'] / r) / X**2 \
      + 3 * E['Erth'] / (r * X**2) + (E['Ea'] + (1 - Lam) * E['Eb']) / r**2
inner = slice(200, -200)
for nm, I_, sc in (('I_t', It, S['Ett'] + S['Etr']), ('I_r', Ir, S['Err'] + S['Etr']), ('I_th', Ith, S['Ea'] + S['Erth'])):
    print('%-5s max |I| / scale = %.2e   (|I| max %.2e)' % (nm, np.max(np.abs(I_[inner]) / (sc[inner] / r[inner]**2 + 1e-300)), np.max(np.abs(I_[inner]))))
