"""Projected linearised Bianchi identities, I_t, I_r, I_theta, evaluated with the coefficient tables on a Schwarzschild
exterior (no matter) for arbitrary metric perturbations ~ e^{lam t}. All three must vanish identically."""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
import dirac_nonradial_coeffs as nc
M, L, lam = 1.0, 2, 0.3 + 0.2j
Lam = L * (L + 1)
r = np.linspace(3., 12., 40001)
al = np.sqrt(1 - 2 * M / r); X = 1 / al
dal = np.gradient(al, r); dX = np.gradient(X, r)
# arbitrary smooth fields with analytic derivatives: F = c * exp(k r)
spec = {'H0': (1.0, 0.31 + 0.1j), 'H1': (0.7, -0.22 + 0.3j), 'H2': (-0.4, 0.17 - 0.05j), 'K': (0.9, -0.13 + 0.2j)}
def vals(name, d):
    c, k = spec[name]
    F = c * np.exp(k * r)
    if d == '': return F
    if d == '_r': return k * F
    if d == '_rr': return k * k * F
    if d == '_t': return lam * F
    if d == '_tr': return lam * k * F
    if d == '_tt': return lam**2 * F
def row(name):
    co = nc.ROWS[name](r, 0.9, 1.0, float(Lam), 0 * r, 0 * r, X, al)
    out = 0 * r + 0j; sc = 0 * r
    for (n, d), c in co.items():
        if n in spec:
            out = out + c * vals(n, d); sc = sc + np.abs(c * vals(n, d))
    return out, sc
E = {}; S = {}
for n in ('Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb'):
    E[n], S[n] = row(n)
G_rtt = al * dal / X**2; G_ttr = dal / al; G_rrr = dX / X
D = lambda f: np.gradient(f, r)
It = -(lam * E['Ett'] - 2 * G_rtt * E['Etr']) / al**2 + (D(E['Etr']) - (G_rrr + G_ttr) * E['Etr']) / X**2 \
     + 2 * E['Etr'] / (r * X**2) - Lam * E['Etth'] / r**2
Ir = -(lam * E['Etr'] - G_rtt * E['Err'] - G_ttr * E['Ett']) / al**2 + (D(E['Err']) - 2 * G_rrr * E['Err']) / X**2 \
     - Lam * E['Erth'] / r**2 + 2 * E['Err'] / (r * X**2) - (2 * E['Ea'] - Lam * E['Eb']) / r**3
Ith = -(lam * E['Etth'] - G_rtt * E['Erth']) / al**2 + (D(E['Erth']) - G_rrr * E['Erth'] - E['Erth'] / r) / X**2 \
      + 3 * E['Erth'] / (r * X**2) + (E['Ea'] + (1 - Lam) * E['Eb']) / r**2
inner = slice(100, -100)
for nm, I_, scl in (('I_t', It, S['Ett'] + S['Etr']), ('I_r', Ir, S['Err'] + S['Etr']), ('I_th', Ith, S['Ea'] + S['Erth'])):
    print('%-5s max |I| / scale = %.2e' % (nm, np.max(np.abs(I_[inner]) / (scl[inner] / r[inner]**2 + 1e-300))))
