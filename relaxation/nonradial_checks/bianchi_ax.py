"""phi component of the projected Bianchi identity for the axial tables: vacuum (arbitrary h0, h1) and interior with the
matter made on-shell by integrating the 8 axial Dirac rows."""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from scipy.integrate import solve_ivp
import dirac_axial_coeffs as ac
from bg import get
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
lam = 0.13 + 0.07j; Lam = L * (L + 1)
MAT = ['cR1', 'cI1', 'cR2', 'cI2', 'dR1', 'dI1', 'dR2', 'dI2']
DR = ['A1re', 'A1im', 'A2re', 'A2im', 'C1re', 'C1im', 'C2re', 'C2im']
spec = {'h0': (1.0, 0.11 + 0.05j), 'h1': (0.7, -0.08 + 0.1j)}
def mvals(n, d, r):
    c, k = spec[n]; F = c * np.exp(k * r)
    return {'': F, '_r': k * F, '_rr': k * k * F, '_t': lam * F, '_tr': lam * k * F, '_tt': lam**2 * F}[d]
def identity(r, args, vals):
    E = {}; S = {}
    for name in ('Etph', 'Erph', 'Ethph'):
        out = 0 * r + 0j; sc = 0 * r
        for (n, d), c in ac.ROWS[name](*args).items():
            v = vals(n, d); out = out + c * v; sc = sc + np.abs(c * v)
        E[name], S[name] = out, sc
    X, al = args[6], args[7]
    dal = np.gradient(al, r); dX = np.gradient(X, r)
    I = -(lam * E['Etph'] - al * dal / X**2 * E['Erph']) / al**2 + (np.gradient(E['Erph'], r) - (dX / X + 1 / r) * E['Erph']) / X**2 \
        + 3 * E['Erph'] / (r * X**2) - (Lam - 2) * E['Ethph'] / r**2
    sc = (S['Etph'] + S['Erph'] + S['Ethph']) / r**2
    return np.max(np.abs(I[200:-200]) / (sc[200:-200] + 1e-300))
# vacuum
r = np.linspace(3., 12., 40001); M = 1.
args = (r, 0.9, 1.0, float(Lam), 0 * r, 0 * r, 1 / np.sqrt(1 - 2 * M / r), np.sqrt(1 - 2 * M / r))
print('L %d vacuum (Schwarzschild), arbitrary h0, h1:  max |I_phi| / scale = %.1e' % (L, identity(r, args, lambda n, d: mvals(n, d, r) if n in spec else 0 * r)))
# interior, on-shell matter
sol, dd = get(0.11); w0 = dd['omega']
def bg(r):
    r = np.atleast_1d(r); f, g, m, ph, w = sol.sol(r)
    return (r, w0, 1.0, float(Lam), g, f, 1 / np.sqrt(1 - 2 * m / r), np.exp(ph))
def rhs(r, y):
    a = bg(r); Mr = np.zeros((8, 8), complex); rest = np.zeros(8, complex)
    for i, name in enumerate(DR):
        for (n, d), c in ac.ROWS[name](*a).items():
            c = c[0]
            if n in MAT:
                j = MAT.index(n)
                if d == '_r': Mr[i, j] += c
                elif d == '': rest[i] += c * y[j]
                elif d == '_t': rest[i] += c * lam * y[j]
                else: raise ValueError(d)
            else:
                rest[i] += c * mvals(n, d, r)
    return -np.linalg.solve(Mr, rest)
s = solve_ivp(rhs, (2., 9.), np.array([0.3, -0.2, 0.15, 0.1, -0.25, 0.05, 0.2, -0.1], complex), rtol=1e-11, atol=1e-13, dense_output=True)
r = np.linspace(2., 9., 14001); Yv = s.sol(r); dY = np.array([rhs(ri, Yv[:, k]) for k, ri in enumerate(r)]).T
def vals(n, d):
    if n in MAT:
        j = MAT.index(n); return {'': Yv[j], '_r': dY[j], '_t': lam * Yv[j], '_tr': lam * dY[j]}[d]
    return mvals(n, d, r)
print('L %d interior (A = 0.11), on-shell matter:      max |I_phi| / scale = %.1e' % (L, identity(r, bg(r), vals)))
