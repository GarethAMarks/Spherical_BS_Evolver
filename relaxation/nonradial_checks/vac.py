import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
import dirac_nonradial_coeffs as nc
r = np.linspace(2.5, 9, 7)
def check(L, M, H, dH, ddH, K, dK, ddK, label):
    X = 1 / np.sqrt(1 - 2 * M / r); al = np.sqrt(1 - 2 * M / r)
    args = (r, 0.9, 1.0, float(L * (L + 1)), 0 * r, 0 * r, X, al)
    vals = {('H0', ''): H, ('H0', '_r'): dH, ('H0', '_rr'): ddH, ('H2', ''): H, ('H2', '_r'): dH, ('H2', '_rr'): ddH,
            ('K', ''): K, ('K', '_r'): dK, ('K', '_rr'): ddK}
    out = []
    for name in ('Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb'):
        row = nc.ROWS[name](*args)
        res = sum(c * vals.get(k, 0 * r) for k, c in row.items())
        sc = sum(np.abs(c * vals.get(k, 0 * r)) for k, c in row.items()) + 1e-300
        out.append('%s %.1e' % (name, np.max(np.abs(res) / sc)))
    print(label, ' | '.join(out))
for L in (2, 3):
    check(L, 0., -2 * r**L, -2 * L * r**(L - 1), -2 * L * (L - 1) * r**(L - 2), -2 * r**L, -2 * L * r**(L - 1), -2 * L * (L - 1) * r**(L - 2), 'flat L=%d' % L)
# Schwarzschild L = 2 static, regular at horizon: H = r (r - 2M); K from the equations (find numerically: try K = r^2 - 2 M^2)
M = 1.
H = r * (r - 2 * M); dH = 2 * r - 2 * M; ddH = 2 + 0 * r
for name, K, dK in [('K=r^2-2M^2', r**2 - 2 * M**2, 2 * r), ('K=r^2', r**2, 2 * r)]:
    check(2, M, H, dH, ddH, K, dK, 2 + 0 * r, 'Schw L=2 ' + name)
