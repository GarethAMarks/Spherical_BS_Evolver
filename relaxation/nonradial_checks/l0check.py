"""End-to-end check of the non-radial coefficient tables: for L = 0 (Lam = 0) with K = H1 = 0 (areal gauge, zero shift)
they must reproduce the radial spectrum of dirac_modes.ModeSolver (independent derivation, validated against evolutions)."""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from scipy.linalg import eig
import dirac_nonradial_coeffs as nc
from dirac_nonradial_modes import NonRadialModes
from dirac_modes import ModeSolver
from dirac_relax import Potential
from bg import get
A = float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 120
sol, d = get(A)
m = NonRadialModes(sol, 0, N=N, regular=False)
m.fields = ['aR1', 'aI1', 'bR1', 'bI1', 'H0', 'H2']; m.idx = {n: k for k, n in enumerate(m.fields)}
rows = ['A1re', 'A1im', 'C1re', 'C1im', 'Etr', 'Err']
Nn = m.N; nf = 6; jR = int(np.argmax(m.r))
Am = np.zeros((nf * Nn, nf * Nn)); Bm = np.zeros_like(Am)
for k, name in enumerate(rows):
    Ar, Br, Cr = m._op(nc.ROWS[name](*m.bgargs))
    assert not Cr.any()
    Am[k * Nn:(k + 1) * Nn] = Ar; Bm[k * Nn:(k + 1) * Nn] = -Br
for k, f in enumerate(['aR1', 'aI1', 'bR1', 'bI1']):
    row = k * Nn + jR; Am[row] = 0; Bm[row] = 0; Am[row, m.idx[f] * Nn + jR] = 1.
row = 5 * Nn + jR; Am[row] = 0; Bm[row] = 0; Am[row, m.idx['H0'] * Nn + jR] = 1.; Am[row, m.idx['H2'] * Nn + jR] = 1.   # alpha X -> 1
lam = eig(Am, Bm, right=False); lam = lam[np.isfinite(lam) & (np.abs(lam) < 0.5) & (np.abs(lam) > 1e-6)]
lam_r = ModeSolver(Potential(lam=0., solitonic=False), sol, N=N).eigenvalues(); lam_r = lam_r[(np.abs(lam_r) < 0.5) & (np.abs(lam_r) > 1e-6)]
def show(l):
    un = l[l.real > 1e-4]; st = np.unique(np.round(np.abs(l[np.abs(l.real) < 1e-4].imag), 6))
    return 'unstable %s  sigma %s' % (np.array2string(un, precision=6), np.array2string(st[:5], precision=6))
print('A %.3f  non-radial tables, L=0:' % A, show(lam))
print('          radial ModeSolver      :', show(lam_r))
