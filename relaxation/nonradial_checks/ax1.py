import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, time
from dirac_axial_modes import AxialModes, classified
from bg import get
from bg_exc import get_exc
A = float(sys.argv[1]); L = int(sys.argv[2]); R = float(sys.argv[3]); exc = len(sys.argv) > 4 and sys.argv[4] == 'exc'
N1, N2 = (110, 150) if len(sys.argv) <= 5 else tuple(int(x) for x in sys.argv[5].split(','))
sol, d = get_exc(A) if exc else get(A)
t0 = time.time()
m = AxialModes(sol, L, N=N2, R=R)
l, V = m.eigen(vectors=True)
z0 = l[np.argsort(np.abs(l))[:3]]
# zero-mode check: smallest |lambda| and its eigenvector's structure
k0 = np.argmin(np.abs(l)); f = m.split(V[:, k0]); phys = {n: m.wts[n][0] * f[n] for n in f}
mx = {n: np.max(np.abs(v)) for n, v in phys.items()}
cl = classified(sol, L, N1, N2, R=R, lam_max=1.0)
good = [(z, rr) for z, c, rr in cl if c and rr < 1e-2]
un = [(z, rr) for z, rr in good if z.real > 5e-4 and z.imag >= 0]
low = sorted([(z, rr) for z, rr in good if abs(z.real) < 1e-4 and z.imag > 1e-6], key=lambda t: t[0].imag)[:5]
cand = sorted([(z, c, rr) for z, c, rr in cl if z.real > 1e-3 and z.imag >= 0], key=lambda t: -t[0].real)[:3]
print('A %.3f omega %.5f %sL %d R %g | smallest |lam| %s | unstable %s | lowest sigma %s | raw %s (%.0fs)' % (A, d['omega'], 'exc ' if exc else '', L, R,
      np.array2string(z0, precision=2), ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in un),
      ' '.join('%.5f' % z.imag for z, rr in low), ' '.join('%.4f%+.4fi[%s,%.0e]' % (z.real, z.imag, 'c' if c else '-', rr) for z, c, rr in cand), time.time() - t0), flush=True)
if L == 1:
    v = V[:, k0]; ff = m.split(v); ph = {n: m.wts[n][0] * ff[n] for n in ff}
    print('   zero-mode content: ' + ' '.join('%s %.1e' % (n, mx[n]) for n in mx))
