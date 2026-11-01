import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, time
from proca_nonradial_modes import classified
from fam import get
A = float(sys.argv[1]); L = int(sys.argv[2]); R = float(sys.argv[3]); N1, N2 = (int(x) for x in sys.argv[4].split(','))
lam_max = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
sol, d = get(A)
t0 = time.time()
cl = classified(sol, L, N1, N2, R=R, lam_max=lam_max)
print('A %.4f omega %.5f M %.5f L %d R %g N %d,%d (%.0fs)' % (A, d['omega'], d['M'], L, R, N1, N2, time.time() - t0))
good = [(z, rr) for z, c, rr in cl if c and rr < 1e-2]
print(' physical, Re > 1e-3:', ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in good if z.real > 1e-3))
print(' physical, lowest |lam|:', ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in sorted(good, key=lambda t: abs(t[0]))[:8] if z.imag >= 0))
cand = [(z, c, rr) for z, c, rr in cl if z.real > 1e-3 and z.imag >= 0]
print(' all Re > 1e-3 candidates (conv, resid):', ' '.join('%.4f%+.4fi[%s,%.0e]' % (z.real, z.imag, 'c' if c else '-', rr) for z, c, rr in cand[:8]))
