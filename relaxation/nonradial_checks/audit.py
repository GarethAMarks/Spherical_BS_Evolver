"""Track eigenvalues with Re > 1e-3 and |lambda| < lam_max across resolutions and domain sizes; a physical mode must
persist with a residual that falls as N grows."""
import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, pickle
from dirac_nonradial_modes import NonRadialModes
from bg import get
from bg_exc import get_exc
A = float(sys.argv[1]); L = int(sys.argv[2]); lam_max = float(sys.argv[3]); exc = len(sys.argv) > 4 and sys.argv[4] == 'exc'
sol, d = get_exc(A) if exc else get(A)
runs = [(60, 110), (60, 150), (60, 190), (80, 150), (80, 190)]
res = {}
for R, N in runs:
    m = NonRadialModes(sol, L, N=N, R=R)
    l, V = m.eigen(vectors=True)
    s = (np.abs(l) < lam_max) & (np.abs(l) > 1e-7); l = l[s]; V = V[:, s]
    out = []
    for k in range(len(l)):
        if l[k].real > 1e-3 or abs(l[k].imag) < 0.25:
            rs = m.residuals(l[k], V[:, k]); out.append((l[k], max(rs['Ett'], rs['Err'], rs['Ea'])))
    res[(R, N)] = out
    un = sorted([o for o in out if o[0].real > 1e-3 and o[0].imag >= 0], key=lambda t: -t[0].real)
    print('R %d N %d: %d candidates Re>1e-3; top:' % (R, N, len(un)), ' '.join('%.4f%+.4fi(%.0e)' % (z.real, z.imag, rr) for z, rr in un[:5]), flush=True)
    lowg = sorted([o for o in out if o[1] < 1e-2 and o[0].imag > 1e-6 and abs(o[0].imag) < 0.25], key=lambda t: t[0].imag)
    print('      low-freq modes with resid<1e-2:', ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in lowg[:6]), flush=True)
pickle.dump(res, open('audit_%s_A%g_L%d.pkl' % ('exc' if exc else 'gs', A, L), 'wb'))
