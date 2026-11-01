import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from dirac_nonradial_modes import NonRadialModes
from bg import get
for A in (0.05, 0.14):
    sol, d = get(A)
    for R, N in [(60, 110), (60, 150), (60, 190), (90, 150), (90, 190)]:
        m = NonRadialModes(sol, 1, N=N, R=R)
        l, V = m.eigen(vectors=True)
        s = (np.abs(l) < 0.6) & (np.abs(l) > 1e-8) & (l.imag >= -1e-12); l = l[s]; V = V[:, s]
        items = []
        for k in np.argsort(-l.real)[:3]:
            rs = m.residuals(l[k], V[:, k]); items.append('%.5f%+.5fi(%.0e)' % (l[k].real, l[k].imag, max(rs['Ett'], rs['Err'], rs['Ea-Eb'])))
        lows = sorted([(l[k], m.residuals(l[k], V[:, k])) for k in range(len(l)) if abs(l[k].real) < 1e-4 and 1e-5 < l[k].imag < 0.12], key=lambda t: t[0].imag)[:3]
        print('A %.2f R %d N %d | top Re: %s | low: %s' % (A, R, N, ' '.join(items), ' '.join('%.5f(%.0e)' % (z.imag, max(r_['Ett'], r_['Err'], r_['Ea-Eb'])) for z, r_ in lows)), flush=True)
