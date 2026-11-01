import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from dirac_nonradial_modes import NonRadialModes
from bg_exc import get_exc
L = int(sys.argv[2]); R = float(sys.argv[3])
for A in [float(a) for a in sys.argv[1].split(',')]:
    sol, d = get_exc(A, R=R)
    line = []
    for N in (130, 170):
        m = NonRadialModes(sol, L, N=N, R=R)
        l, V = m.eigen(vectors=True)
        s = (np.abs(l) < 0.02) & (np.abs(l) > 1e-8) & (l.imag >= -1e-12); l = l[s]; V = V[:, s]
        items = []
        for k in np.argsort(np.abs(l)):
            rs = m.residuals(l[k], V[:, k]); items.append('%.5f%+.5fi(%.0e)' % (l[k].real, l[k].imag, max(rs['Ett'], rs['Err'], rs['Ea'])))
        line.append('N%d: %s' % (N, ' '.join(items)))
    print('A %.4f omega %.5f | %s' % (A, d['omega'], ' | '.join(line)), flush=True)
