import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from dirac_axial_modes import AxialModes
from bg_exc import get_exc
from bg import get
A = float(sys.argv[1]); L = int(sys.argv[2]); exc = sys.argv[3] == 'exc'
sol, d = get_exc(A) if exc else get(A)
for R, N in [(80, 110), (80, 150), (100, 150)]:
    m = AxialModes(sol, L, N=N, R=R)
    l, V = m.eigen(vectors=True)
    o = np.argsort(np.abs(l))[:6]
    print('R %d N %d smallest:' % (R, N), ' '.join('%.2e%+.2ei' % (z.real, z.imag) for z in l[o]))
    for k in o[:3]:
        f = m.split(V[:, k]); ph = {n: np.max(np.abs(m.wts[n][0] * f[n])) for n in f}
        nrm = max(ph.values())
        res = m.residuals(l[k], V[:, k])
        print('    %.1e: ' % abs(l[k]) + ' '.join('%s %.0e' % (n, v / nrm) for n, v in ph.items()), ' resid', {kk: '%.0e' % vv for kk, vv in res.items()})
