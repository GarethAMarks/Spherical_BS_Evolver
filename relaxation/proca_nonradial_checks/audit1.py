import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from proca_nonradial_modes import ProcaModes
from fam import get
A = float(sys.argv[1]); L = int(sys.argv[2]); lam_max = float(sys.argv[3])
sol, d = get(A)
print('A', A, 'omega %.5f' % d['omega'], 'L', L)
for R, N in [(60, 110), (60, 150), (90, 150)]:
    m = ProcaModes(sol, L, N=N, R=R)
    l, V = m.eigen(vectors=True)
    s = (np.abs(l) < lam_max) & (np.abs(l) > 1e-8) & (l.imag >= -1e-12); l = l[s]; V = V[:, s]
    items = []
    for k in np.argsort(-l.real)[:5]:
        rs = m.residuals(l[k], V[:, k]); items.append('%.5f%+.5fi(%s)' % (l[k].real, l[k].imag, ','.join('%.0e' % v for v in rs.values())))
    print(' R %d N %d top Re: %s' % (R, N, ' '.join(items)), flush=True)
