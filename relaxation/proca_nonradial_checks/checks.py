import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from proca_nonradial_modes import ProcaModes
from fam import get
jobs = [(0.01, 1, [(90, 150), (90, 190), (120, 190)]), (0.03, 1, [(90, 150)]), (0.02, 2, [(90, 150)]),
        (0.003, 2, [(150, 150), (150, 190)]), (0.039, 2, [(60, 150)]), (0.0395, 2, [(60, 150)])]
for A, L, runs in jobs:
    sol, d = get(A)
    for R, N in runs:
        m = ProcaModes(sol, L, N=N, R=R)
        l, V = m.eigen(vectors=True)
        s = (np.abs(l) < 0.3) & (np.abs(l) > 1e-8) & (l.imag >= -1e-12); l = l[s]; V = V[:, s]
        good = []
        for k in range(len(l)):
            if l[k].real > 5e-4 or abs(l[k]) < 0.005:
                rs = m.residuals(l[k], V[:, k]); r_ = max(rs.values())
                if r_ < 5e-2: good.append((l[k], r_))
        good.sort(key=lambda t: -t[0].real)
        print('A %.4f omega %.5f L %d R %d N %d: %s' % (A, d['omega'], L, R, N, ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in good[:4])), flush=True)
