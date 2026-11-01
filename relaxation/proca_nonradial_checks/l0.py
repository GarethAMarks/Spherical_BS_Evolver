import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from proca_nonradial_modes import ProcaModes
from fam import get
for A in [float(a) for a in sys.argv[1].split(',')]:
    sol, d = get(A)
    res = []
    for N in (90, 120):
        m = ProcaModes(sol, 0, N=N, R=60., fields=['aR0', 'aI0', 'aR1', 'aI1', 'H0', 'H2'],
                       rows=['Ptre', 'Ptim', 'Prre', 'Prim', 'Etr', 'Err'],
                       bc=dict(Ptre='aR0', Ptim='aI0', Prre='aR1', Prim='aI1', Err='H0'))
        l = m.eigen(); l = l[(np.abs(l) < 0.3) & (np.abs(l) > 1e-6)]
        un = np.sort(l[l.real > 1e-4].real)[::-1]
        st = np.unique(np.round(np.abs(l[np.abs(l.real) < 1e-4].imag), 5))
        res.append('N%d: unstable %s sigma %s' % (N, np.array2string(un[:2], precision=5), np.array2string(st[:3], precision=5)))
    print('A %.4f omega %.5f M %.5f | %s' % (A, d['omega'], d['M'], ' | '.join(res)), flush=True)
