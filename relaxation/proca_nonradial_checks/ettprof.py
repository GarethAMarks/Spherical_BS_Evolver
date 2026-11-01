import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np
from proca_nonradial_modes import ProcaModes
import proca_nonradial_coeffs as pc
from fam import get
A = 0.02; sol, d = get(A)
for R, N in [(60, 150), (90, 150)]:
    m = ProcaModes(sol, 1, N=N, R=R)
    l, V = m.eigen(vectors=True)
    k = np.argmin(np.abs(l - (0.00282 + 0.02653j))); z = l[k]; v = V[:len(m.fields) * m.N, k]
    Ar, Br, Cr = m._op(pc.ROWS['Ett'](*m.bgargs)); res = (Ar + z * Br + z**2 * Cr) @ v
    sc = (np.abs(Ar) + abs(z) * np.abs(Br)) @ np.abs(v)
    f = m.split(v); phys = {n: m.wts[n][0] * f[n] for n in f}
    o = np.argsort(m.r)
    print('R %d N %d lam %s' % (R, N, np.round(z, 6)))
    for j in o[::12]:
        print('  r %6.2f |Ett res| %.1e scale %.1e | |aR1| %.1e |H0| %.1e |H1| %.1e |H2| %.1e' % (m.r[j], abs(res[j]), sc[j], abs(phys['aR1'][j]), abs(phys['H0'][j]), abs(phys['H1'][j]), abs(phys['H2'][j])))
