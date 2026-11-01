import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, pickle, time
from dirac_nonradial_modes import classified
from bg_exc import get_exc
As = [float(a) for a in sys.argv[1].split(',')]; L = int(sys.argv[2]); R = float(sys.argv[3]); N1, N2 = int(sys.argv[4]), int(sys.argv[5])
out = {}
for A in As:
    t0 = time.time()
    sol, d = get_exc(A, R=max(80., R))
    cl = classified(sol, L, N1, N2, R=R, lam_max=0.08)
    g = [(z, rr) for z, c, rr in cl if c and rr < 1e-2 and (z.imag >= 0)]
    g.sort(key=lambda t: abs(t[0]))
    out[A] = (d, cl)
    print('A %.4f omega %.5f M %.5f | near-zero physical: %s (%.0fs)' % (A, d['omega'], d['M'],
          ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in g[:5]), time.time() - t0), flush=True)
    pickle.dump(out, open('track0_L%d_R%g.pkl' % (L, R), 'wb'))
