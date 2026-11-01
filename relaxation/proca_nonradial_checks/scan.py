import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, pickle, time, os
from proca_nonradial_modes import classified
from fam import get
As = [float(a) for a in sys.argv[1].split(',')]; L = int(sys.argv[2]); N1, N2 = (int(x) for x in sys.argv[3].split(','))
out = sys.argv[4]; Rfix = float(sys.argv[5]) if len(sys.argv) > 5 else None
res = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
for A in As:
    sol, d = get(A)
    R = Rfix or sol.x[-1]
    t0 = time.time()
    cl = classified(sol, L, N1, N2, R=R, lam_max=0.6)
    res[(A, L)] = dict(cl=cl, d=d, R=R)
    good = [(z, rr) for z, c, rr in cl if c and rr < 1e-2]
    un = [(z, rr) for z, rr in good if z.real > 1e-3 and z.imag >= 0]
    low = sorted([(z, rr) for z, rr in good if z.imag > 1e-6 and abs(z.real) < 1e-4], key=lambda t: t[0].imag)[:3]
    cand = [(z, c, rr) for z, c, rr in cl if z.real > 1e-3 and z.imag >= 0 and not (c and rr < 1e-2)]
    print('A %.4f omega %.5f M %.5f L %d R %g | unstable %s | lowest sigma %s | rejected %s (%.0fs)' % (
        A, d['omega'], d['M'], L, R, ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in un),
        ' '.join('%.5f' % z.imag for z, rr in low), ' '.join('%.4f%+.4fi[%s,%.0e]' % (z.real, z.imag, 'c' if c else '-', rr) for z, c, rr in cand[:3]),
        time.time() - t0), flush=True)
    pickle.dump(res, open(out, 'wb'))
