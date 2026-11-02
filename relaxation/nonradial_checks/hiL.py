import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, time
import dirac_nonradial_modes as pm
import dirac_axial_modes as am
from bg import get
from bg_exc import get_exc
A = float(sys.argv[1]); exc = sys.argv[2] == 'exc'; sector = sys.argv[3]; Ls = [int(x) for x in sys.argv[4].split(',')]
R = float(sys.argv[5]); rs = float(sys.argv[6]) if len(sys.argv) > 6 else 2.
N1, N2 = (110, 150) if len(sys.argv) <= 7 else tuple(int(x) for x in sys.argv[7].split(','))
sol, d = get_exc(A) if exc else get(A)
mod = pm if sector == 'polar' else am
for L in Ls:
    t0 = time.time()
    cl = mod.classified(sol, L, N1, N2, R=R, lam_max=1.0, r_s=rs)
    good = [(z, rr) for z, c, rr in cl if c and rr < 1e-2]
    un = [(z, rr) for z, rr in good if z.real > 5e-4 and z.imag >= 0]
    low = sorted([(z, rr) for z, rr in good if abs(z.real) < 1e-4 and z.imag > 1e-6], key=lambda t: t[0].imag)[:3]
    cand = sorted([(z, c, rr) for z, c, rr in cl if z.real > 5e-4 and z.imag >= 0], key=lambda t: -t[0].real)[:3]
    print('A %.3f %s %s L %2d R %g rs %g | unstable %s | lowest sigma %s | n_conv %d | raw Re>5e-4 %s (%.0fs)' % (
        A, 'exc' if exc else 'gs ', sector, L, R, rs, ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in un),
        ' '.join('%.5f' % z.imag for z, rr in low), sum(1 for z, c, rr in cl if c), ' '.join('%.4f%+.4fi[%s,%.0e]' % (z.real, z.imag, 'c' if c else '-', rr) for z, c, rr in cand),
        time.time() - t0), flush=True)
