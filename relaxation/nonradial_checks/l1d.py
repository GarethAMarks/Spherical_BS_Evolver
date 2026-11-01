import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, numpy as np, time
from dirac_nonradial_modes import NonRadialModes
from bg import get
from bg_exc import get_exc
exc = sys.argv[3] == 'exc' if len(sys.argv) > 3 else False
for A in [float(a) for a in sys.argv[1].split(',')]:
    sol, d = get_exc(A) if exc else get(A)
    R = float(sys.argv[2])
    out = []
    t0 = time.time()
    for N in (110, 150):
        m = NonRadialModes(sol, 1, N=N, R=R)
        l, V = m.eigen(vectors=True)
        s = (np.abs(l) < 0.6) & (np.abs(l) > 1e-8) & (l.imag >= -1e-12); l = l[s]; V = V[:, s]
        items = []
        for k in range(len(l)):
            rs = m.residuals(l[k], V[:, k]); items.append((l[k], max(rs['Ett'], rs['Err'], rs['Ea-Eb'])))
        out.append(items)
    # physical: converged between N and small residual
    l1 = np.array([z for z, _ in out[0]])
    phys = [(z, rr) for z, rr in out[1] if rr < 1e-2 and len(l1) and np.min(np.abs(l1 - z)) < 2e-3 * abs(z) + 3e-4]
    un = sorted([p for p in phys if p[0].real > 5e-4], key=lambda t: -t[0].real)
    low = sorted([p for p in phys if abs(p[0].real) < 1e-4 and p[0].imag > 1e-6], key=lambda t: t[0].imag)[:4]
    cand = sorted([(z, rr) for z, rr in out[1] if z.real > 1e-3], key=lambda t: -t[0].real)[:3]
    print('A %.3f omega %.5f %s| unstable %s | lowest sigma %s | raw Re>1e-3 %s (%.0fs)' % (A, d['omega'], 'exc ' if exc else '',
          ' '.join('%.5f%+.5fi(%.0e)' % (z.real, z.imag, rr) for z, rr in un),
          ' '.join('%.5f(%.0e)' % (z.imag, rr) for z, rr in low), ' '.join('%.4f%+.4fi(%.0e)' % (z.real, z.imag, rr) for z, rr in cand), time.time() - t0), flush=True)
