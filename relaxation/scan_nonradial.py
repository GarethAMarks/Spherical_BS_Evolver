"""
Non-radial spectra (polar, M = 0, L >= 1) of free spherical Dirac stars found by shooting (shoot_modes.background), using
dirac_nonradial_modes.NonRadialModes with all Einstein components imposed (overdetermined pencil, random projections).
An eigenvalue is kept if it agrees between two resolutions N1, N2.

Usage: python3 scan_nonradial.py <binary> <out.npz> <A1,A2,...> <L1,L2,...> [R] [N1] [N2] [eigen]
Records, per star and L: the kept eigenvalues with |lambda| < lam_max.
"""
import sys, os, time, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_nonradial_modes import NonRadialModes, classified
from shoot_modes import background


def robust_full(sol, L, N1=80, N2=110, R=None, lam_max=2.5, tol=3e-3, atol=2e-4):
    ls = []
    for N in (N1, N2):
        m = NonRadialModes(sol, L, N=N, R=R)
        l = m.eigen_projected(tol=2e-3, full=True)
        ls.append(l[np.abs(l) < lam_max])
    keep = []
    for z in ls[0]:
        d = np.abs(ls[1] - z)
        if len(d) and d.min() < tol * abs(z) + atol:
            keep.append(0.5 * (z + ls[1][d.argmin()]))
    keep = np.array(keep)
    return keep[np.argsort(-keep.real)] if len(keep) else keep, ls


if __name__ == '__main__':
    binary, out = sys.argv[1], sys.argv[2]
    As = [float(a) for a in sys.argv[3].split(',')]
    Ls = [int(a) for a in sys.argv[4].split(',')]
    R = float(sys.argv[5]) if len(sys.argv) > 5 else 60.
    N1 = int(sys.argv[6]) if len(sys.argv) > 6 else 80
    N2 = int(sys.argv[7]) if len(sys.argv) > 7 else 110
    eigen = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    res = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
    for A in As:
        try:
            s, sol, d = background(binary, A, 0., False, 0.2, R, 'shoot_nr_%g' % A)
        except Exception as e:
            print('A %g: background failed: %s' % (A, e), flush=True); continue
        for L in Ls:
            t0 = time.time()
            cl = classified(sol, L, N1, N2, R=R)
            res[(A, L)] = dict(cl=cl, omega=d['omega'], M=d['M'], C=d['C'], R=R, N=(N1, N2))
            good = [z for z, c, rr in cl if c and rr < 1e-2]
            g = np.array(good)
            un = g[g.real > 1e-4] if len(g) else g
            cand = [(z, c, rr) for z, c, rr in cl if z.real > 1e-3]
            low = np.unique(np.round(np.abs(g.imag), 5)) if len(g) else []
            print('A %.4f omega %.5f M %.5f C %.4f L %d | physical unstable %s | lowest |Im| %s | %d rejected Re>1e-3 '
                  '(min resid %.0e) (%.0fs)' % (A, d['omega'], d['M'], d['C'], L, np.array2string(un[:4], precision=5),
                  np.array2string(np.array(low[:5]), precision=5), len(cand), min([rr for _, _, rr in cand]) if cand else np.nan,
                  time.time() - t0), flush=True)
            pickle.dump(res, open(out, 'wb'))
