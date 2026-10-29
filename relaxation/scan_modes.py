"""
Radial-mode spectrum along a family computed by families.py (run with profile_stride 1 or small), using dirac_modes.py.
For each stored profile (every <stride>-th, with A/A_vac in [a_lo, a_hi]) the eigenvalues are computed at two
resolutions and only those that agree are kept. Recorded per model: the most unstable eigenvalue (largest Re lambda,
with its Im), whether it is real or oscillatory, and the lowest stable frequencies sigma (lambda = -i sigma).
The collocation map mixes in the tortoise coordinate (resolving redshifted cores and walls); the resolution is raised for
small central lapse. Weakly unstable high-frequency modes (Re lambda ~ 0.01, |Im lambda| above the continuum threshold) can
be artefacts of the reflecting outer boundary: their growth rate falls as the domain grows (check with dirac_modes.ModeSolver
and a larger R).

Usage: python3 scan_modes.py <dir> <sigma> <out.npz> [stride] [a_lo] [a_hi] [N1] [N2]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from export_profile import load, resolve
from dirac_modes import ModeSolver


def wall(sol):
    x = np.linspace(sol.x[0], sol.x[-1], 20001)
    g = sol.sol(x)[1]
    dg = np.abs(np.gradient(g, x))
    j = dg.argmax()
    return x[j], np.abs(g).max() / dg[j]


def spectrum(pot, sol, N1, N2, lam_max=2.5, tol=2e-3, alpha0=1.):
    """eigenvalues that agree between N1 and N2 collocation points (tortoise-coordinate map, which resolves redshifted
    cores and walls); N is raised for deep potential wells (small central lapse)"""
    fac = 1. if alpha0 > 0.3 else (1.4 if alpha0 > 0.05 else 1.8)
    kw = dict(tortoise=0.8)
    ls = []
    for N in (int(N1 * fac), int(N2 * fac)):
        l = ModeSolver(pot, sol, N=N, **kw).eigenvalues()
        ls.append(l[(np.abs(l) < lam_max) & (np.abs(l) > 1e-4)])
    good = np.array([z for z in ls[0] if len(ls[1]) and np.min(np.abs(ls[1] - z)) < tol * max(1., abs(z))])
    # keep only frequencies resolved by the coarser grid (about 6 points per wavelength of the outer, flat-space region)
    k_res = 2. * int(N1 * fac) / (3. * sol.x[-1])
    good = good[np.abs(good.imag) < k_res] if len(good) else good
    return good, kw


if __name__ == '__main__':
    work, sg, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    a_lo = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    a_hi = float(sys.argv[6]) if len(sys.argv) > 6 else 10.
    N1 = int(sys.argv[7]) if len(sys.argv) > 7 else 120
    N2 = int(sys.argv[8]) if len(sys.argv) > 8 else 170
    prof, fam, count = load(work, sg)
    Av = sg / np.sqrt(2)
    stored_stride = max(1, round(len(fam['A']) / len(prof)))
    rows = []
    for k in range(0, len(prof), stride):
        A, x, y = prof[k]
        if not (a_lo <= A / Av <= a_hi):
            continue
        i = min(k * stored_stride, len(fam['A']) - 1)
        t0 = time.time()
        try:
            s, sol, d = resolve(sg, A, x, y)
            lam, kw = spectrum(s.pot, sol, N1, N2, alpha0=d['alpha0'])
        except Exception as e:
            print('profile %d failed: %s' % (k, e), flush=True); continue
        un = lam[lam.real > 1e-5]
        top = un[np.argmax(un.real)] if len(un) else 0j
        st = np.unique(np.round(np.abs(lam[np.abs(lam.real) <= 1e-5].imag), 7))
        st = st[st > 1e-4]
        real_un = un[np.abs(un.imag) < 1e-6]
        osc_un = un[np.abs(un.imag) >= 1e-6]
        top_osc = osc_un[np.argmax(osc_un.real)] if len(osc_un) else 0j
        rows.append([A, A / Av, d['omega'], d['M'], d['N'], d['C'], d['Cmax'], count[i], top.real, abs(top.imag),
                     len(un[un.imag >= 0]), *(list(st[:3]) + [np.nan] * (3 - len(st[:3]))),
                     real_un.real.max() if len(real_un) else 0., top_osc.real, abs(top_osc.imag), d['alpha0']])
        print('k %3d A/Avac %.4f M %.4f C %.4f TP-modes %d | max Re lambda %.5f (Im %.4f), %d unstable | sigma %s  %s %.0fs' % (
            k, A / Av, d['M'], d['C'], count[i], top.real, abs(top.imag), len(un[un.imag >= 0]),
            np.array2string(st[:3], precision=4), '', time.time() - t0), flush=True)
        np.savez(out, sigma=sg, data=np.array(rows),
                 columns='A A/Avac omega M N C Cmax TP_modes Re_lam_max Im_lam_max n_unstable sigma1 sigma2 sigma3 '
                         'real_lam_max osc_Re_max osc_Im alpha0')
