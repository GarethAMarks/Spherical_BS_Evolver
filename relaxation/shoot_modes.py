"""
Radial-mode spectrum of individual Dirac stars that the shooting solver can find (free, quartic, and the degenerate
potential for sigma >~ 0.3): each star is shot with the main binary (cycle_only = 0, stop_time = 0), re-solved with
dirac_relax.py at tolerance 1e-8, and its spectrum computed with scan_modes.spectrum (two resolutions, tortoise map).

Usage: python3 shoot_modes.py <binary> <out.npy> <A1,A2,...> [lambda] [solitonic] [sigma] [R] [workdir]
Rows of <out.npy>: A omega M N C real_lam_max osc_Re_max sigma1 sigma2   (lambda = -i sigma; zeros when stable)
"""
import sys, os, re, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_relax import Potential, DiracRelax, seed_from_DSdata
from scan_modes import spectrum

PAR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'BSParams.par')


def shoot(binary, work, A, lam=0., solitonic=0, sigma=0.2, R=60., n=12001):
    os.makedirs(work, exist_ok=True)
    par = open(PAR).read()
    sub = {'dirac': 1, 'lambda': lam, 'solitonic': solitonic, 'sigma': sigma, 'A_central': A, 'cycle_only': 0, 'stop_time': 0,
           'R': R, 'n_gridpoints': n, 'frequency_guess': 20, 'perturb': 0, 'isotropic': 0, 'BS_resolution_factor': 1,
           'dirac_profile': 0}
    for k, v in sub.items():
        par = re.sub(r'^%s = [^\n#]*' % k, '%s = %s ' % (k, v), par, flags=re.M)
    open(os.path.join(work, 'BSParams.par'), 'w').write(par)
    out = subprocess.run([os.path.abspath(binary)], cwd=work, capture_output=True, text=True).stdout
    w = float(re.search(r'Final frequency after enforcing lapse condition is ([0-9.eE+-]+)', out).group(1))
    return os.path.join(work, 'DSdata.dat'), w


def background(binary, A, lam=0., solitonic=False, sigma=0.2, R=60., work='shoot_work'):
    """shooting solution re-solved by relaxation: returns (DiracRelax, solve_bvp solution, diagnostics)"""
    path, w = shoot(binary, work, A, lam, int(solitonic), sigma, R)
    s = DiracRelax(Potential(lam=lam, solitonic=solitonic, sigma=sigma), R=R, tol=1e-8)
    r = np.linspace(s.r0, R, 6001)
    y = seed_from_DSdata(path, r); y[4] = w
    sol, _ = s.solve_fixed_A(A, r, y)
    assert sol.status == 0, sol.message
    return s, sol, s.diagnostics(sol, A)


if __name__ == '__main__':
    binary, out, As = sys.argv[1], sys.argv[2], [float(a) for a in sys.argv[3].split(',')]
    lam = float(sys.argv[4]) if len(sys.argv) > 4 else 0.
    sol_on = bool(int(sys.argv[5])) if len(sys.argv) > 5 else False
    sg = float(sys.argv[6]) if len(sys.argv) > 6 else 0.2
    R = float(sys.argv[7]) if len(sys.argv) > 7 else 60.
    work = sys.argv[8] if len(sys.argv) > 8 else 'shoot_work'
    rows = []
    for A in As:
        try:
            s, sol, d = background(binary, A, lam, sol_on, sg, R, work)
            ev, _ = spectrum(s.pot, sol, 110, 150, alpha0=d['alpha0'])
        except Exception as e:
            print('A %g failed: %s' % (A, e), flush=True); continue
        un = ev[ev.real > 1e-5]
        real_un = un[np.abs(un.imag) < 1e-6]; osc = un[np.abs(un.imag) >= 1e-6]
        st = np.unique(np.round(np.abs(ev[np.abs(ev.real) <= 1e-5].imag), 7)); st = st[st > 1e-4]
        rows.append([A, d['omega'], d['M'], d['N'], d['C'], real_un.real.max() if len(real_un) else 0.,
                     osc.real.max() if len(osc) else 0., st[0] if len(st) else np.nan, st[1] if len(st) > 1 else np.nan])
        print('A %.4f omega %.5f M %.5f C %.4f | real unstable %s, oscillatory unstable %s | sigma %s' % (
            A, d['omega'], d['M'], d['C'], np.array2string(real_un.real, precision=5), np.array2string(osc[:2], precision=4),
            np.array2string(st[:3], precision=5)), flush=True)
        np.save(out, np.array(rows))
