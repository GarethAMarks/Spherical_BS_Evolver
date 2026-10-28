"""
Families of Dirac stars with the degenerate-vacuum potential mu S (1 - S/S0)^2 (S0 = sigma^2 / 4 pi), by pseudo-arclength
continuation with the relaxation solver (dirac_relax.py), seeded by one shooting solution of the C++ DiracStar solver at a
low central amplitude (where shooting is reliable).

Usage: python3 families.py <sigma> <work_dir> <path to the BosonStarStability binary> [A_seed/A_vac] [A_max/A_vac]
Writes <work_dir>/family_<sigma>.npz with the diagnostics along the family (in continuation order).
"""
import sys, os, re, subprocess, time, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_relax import Potential, DiracRelax, seed_from_DSdata, continue_family, equidistribute

PAR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'BSParams.par')


def shoot_seed(sigma, A, work, binary, R=100., n=10001):
    """one DiracStar shooting solve (cycle_only = 0, stop_time = 0) at central amplitude A; returns (DSdata path, omega)"""
    os.makedirs(work, exist_ok=True)
    par = open(PAR).read()
    sub = {'dirac': 1, 'lambda': 0, 'solitonic': 1, 'sigma': sigma, 'A_central': A, 'cycle_only': 0, 'stop_time': 0,
           'R': R, 'n_gridpoints': n, 'frequency_guess': 20, 'perturb': 0, 'isotropic': 0}
    for k, v in sub.items():
        par = re.sub(r'^%s = [^\n#]*' % k, '%s = %s ' % (k, v), par, flags=re.M)
    open(os.path.join(work, 'BSParams.par'), 'w').write(par)
    out = subprocess.run([binary], cwd=work, capture_output=True, text=True).stdout
    w = float(re.search(r'Final frequency after enforcing lapse condition is ([0-9.eE+-]+)', out).group(1))
    return os.path.join(work, 'DSdata.dat'), w


def family(sigma, work, binary, seed_ratio=0.5, A_max_ratio=4., A_min_ratio=1e-4, ds_max=0.05, n_steps=600, verbose=True):
    A_vac = sigma / np.sqrt(2.)
    A_seed = seed_ratio * A_vac
    path, w = shoot_seed(sigma, A_seed, os.path.join(work, 'seed_%g' % sigma), binary)
    pot = Potential(sigma=sigma)
    s = DiracRelax(pot, R=100.)
    r = np.linspace(s.r0, s.R, 4001)
    y = seed_from_DSdata(path, r); y[4] = w
    sol, A = s.solve_fixed_A(A_seed, r, y)
    assert sol.status == 0, sol.message
    # outer radius: beyond the star by 25 decay lengths of the spinor tail, exp(-sqrt(mu^2 - omega^2) r)
    R_rule = lambda d: float(1.2 * d['r99'] + 25. / np.sqrt(max(1. - d['omega']**2, 1e-6)))
    branches = []
    for direction, stop, dsm in ((+1, lambda d: d['A'] > A_max_ratio * A_vac, ds_max),
                                 (-1, lambda d: d['A'] < A_min_ratio * A_vac or d['omega'] > 0.99995, 0.2 * ds_max)):
        t = time.time()
        s.R = 100.
        out = continue_family(s, sol, A, n_steps=n_steps, ds=0.01, direction=direction, scales=(A_vac, 1.), R_rule=R_rule,
                              stop=stop, verbose=verbose, ds_max=dsm)
        print('sigma %g direction %+d: %d points, A %.4g .. %.4g, %.0fs' % (sigma, direction, len(out), out[0]['A'], out[-1]['A'],
                                                                         time.time() - t), flush=True)
        branches.append(out)
    fam = branches[1][::-1] + branches[0][1:]          # ordered from low to high amplitude along the arc
    keys = ['A', 'omega', 'M', 'N', 'Eb', 'r99', 'C', 'Cmax', 'r_Cmax', 'alpha0', 'nodes', 'wec', 'dec']
    data = {k: np.array([d[k] for d in fam]) for k in keys}
    data['n_light_rings'] = np.array([len(d['light_rings']) for d in fam])
    data['light_rings'] = np.array([np.pad(d['light_rings'][:2], (0, 2 - min(2, len(d['light_rings']))), constant_values=np.nan)
                                    for d in fam])
    data['sigma'] = sigma
    np.savez(os.path.join(work, 'family_%g.npz' % sigma), **data)
    # keep a few profiles (every 10th point) for inspection
    prof = [(d['A'], d['sol'].x, d['sol'].y) for d in fam[::10]]
    pickle.dump(prof, open(os.path.join(work, 'profiles_%g.pkl' % sigma), 'wb'))
    return data


if __name__ == '__main__':
    sigma, work, binary = float(sys.argv[1]), sys.argv[2], sys.argv[3]
    seed_ratio = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    A_max_ratio = float(sys.argv[5]) if len(sys.argv) > 5 else 2.5
    family(sigma, work, binary, seed_ratio=seed_ratio, A_max_ratio=A_max_ratio)
