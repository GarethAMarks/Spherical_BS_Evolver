"""
Export a star computed by families.py for evolution with the C++ code (DiracStar, dirac_profile = 1).

families.py stores every 10th solution of a family in profiles_<sigma>.pkl. This script lists them with their properties,
or re-solves one (at its central amplitude, from the stored solution) and writes it on a uniform radial grid, using the
collocation interpolant, in the format read by DiracStar::read_profile_solution:
    '#' header lines including  omega = , sigma = , lambda = , mu = , solitonic =  (and M, N, C, ... for reference),
    then the columns  r f g m phi  from r = 0 (f, g the physical spinor amplitudes, m the aspect mass, phi = ln alpha).

Usage:
    python3 export_profile.py <dir> <sigma> list
    python3 export_profile.py <dir> <sigma> <profile index> <out file> [dr = 0.002]
"""
import sys, os, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dirac_relax import Potential, DiracRelax, equidistribute
from stability import turning_points, mode_count


def load(work, sigma):
    prof = pickle.load(open(os.path.join(work, 'profiles_%g.pkl' % sigma), 'rb'))
    fam = np.load(os.path.join(work, 'family_%g.npz' % sigma))
    M, N, w = fam['M'], fam['N'], fam['omega']
    count, _ = mode_count(M, N, turning_points(M, N), w)
    return prof, fam, count


def resolve(sigma, A, x, y, n_mesh=6000):
    s = DiracRelax(Potential(sigma=sigma), R=x[-1], tol=1e-8)
    class _S: pass
    tmp = _S(); tmp.x = x; tmp.sol = lambda xx: np.vstack([np.interp(xx, x, yy) for yy in y])
    xx, yy = equidistribute(tmp, n_mesh)
    sol, _ = s.solve_fixed_A(A, xx, yy)
    assert sol.status == 0, sol.message
    return s, sol, s.diagnostics(sol, A)


if __name__ == '__main__':
    work, sigma = sys.argv[1], float(sys.argv[2])
    prof, fam, count = load(work, sigma)
    if sys.argv[3] == 'list':
        print(' idx  A/A_vac   omega      M        C      Cmax    LR  modes   wec      dec')
        for k in range(len(prof)):
            i = 10 * k
            print('%4d  %.4f  %.6f  %8.4f  %.4f  %.4f  %d   %d   %+.2e  %+.2e' % (
                k, fam['A'][i] / (sigma / np.sqrt(2)), fam['omega'][i], fam['M'][i], fam['C'][i], fam['Cmax'][i],
                fam['n_light_rings'][i], count[i], fam['wec'][i], fam['dec'][i]))
        sys.exit(0)
    k, out = int(sys.argv[3]), sys.argv[4]
    dr = float(sys.argv[5]) if len(sys.argv) > 5 else 0.002
    A, x, y = prof[k]
    s, sol, d = resolve(sigma, A, x, y)
    r = np.arange(0., sol.x[-1], dr)
    yy = sol.sol(np.clip(r, sol.x[0], None))
    f, g, m, phi = yy[0], yy[1], yy[2], yy[3]
    f[0], m[0] = 0., 0.
    hdr = ('Dirac star (degenerate-vacuum potential) from relaxation/families.py, profile %d of the family\n'
           'omega = %.15g\nsigma = %.15g\nlambda = 0\nmu = 1\nsolitonic = 1\n'
           'A_central = %.12g  M = %.12g  N = %.12g  r_99 = %.10g  C = %.6f  Cmax = %.6f  C4 = %.6f  light rings at %s\n'
           'unstable modes (turning points) = %d  wec = %.3e  dec = %.3e\n'
           'r f g m phi' % (k, d['omega'], sigma, A, d['M'], d['N'], d['r99'], d['C'], d['Cmax'], d['C4'],
                            np.array2string(np.round(d['light_rings'], 4)), count[10 * k], d['wec'], d['dec']))
    np.savetxt(out, np.column_stack([r, f, g, m, phi]), header=hdr, fmt='%.15e')
    print('wrote %s: A = %.8g omega = %.10g M = %.8g C = %.4f Cmax = %.4f LR %s modes %d, %d points to r = %.1f' % (
        out, A, d['omega'], d['M'], d['C'], d['Cmax'], np.round(d['light_rings'], 3), count[10 * k], len(r), r[-1]))
