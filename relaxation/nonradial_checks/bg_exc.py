import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, os, re, subprocess, pickle, numpy as np
from dirac_relax import Potential, DiracRelax, seed_from_DSdata
PAR = _os.path.join(_HERE, '..', '..', 'BSParams.par')
def shoot(A, eigen, R=80., n=16001, work=None, fguess=20):
    work = work or _os.environ.get('SBSE_WORK', '.') + '/exc_%d_%g' % (eigen, A)
    os.makedirs(work, exist_ok=True)
    par = open(PAR).read()
    sub = {'dirac': 1, 'lambda': 0, 'solitonic': 0, 'sigma': 0.2, 'A_central': A, 'cycle_only': 0, 'stop_time': 0, 'eigen': eigen,
           'R': R, 'n_gridpoints': n, 'frequency_guess': fguess, 'perturb': 0, 'isotropic': 0, 'BS_resolution_factor': 1, 'dirac_profile': 0}
    for k, v in sub.items():
        par = re.sub(r'^%s = [^\n#]*' % k, '%s = %s ' % (k, v), par, flags=re.M)
    open(os.path.join(work, 'BSParams.par'), 'w').write(par)
    out = subprocess.run([_os.environ.get('SBSE_BINARY', _os.path.join(_HERE, '..', '..', 'BosonStarStability'))], cwd=work, capture_output=True, text=True).stdout
    w = float(re.search(r'Final frequency after enforcing lapse condition is ([0-9.eE+-]+)', out).group(1))
    return os.path.join(work, 'DSdata.dat'), w
def get_exc(A, eigen=1, R=80.):
    fn = _os.environ.get('SBSE_WORK', '.') + '/bgexc%d_A%.4f_R%g.pkl' % (eigen, A, R)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    Rs = min(R, 80.)
    path, w = shoot(A, eigen, Rs)
    s = DiracRelax(Potential(lam=0., solitonic=False), R=Rs, tol=1e-8)
    r = np.linspace(s.r0, Rs, 8001)
    y = seed_from_DSdata(path, r); y[4] = w
    sol, _ = s.solve_fixed_A(A, r, y)
    assert sol.status == 0, sol.message
    if R > Rs:
        from dirac_relax import extend
        x, y = extend(sol, R, s)
        s.R = R
        sol, _ = s.solve_fixed_A(A, x, y)
        assert sol.status == 0, sol.message
    d = s.diagnostics(sol, A)
    pickle.dump((sol, d), open(fn, 'wb'))
    return sol, d
if __name__ == '__main__':
    for A in [float(a) for a in sys.argv[1].split(',')]:
        try:
            sol, d = get_exc(A)
            print('A %.3f omega %.5f M %.5f N %.5f C %.4f g_nodes %d alpha0 %.3f' % (A, d['omega'], d['M'], d['N'], d['C'], d['g_nodes'], d['alpha0']))
        except Exception as e:
            print('A', A, 'failed', e)
