import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, pickle, os, numpy as np
from proca_relax import background
def get(A, R=None):
    R = R or (120. if A < 0.006 else (90. if A < 0.015 else 60.))
    fn = _os.environ.get('SBSE_WORK', '.') + '/bg_A%.5f_R%g.pkl' % (A, R)
    if os.path.exists(fn): return pickle.load(open(fn, 'rb'))
    s, sol, d = background(_os.environ.get('SBSE_BINARY', _os.path.join(_HERE, '..', '..', 'BosonStarStability')), _os.path.join(_HERE, '..', '..', 'BSParams.par'), A, R=R, work='sh_%g' % A)
    pickle.dump((sol, d), open(fn, 'wb')); return sol, d
if __name__ == '__main__':
    for A in [float(a) for a in sys.argv[1].split(',')]:
        try:
            sol, d = get(A)
            print('A %.4f omega %.5f M %.5f C %.4f alpha0 %.3f r99 %.1f' % (A, d['omega'], d['M'], d['C'], d['alpha0'], d['r99']), flush=True)
        except Exception as e:
            print(A, 'failed', e, flush=True)
