import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import sys, os, pickle, numpy as np
from shoot_modes import background
def get(A, R=60., eigen=0):
    fn = _os.environ.get('SBSE_WORK', '.') + '/bg_A%.4f_R%g.pkl' % (A, R)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    s, sol, d = background(_os.environ.get('SBSE_BINARY', _os.path.join(_HERE, '..', '..', 'BosonStarStability')), A, 0., False, 0.2, R, _os.environ.get('SBSE_WORK', '.') + '/shoot_%g' % A)
    d = {k: v for k, v in d.items() if k != 'sol'}
    pickle.dump((sol, d), open(fn, 'wb'))
    return sol, d
