import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import numpy as np, pickle, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
ink, mute = '#2b2b2a', '#6f6e69'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': mute, 'axes.labelcolor': ink, 'xtick.color': mute, 'ytick.color': mute,
                     'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
# (a) ground state: signed lambda^2 (positive = unstable, -sigma^2 = stable)
rad = np.load('radial_free.npy')
A = rad[:, 0]; s_rad = np.where(rad[:, 5] > 0, rad[:, 5]**2, -rad[:, 7]**2)
res = pickle.load(open('scan_ground.pkl', 'rb'))
AL, sL = [], []
for (a, L), v in sorted(res.items()):
    if L != 2: continue
    edge = 1 - v['omega']
    b = sorted(abs(z.imag) for z, c, rr in v['cl'] if c and rr < 1e-2 and abs(z.real) < 1e-4 and 1e-4 < abs(z.imag) < edge)
    if b: AL.append(a); sL.append(-b[0]**2)
a0 = ax[0]
a0.axhline(0, color=mute, lw=0.8); a0.axvline(0.23, color=mute, lw=0.8, ls=':')
a0.plot(A, s_rad * 1e3, '-o', color=C1, lw=2, ms=5, label='radial (L = 0) fundamental')
a0.plot(AL, np.array(sL) * 1e3, '-s', color=C2, lw=2, ms=5, label='L = 2, lowest bound mode')
a0.text(0.232, a0.get_ylim()[1] * 0.85 if False else 3.2, 'max mass', color=mute, fontsize=9)
a0.set_xlabel(r'central amplitude $A$'); a0.set_ylabel(r'$\lambda^2\ \ [10^{-3}\,\mu^2]$   (> 0 unstable)')
a0.set_title('Ground state (free Dirac stars)', color=ink, fontsize=11, loc='left')
a0.legend(frameon=False, fontsize=9, loc='lower left')
# (b) excited state growth rates
exc = pickle.load(open('radial_exc.pkl', 'rb'))
Ae = [r[0] for r in exc]; osc = [max([z.real for z in r[4] if abs(z.imag) > 1e-6] or [0]) for r in exc]
L2A = [0.025, 0.03, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.21, 0.22, 0.23, 0.24, 0.245, 0.25, 0.2521]
L2l = [0.00244, 0.00291, 0.00474, 0.00719, 0.00918, 0.01054, 0.01107, 0.01044, 0.00987, 0.00904, 0.00785, 0.00607, 0.00476, 0.00267, 0.0]
a1 = ax[1]
a1.plot(Ae, osc, '-o', color=C1, lw=2, ms=5, label='radial (L = 0), oscillatory')
a1.plot(L2A, L2l, '-s', color=C2, lw=2, ms=5, label='L = 2, real (non-oscillatory)')
a1.plot([0.2521], [0], 'D', color=ink, ms=6)
a1.annotate(r'L = 2 zero mode' '\n' r'$A\approx0.252,\ \omega\approx0.850$', (0.2521, 0), (0.27, 0.012), color=ink, fontsize=9,
            arrowprops=dict(arrowstyle='-', color=mute, lw=0.8))
a1.axvline(0.2, color=mute, lw=0.8, ls=':'); a1.text(0.203, 0.041, 'max mass', color=mute, fontsize=9)
a1.set_xlabel(r'central amplitude $A$'); a1.set_ylabel(r'growth rate Re $\lambda\ \ [\mu]$')
a1.set_title('First excited state (one node)', color=ink, fontsize=11, loc='left')
a1.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout()
fig.savefig('dirac_star_modes.png', dpi=160)
