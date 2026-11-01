import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..')); _sys.path.insert(0, _HERE)
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C1, C2 = '#2a78d6', '#eb6834'; ink, mute = '#2b2b2a', '#6f6e69'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': mute, 'axes.labelcolor': ink, 'xtick.color': mute, 'ytick.color': mute,
                     'axes.spines.top': False, 'axes.spines.right': False})
w2 = [0.96711, 0.95016, 0.93339, 0.91728, 0.90443, 0.89377, 0.88472, 0.87991, 0.87549, 0.87273, 0.87140, 0.87129]
l2 = [0.00174, 0.00309, 0.00448, 0.00566, 0.00629, 0.00638, 0.00584, 0.00511, 0.00387, 0.00239, 0.00069, 0.0]
w1 = [0.92983, 0.92645, 0.91728, 0.90443, 0.89377, 0.88472, 0.87549, 0.87011, 0.85882]
l1 = [0.0, 0.00080, 0.00175, 0.00282, 0.00366, 0.00434, 0.00499, 0.00535, 0.00596]
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.axvline(0.874, color=mute, lw=0.8, ls=':'); ax.text(0.8735, 0.0068, 'max mass', color=mute, fontsize=9)
ax.plot(w2, l2, '-s', color=C2, lw=2, ms=5, label='L = 2, real (quadrupolar, non-oscillatory)')
ax.plot(w1, l1, '-o', color=C1, lw=2, ms=5, label='L = 1, Re of oscillatory quartet')
ax.plot([0.87129], [0], 'D', color=ink, ms=6)
ax.annotate('L = 2 zero mode, ' r'$\omega = 0.8713$' '\n(2501.05342: 0.8713)', (0.87129, 0), (0.886, 0.0010), color=ink, fontsize=9,
            arrowprops=dict(arrowstyle='-', color=mute, lw=0.8))
ax.annotate('L = 1 onset\n' r'$\omega \approx 0.928$', (0.928, 0), (0.945, 0.0008), color=ink, fontsize=9,
            arrowprops=dict(arrowstyle='-', color=mute, lw=0.8))
ax.set_xlim(0.975, 0.85)
ax.set_xlabel(r'frequency $\omega/\mu$   (Newtonian limit on the left)'); ax.set_ylabel(r'growth rate Re $\lambda\ \ [\mu]$')
ax.set_title('Axisymmetric instabilities of spherical free Proca stars', color=ink, fontsize=11, loc='left')
ax.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig('proca_nonradial.png', dpi=160)
