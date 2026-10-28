"""
Plot families from families.py: compactness against frequency and mass, stable segments (turning-point analysis,
stability.py) solid, unstable dashed, the ends of the stable branches (maximum-mass turning points) marked, and the onset
of light rings. Usage: python3 plot_families.py <work_dir> <out.png> [sigma ...]
"""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stability import turning_points, mode_count, stable_segments

work, out = sys.argv[1], sys.argv[2]
sigmas = [float(x) for x in sys.argv[3:]] or sorted(float(os.path.basename(f)[7:-4]) for f in glob.glob(os.path.join(work, 'family_*.npz')))
sigmas = sorted(sigmas, reverse=True)
cmap = plt.get_cmap('viridis')
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
rows = []
for k, sg in enumerate(sigmas):
    d = np.load(os.path.join(work, 'family_%g.npz' % sg))
    M, N, w, C, A, nlr = d['M'], d['N'], d['omega'], d['C'], d['A'], d['n_light_rings']
    idx = turning_points(M, N)
    cnt, ch = mode_count(M, N, idx, w)
    col = cmap(k / max(len(sigmas) - 1, 1))
    lab = r'$\sigma = %g$' % sg
    for a_, (x, y) in enumerate([(w, C), (M * sg**2, C), (w, M * sg**2)]):
        axx = ax[a_]
        # unstable: dashed, stable: solid (split into runs of constant stability)
        st = cnt == 0
        for i0 in range(len(x) - 1):
            axx.plot(x[i0:i0 + 2], y[i0:i0 + 2], '-' if (st[i0] and st[i0 + 1]) else ':', color=col, lw=1.8 if st[i0] and st[i0 + 1] else 1.)
        if a_ == 0:
            axx.plot([], [], '-', color=col, lw=1.8, label=lab)
        if st[-1]:
            # the continuation stopped on the stable branch before reaching the maximum: open marker
            axx.plot(x[-1], y[-1], 'o', mfc='none', mec=col, ms=8, mew=1.5)
        segs_ = stable_segments(cnt)
        main_ = max(segs_, key=lambda ab: M[ab[1]])
        for a, b in segs_:
            if b < len(x) - 1 and b in idx:
                axx.plot(x[b], y[b], 'o', color=col, ms=8 if (a, b) == main_ else 4, mec='k')
        lr = np.where(nlr > 0)[0]
        if len(lr):
            axx.plot(x[lr[0]], y[lr[0]], '^', color=col, ms=6, mec='k')
        if 'dec' in d.files:
            for key, thr, mk in (('dec', 0.01, 'x'), ('wec', -1e-3, '+')):
                v = d[key]
                bad = np.where(v > thr)[0] if key == 'dec' else np.where(v < thr)[0]
                if len(bad):
                    axx.plot(x[bad[0]], y[bad[0]], mk, color='k', ms=7, mew=1.5)
    segs = stable_segments(cnt)
    main = max(segs, key=lambda ab: M[ab[1]])
    a, b = main
    lr_st = [i for i in range(a, b + 1) if nlr[i] > 0]
    # energy conditions counted as violated beyond 1% (dominant) / 0.1% (weak) of the maximum energy density
    dec_st = [i for i in range(a, b + 1) if 'dec' in d.files and d['dec'][i] > 0.01]
    wec_st = [i for i in range(a, b + 1) if 'wec' in d.files and d['wec'][i] < -1e-3]
    rows.append((sg, M[b], M[b] * sg**2, w[b], C[b], d['Cmax'][b], A[b] / (sg / np.sqrt(2)), M[lr_st[0]] if lr_st else np.nan,
                 C[lr_st[0]] if lr_st else np.nan, C[dec_st[0]] if dec_st else np.nan, C[wec_st[0]] if wec_st else np.nan, C.max()))
ax[0].set_xscale('log'); ax[0].set_xlabel(r'$\omega/\mu$'); ax[0].set_ylabel(r'$C = M / r_{99}$')
ax[1].set_xscale('log'); ax[1].set_xlim(1e-5, 0.5); ax[1].set_xlabel(r'$M \sigma^2$  ($G = \mu = 1$)'); ax[1].set_ylabel(r'$C = M / r_{99}$')
for axx in ax[:2]:
    axx.set_ylim(0, 0.47)
    axx.axhline(4 / 9, color='grey', lw=0.6); axx.text(axx.get_xlim()[0] * 1.1, 4 / 9 + 0.005, 'Buchdahl 4/9', color='grey', fontsize=8)
ax[0].legend(fontsize=8, loc='upper right')
ax[2].set_xscale('log'); ax[2].set_yscale('log'); ax[2].set_xlabel(r'$\omega/\mu$'); ax[2].set_ylabel(r'$M \sigma^2$')
ax[2].set_title('mass: stability changes at the extrema of M (= extrema of N)', fontsize=9)
ax[0].set_title('solid: stable (turning points), dotted: unstable; filled o: end of stable branch (open: not reached),\n^: light rings appear, x: dominant, +: weak energy condition first violated', fontsize=9)
fig.suptitle(r'Dirac stars, degenerate-vacuum potential $W(S) = \mu S (1 - S/S_0)^2$, $S_0 = \sigma^2/4\pi$')
fig.tight_layout(); fig.savefig(out, dpi=130)
print('end of the stable branch (maximum mass), and the onsets along it of light rings and of dominant / weak energy condition violation')
print('sigma   M_max   M_max s^2   omega   C(M_max)  Cmax(M_max)  A/A_vac  M at LR onset  C at LR onset  C at DEC  C at WEC  max C on family')
for r in rows:
    print('%-6g  %8.3f  %.4f  %.5f  %.4f  %.4f  %.3f  %8.3f  %.4f  %.4f  %.4f  %.4f' % r)
