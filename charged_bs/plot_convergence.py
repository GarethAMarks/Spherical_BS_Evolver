#!/usr/bin/env python3
"""
Constraint convergence plots for SBSE runs at three resolutions dr, dr/2, dr/4 (same star, same R).

  python3 plot_convergence.py --runs low/ mid/ high/ --order 4 --label "stable star" --out conv.png

Reads constraint_norms.dat in each run directory and plots the L2 norms of the Hamiltonian, momentum and (if present)
Gauss constraints, with the medium and high resolutions multiplied by 2^order and 4^order: for convergence at that
order the three curves coincide. Also plots the convergence factor log2(||C_low|| / ||C_mid||) and
log2(||C_mid|| / ||C_high||) against time, and the central amplitude, lapse and apparent-horizon radius.
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(d):
    f = os.path.join(d, 'constraint_norms.dat')
    with open(f) as fh:
        header = fh.readline()
    dat = np.genfromtxt(f, comments='#', invalid_raise=False)
    cols = dict(t=0, ham=1, mom=2, A=4, ah=5, lapse=10)
    out = {k: dat[:, i] for k, i in cols.items()}
    if 'Gauss_L2' in header:
        out['gauss'] = dat[:, 19]
        out['Q'] = dat[:, 20]
    return out


def common(t_ref, t, y):
    return np.interp(t_ref, t, y, left=np.nan, right=np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs=3, required=True, help='low, medium, high resolution run directories')
    ap.add_argument('--order', type=float, default=4.0)
    ap.add_argument('--label', default='')
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--out', default='convergence.png')
    a = ap.parse_args()

    R = [load(d) for d in a.runs]
    tmax = min(r['t'][-1] for r in R) if a.tmax is None else a.tmax
    t = R[0]['t'][R[0]['t'] <= tmax]
    f = 2.0**a.order
    keys = [('ham', 'Hamiltonian'), ('mom', 'momentum')] + ([('gauss', 'Gauss')] if 'gauss' in R[0] else [])

    fig, ax = plt.subplots(3, len(keys), figsize=(4.6 * len(keys), 10), sharex=True)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    names = ['dr', 'dr/2 (x%g)' % f, 'dr/4 (x%g)' % f**2]
    for c, (k, title) in enumerate(keys):
        ys = [common(t, r['t'], r[k]) for r in R]
        for i in range(3):
            ax[0, c].semilogy(t, ys[i] * f**i, color=colors[i], label=names[i], lw=1.2)
        ax[0, c].set_title('%s constraint (L2)' % title)
        ax[0, c].legend(fontsize=8)
        with np.errstate(divide='ignore', invalid='ignore'):
            ax[1, c].plot(t, np.log2(ys[0] / ys[1]), color=colors[1], label='log2(low/mid)', lw=1)
            ax[1, c].plot(t, np.log2(ys[1] / ys[2]), color=colors[2], label='log2(mid/high)', lw=1)
        ax[1, c].axhline(a.order, color='k', ls='--', lw=0.8)
        ax[1, c].set_ylim(0, 2 * a.order)
        ax[1, c].set_title('convergence order')
        ax[1, c].legend(fontsize=8)
    # bottom row: central amplitude, lapse, AH radius / charge
    for i, r in enumerate(R):
        m = r['t'] <= tmax
        ax[2, 0].plot(r['t'][m], r['A'][m], color=colors[i], lw=1)
        ax[2, 1].plot(r['t'][m], r['lapse'][m], color=colors[i], lw=1)
        if len(keys) > 2:
            ah = np.where(r['ah'][m] > 0, r['ah'][m], np.nan)
            ax[2, 2].plot(r['t'][m], ah, color=colors[i], lw=1)
    ax[2, 0].set_title('central |phi|')
    ax[2, 1].set_title('central lapse')
    if len(keys) > 2:
        ax[2, 2].set_title('apparent horizon radius')
    for x in ax[2]:
        x.set_xlabel('t')
    fig.suptitle(a.label)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
