#!/usr/bin/env python3
"""
Charged thin-shell (solitonic) boson stars as initial data for SBSE.

Uses the static Einstein-Maxwell-Klein-Gordon solver in this directory (emkg.py / emkg2.py / arc2.py: pseudo-arclength
continuation in the scaled central local frequency Wc = Omega_c / sigma0), then re-integrates the metric and Maxwell
equations to high accuracy for the converged scalar profile and writes the star at the nodes of the SBSE evolution grid.

Units G = mu = 1, Gaussian Maxwell, V = |phi|^2 (1 - 2|phi|^2/sigma0^2)^2 (SBSE: solitonic = 1, sigma = sigma0, mu = 1,
lambda = 0, D = 4); gauge coupling q = qhat * sigma0.

Subcommands
-----------
  trace   trace a family at fixed (sigma0, qhat) and keep every profile:
            python3 make_charged_id.py trace --s0 0.1 --qhat 1.0 --rmax 100 --out fam.pkl
          (both directions from Wc^2 = --W2a; prints M, N, Q/M along the family and its N extrema)
  pairs   list same-N pairs (one member on each side of an N maximum) in a traced family:
            python3 make_charged_id.py pairs --family fam.pkl --nfrac 0.98
  write   write one star at the grid nodes of an SBSE run (R, n_gridpoints taken from BSParams.par unless given):
            python3 make_charged_id.py write --family fam.pkl --index 123 [--R 200 --n 4001] --out charged_BSdata.dat

Output file (read by BosonStar::read_charged_data): header lines "# key = value", then columns
  r   A   X   ln(alpha)   E^z   Phi_e   Omega
with A the scalar amplitude (phi = A e^{+i omega t} in SBSE's convention), X^2 = g_rr (polar-areal), alpha normalised
to 1 at infinity, E^z = Q(r) / (r^2 X), Phi_e = -A_t / alpha (A_t -> 0 at infinity), Omega = (omega + q A_t) / alpha.
The Hamiltonian and Gauss constraints are satisfied by construction up to the ODE tolerance (1e-12): only the scalar
profile A(r) comes from the second-order solver (trapezoidal, ~h^2 error), which affects how stationary the star is,
not the constraints.
"""
import argparse
import os
import pickle
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PI = np.pi


# ----------------------------------------------------------------------------------------------------------------------
# family tracing
# ----------------------------------------------------------------------------------------------------------------------
def trace(s0, qh, rmax, W2a=25.0, fine=False, nsteps=4000, both=True, verbose=True):
    import arc2
    kw = dict(hc=0.05, fac=0.04, hmax=0.02) if fine else {}
    fam = []

    def run(direction):
        rows = []

        def snap(o, mm, pt):
            Y, D, lam = pt
            rows.append(dict(obs=o, r=mm.r_of(D).copy(), Y=Y.copy(), lam=lam))

        Mmax = [0.0]

        def stop(o, out):
            # stop well past the mass maximum, or at a near-extremal end
            Mmax[0] = max(Mmax[0], o['Mhat'])
            return o['Mhat'] < 0.9 * Mmax[0] or arc2.fmin_ext(o) < 1e-4

        arc2.arc(s0, qh, rmax, W2a=W2a, nsteps=nsteps, stop=stop, verbose=verbose, save=False, direction=direction,
                 snap=snap, **kw)
        return rows

    fwd = run(+1)
    bwd = run(-1) if both else []
    fam = bwd[::-1] + fwd
    return dict(s0=s0, qh=qh, rmax=rmax, fine=fine, rows=fam)


def family_table(fam):
    o = [r['obs'] for r in fam['rows']]
    s0 = fam['s0']
    Nhat = np.array([x['N'] for x in o]) * s0**3
    return dict(Wc2=np.array([x['Wc2'] for x in o]), Mhat=np.array([x['Mhat'] for x in o]),
                Nhat=Nhat, QM=np.array([x['QM'] for x in o]), C=np.array([x['C'] for x in o]),
                M=np.array([x['M'] for x in o]), R=np.array([x['R'] for x in o]),
                alpha_c=np.array([x['alpha_c'] for x in o]))


def extrema(y):
    d = np.diff(y)
    return [i + 1 for i in range(len(d) - 1) if d[i] * d[i + 1] < 0]


# ----------------------------------------------------------------------------------------------------------------------
# high-accuracy re-integration of the metric + Maxwell field for a converged scalar profile
# ----------------------------------------------------------------------------------------------------------------------
def build_star(s0, q, r_sol, Y, r_nodes, rtol=1e-12, max_step=0.05):
    """Given the solver's profile (r_sol, Y = [a, b, mt, L, wt, Qt]) return dict of node arrays for SBSE.
    A(r) is represented by an even quintic spline through the solver values; m, ln alpha, omega + q A_t, Q (and the
    Noether charge) are integrated with DOP853 for that A(r)."""
    v = s0 / np.sqrt(2.0)
    A_sol = v * Y[:, 0]
    w_c = s0 * Y[0, 4]                       # central omega + q A_t (unnormalised time, alpha(0) = 1)
    # even extension about r = 0 and exponential tail beyond the solver's outer radius
    rN = r_sol[-1]
    Om2N = s0**2 * Y[-1, 4]**2 * np.exp(-2 * Y[-1, 3])
    kap = np.sqrt(max(1 - Om2N, 1e-12))
    r_out = max(r_nodes[-1], rN) + 1.0
    r_tail = np.arange(rN + 0.5, r_out + 0.5, 0.5)
    A_tail = A_sol[-1] * (rN / r_tail) * np.exp(-kap * (r_tail - rN))
    rr = np.concatenate([-r_sol[::-1], r_sol, r_tail])
    AA = np.concatenate([A_sol[::-1], A_sol, A_tail])
    keep = np.concatenate([[True], np.diff(rr) > 1e-10])
    splA = make_interp_spline(rr[keep], AA[keep], k=5)
    dsplA = splA.derivative()

    def rhs(r, y):
        m, L, w, Q, N = y
        A = splA(r); Ap = dsplA(r)
        X2 = 1.0 / (1.0 - 2.0 * m / r)
        X = np.sqrt(X2)
        Om = w * np.exp(-L)
        V = A**2 * (1 - 2 * A**2 / s0**2)**2
        kin = Om**2 * A**2 + Ap**2 / X2
        rho = 0.5 * (kin + V)
        p = 0.5 * (kin - V)
        return [4 * PI * r**2 * rho + Q**2 / (2 * r**2),
                X2 * (m + 4 * PI * r**3 * p - Q**2 / (2 * r)) / r**2,
                q * np.exp(L) * X * Q / r**2,
                4 * PI * r**2 * X * q * Om * A**2,
                4 * PI * r**2 * X * Om * A**2]

    # series start at small r
    A_c = float(splA(0.0)); Om_c = w_c
    V_c = A_c**2 * (1 - 2 * A_c**2 / s0**2)**2
    rho_c = 0.5 * (Om_c**2 * A_c**2 + V_c); p_c = 0.5 * (Om_c**2 * A_c**2 - V_c)
    r1 = 1e-4
    y1 = [4 * PI / 3 * rho_c * r1**3, 2 * PI / 3 * (rho_c + 3 * p_c) * r1**2,
          w_c + 2 * PI / 3 * q * q * Om_c * A_c**2 * r1**2,
          4 * PI / 3 * q * Om_c * A_c**2 * r1**3, 4 * PI / 3 * Om_c * A_c**2 * r1**3]
    t_eval = r_nodes[r_nodes > r1]
    sol = solve_ivp(rhs, (r1, r_nodes[-1]), y1, method='DOP853', rtol=rtol, atol=1e-14, t_eval=t_eval,
                    max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    nz = len(r_nodes) - len(t_eval)                 # nodes at r <= r1 (the origin)
    m = np.concatenate([np.zeros(nz), sol.y[0]])
    L = np.concatenate([np.zeros(nz), sol.y[1]])
    w = np.concatenate([np.full(nz, w_c), sol.y[2]])
    Q = np.concatenate([np.zeros(nz), sol.y[3]])
    Nr = np.concatenate([np.zeros(nz), sol.y[4]])
    r = r_nodes
    A = splA(r)
    X = 1.0 / np.sqrt(1.0 - 2.0 * np.where(r > 0, m / np.maximum(r, 1e-300), 0.0))
    # normalise to Reissner-Nordstrom at the outer node
    rE = r[-1]
    M = m[-1] + Q[-1]**2 / (2 * rE)
    f = 1 - 2 * M / rE + Q[-1]**2 / rE**2
    ainf = np.exp(L[-1]) / np.sqrt(f)
    lnalpha = L - np.log(ainf)
    alpha = np.exp(lnalpha)
    omega = w[-1] / ainf + q * Q[-1] / rE
    Omega = w * np.exp(-L)                          # (omega + q A_t) / alpha, normalisation independent
    A_t = (w / ainf - omega) / q if q != 0 else np.zeros_like(r)
    Phi_e = -A_t / alpha
    Ez = np.where(r > 0, Q / np.maximum(r, 1e-300)**2 / X, 0.0)
    Ntot = Nr[-1]
    # radii: wall (A^2 = sigma0^2/4, i.e. a^2 = 1/2), 99% of N, and the constraint-norm radius (N within 1e-6 of total)
    a2 = (A / v)**2
    iw = np.where(a2 > 0.5)[0]
    R_wall = r[iw[-1]] if len(iw) else np.nan
    r_99 = float(np.interp(0.99 * Ntot, Nr, r))
    r_norm = float(np.interp((1 - 1e-6) * Ntot, Nr, r))
    return dict(r=r, A=A, X=X, lnalpha=lnalpha, Ez=Ez, Phi_e=Phi_e, Omega=Omega,
                meta=dict(sigma=s0, mu=1.0, q=q, qhat=q / s0, omega=omega, M=M, Q=Q[-1], N=Ntot, QoverM=Q[-1] / M,
                          R_wall=R_wall, r_99=r_99, r_norm=r_norm, C=M / R_wall, alpha_c=alpha[0], Wc=Y[0, 4],
                          Mhat=s0**2 * M, Nhat=s0**3 * Ntot))


def write_star(star, fname, extra=None):
    meta = dict(star['meta'])
    if extra:
        meta.update(extra)
    with open(fname, 'w') as fh:
        fh.write('# charged solitonic boson star for SBSE (charged = 1); written by charged_bs/make_charged_id.py\n')
        for k, val in meta.items():
            fh.write('# %s = %.16g\n' % (k, val))
        fh.write('# columns: r A X ln(alpha) E^z Phi_e Omega\n')
        dat = np.column_stack([star['r'], star['A'], star['X'], star['lnalpha'], star['Ez'], star['Phi_e'], star['Omega']])
        np.savetxt(fh, dat, fmt='%.17g')


def read_par(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split('#')[0]
            if '=' in line:
                k, val = line.split('=', 1)
                out[k.strip()] = val.strip()
    return out


def refine_row(fam, row, tol=1e-11):
    """Re-converge a stored profile at its own parameter with a tight tolerance (returns r, Y)."""
    from emkg2 import MM, make_xi
    s0, qh, rmax = fam['s0'], fam['qh'], fam['rmax']
    o = row['obs']
    kw = dict(hc=0.05, fac=0.04, hmax=0.02) if fam['fine'] else dict(hc=0.1, fac=0.08, hmax=0.03)
    u = np.sqrt(max(o['fR_out'], 1e-12))
    hf = min(kw['hmax'], kw['fac'] * u)
    xi, kwall = make_xi(o['R'], hf, rmax, hc=kw['hc'])
    mm = MM(s0, qh * s0, xi, kwall)
    Y0 = np.stack([np.interp(xi, row['r'], row['Y'][:, k]) for k in range(row['Y'].shape[1])], axis=1)
    Y, D, ok, _ = mm.solve(Y0, 0.0, row['lam'], tol=tol)
    if not ok:
        print('  tight re-solve failed; using the stored profile')
        return row['r'], row['Y']
    return mm.r_of(D), Y


# ----------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('trace')
    t.add_argument('--s0', type=float, required=True)
    t.add_argument('--qhat', type=float, required=True)
    t.add_argument('--rmax', type=float, required=True)
    t.add_argument('--W2a', type=float, default=25.0)
    t.add_argument('--fine', action='store_true')
    t.add_argument('--one-way', action='store_true', help='only trace towards smaller Wc (more compact)')
    t.add_argument('--out', required=True)
    p = sub.add_parser('pairs')
    p.add_argument('--family', required=True)
    p.add_argument('--nfrac', type=float, nargs='+', default=[0.99, 0.97, 0.95])
    w = sub.add_parser('write')
    w.add_argument('--family', required=True)
    w.add_argument('--index', type=int, required=True)
    w.add_argument('--par', default='BSParams.par')
    w.add_argument('--R', type=float)
    w.add_argument('--n', type=int)
    w.add_argument('--no-refine', action='store_true')
    w.add_argument('--out', default='charged_BSdata.dat')
    a = ap.parse_args()

    if a.cmd == 'trace':
        fam = trace(a.s0, a.qhat, a.rmax, W2a=a.W2a, fine=a.fine, both=not a.one_way)
        pickle.dump(fam, open(a.out, 'wb'))
        T = family_table(fam)
        print('%5s %10s %10s %10s %10s %8s %8s' % ('idx', 'Wc2', 'Mhat', 'Nhat', 'Q/M', 'C', 'alpha_c'))
        for i in range(0, len(T['Wc2']), max(1, len(T['Wc2']) // 60)):
            print('%5d %10.5f %10.6f %10.6f %10.6f %8.4f %8.4f' % (i, T['Wc2'][i], T['Mhat'][i], T['Nhat'][i], T['QM'][i],
                                                                 T['C'][i], T['alpha_c'][i]))
        print('N extrema at indices', extrema(T['Nhat']), ' M extrema at', extrema(T['Mhat']))
    elif a.cmd == 'pairs':
        fam = pickle.load(open(a.family, 'rb'))
        T = family_table(fam)
        N = T['Nhat']
        for k in extrema(N):
            for fr in a.nfrac:
                target = fr * N[k]
                lo = [i for i in range(1, k + 1) if (N[i - 1] - target) * (N[i] - target) <= 0]
                hi = [i for i in range(k, len(N)) if (N[i - 1] - target) * (N[i] - target) <= 0]
                if lo and hi:
                    i, j = lo[-1], hi[0]
                    print('extremum %d (Nhat=%.6f): N/Next=%.3f  side A idx %d (Wc2=%.4f C=%.4f Q/M=%.5f)  side B idx %d '
                          '(Wc2=%.4f C=%.4f Q/M=%.5f)' % (k, N[k], fr, i, T['Wc2'][i], T['C'][i], T['QM'][i],
                                                          j, T['Wc2'][j], T['C'][j], T['QM'][j]))
    elif a.cmd == 'write':
        fam = pickle.load(open(a.family, 'rb'))
        row = fam['rows'][a.index]
        par = read_par(a.par) if os.path.exists(a.par) else {}
        R = a.R if a.R else float(par['R'])
        n = a.n if a.n else int(par['n_gridpoints'])
        if par:
            s_par = float(par.get('sigma', fam['s0']))
            if abs(s_par - fam['s0']) > 1e-12 or par.get('solitonic', '1') != '1':
                print('WARNING: BSParams.par has sigma = %s, solitonic = %s; the star has sigma0 = %g (SBSE will use the '
                      'star\'s values)' % (par.get('sigma'), par.get('solitonic'), fam['s0']))
        r_sol, Y = (row['r'], row['Y']) if a.no_refine else refine_row(fam, row)
        r_nodes = np.linspace(0.0, R, n)
        star = build_star(fam['s0'], fam['qh'] * fam['s0'], r_sol, Y, r_nodes)
        write_star(star, a.out, extra=dict(family_index=a.index, Wc2=row['lam']**2))
        mt = star['meta']
        print('wrote %s: sigma0=%g qhat=%g Wc2=%.5f  M=%.6g Q=%.6g Q/M=%.6f N=%.6g omega=%.6g R_wall=%.3f C=%.4f '
              'alpha_c=%.4g r_norm=%.2f' % (a.out, fam['s0'], fam['qh'], row['lam']**2, mt['M'], mt['Q'], mt['QoverM'],
                                           mt['N'], mt['omega'], mt['R_wall'], mt['C'], mt['alpha_c'], mt['r_norm']))
        o = row['obs']
        print('   solver values:  M=%.6g Q=%.6g N=%.6g omega=%.6g' % (o['M'], o['Q'], o['N'], o['omega']))


if __name__ == '__main__':
    main()
