"""
Relaxation (collocation) solver for spherical Dirac stars, for the regimes where shooting fails (thin walls, ultracompact
stars of the degenerate-vacuum potential). Same physics and conventions as DiracStar (DiracStar.h, DiracPotential.h),
signature (-,+,+,+), G = 1:

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,   alpha = exp(phi),   X^2 = 1 / (1 - 2 m / r),
    two spinors psi_m = e^{-i omega t} ( g(r) Omega_{-1,m}, i f(r) Omega_{+1,m} )   (f, g the physical amplitudes of DiracStar:
                                                                                      F = sqrt(alpha) r f, G = sqrt(alpha) r g),
    S = (g^2 - f^2) / (2 pi),   n = (g^2 + f^2) / (2 pi),   mu_eff = mu + U'(S),   P = S U' - U,
    U(S) = (lambda/2) S^2 + [solitonic] mu S ((1 - S/S0)^2 - 1),   S0 = sigma^2 / (4 pi)   (W = mu S + U = mu S (1 - S/S0)^2).

The first-order system (from DiracStar's equations for F, G with F = sqrt(alpha) r f, G = sqrt(alpha) r g):
    m'     = 4 pi r^2 rho,                       rho = (omega/alpha) n - P(S)
    phi'   = (X^2 - 1) / (2 r) + 4 pi r X^2 p_r,     p_r = (omega/alpha) n - mu S + f g / (pi r) - U(S)
    g'     = (X - 1) g / r - phi' g / 2 + X (mu_eff + omega/alpha) f
    f'     = -(1 + X) f / r - phi' f / 2 + X (mu_eff - omega/alpha) g
    omega' = 0                                   (the frequency promoted to a grid variable)
on r in [r0, R], with the boundary conditions
    g(r0) = A,   f(r0) = (mu_eff(A^2/2pi) - omega/alpha(r0)) A r0 / 3,   m(r0) = 0          (regularity),
    phi(R) = ln(1 - 2 m(R)/R) / 2                                                           (alpha -> 1 at infinity),
    f(R) + sqrt((mu - omega/alpha)/(mu + omega/alpha)) g(R) = 0                              (exponential decay),
and one more condition fixing the point on the family: A = A_target, or a pseudo-arclength condition in (A, omega) with A
an unknown parameter, which follows the family through turning points in A. scipy's solve_bvp (4th-order collocation with
adaptive mesh refinement) is used, which resolves thin walls automatically.

Diagnostics: M = m(R), N = int 2 r^2 X (f^2 + g^2) dr (both spinors), r_99 (99% of M), C = M / r_99, C_max = max m/r, the
light rings (zeros of d/dr (alpha / r), i.e. r phi' = 1), C4 = M / r4 (BosonStar::get_c4), omega / alpha(0), the energy conditions (wec = min rho / max rho, negative if the
weak energy condition fails; dec = max(|p_r|, |p_t|) - rho, over max rho, positive if the dominant one fails), and the
first-law residual along a family (dM = omega dN).
"""
import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid


class Potential:
    def __init__(self, mu=1., lam=0., solitonic=True, sigma=0.2):
        self.mu, self.lam, self.solitonic, self.sigma = mu, lam, solitonic, sigma
        self.S0 = sigma**2 / (4 * np.pi)

    def U(self, S):
        u = 0.5 * self.lam * S**2
        if self.solitonic:
            y = S / self.S0
            u = u + self.mu * S * (y * y - 2 * y)
        return u

    def dU(self, S):
        du = self.lam * S
        if self.solitonic:
            y = S / self.S0
            du = du + self.mu * (3 * y * y - 4 * y)
        return du

    def ddU(self, S):
        ddu = self.lam + 0. * S
        if self.solitonic:
            y = S / self.S0
            ddu = ddu + self.mu * (6 * y - 4) / self.S0
        return ddu

    def mu_eff(self, S):
        return self.mu + self.dU(S)

    def P(self, S):
        return S * self.dU(S) - self.U(S)


class DiracRelax:
    def __init__(self, pot, r0=1e-4, R=60., tol=1e-6, max_nodes=300000):
        self.pot, self.r0, self.R, self.tol, self.max_nodes = pot, r0, R, tol, max_nodes

    # --- equations ---
    def _terms(self, r, y):
        f, g, m, phi, w = y
        pot = self.pot
        X = 1. / np.sqrt(np.maximum(1. - 2. * m / r, 1e-300))
        wa = w * np.exp(-phi)
        S = (g * g - f * f) / (2 * np.pi)
        n = (g * g + f * f) / (2 * np.pi)
        rho = wa * n - pot.P(S)
        pr = wa * n - pot.mu * S + f * g / (np.pi * r) - pot.U(S)
        dphi = (X * X - 1.) / (2. * r) + 4. * np.pi * r * X * X * pr
        return X, wa, S, rho, pr, dphi

    def rhs(self, r, y, p=None):
        f, g, m, phi, w = y
        X, wa, S, rho, pr, dphi = self._terms(r, y)
        me = self.pot.mu_eff(S)
        dm = 4. * np.pi * r * r * rho
        dg = (X - 1.) * g / r - 0.5 * dphi * g + X * (me + wa) * f
        df = -(1. + X) * f / r - 0.5 * dphi * f + X * (me - wa) * g
        return np.vstack([df, dg, dm, dphi, np.zeros_like(r)])

    def _bc_common(self, ya, yb, A):
        pot, r0, R = self.pot, self.r0, self.R
        S0c = A * A / (2 * np.pi)
        wa0 = ya[4] * np.exp(-ya[3])
        alR = np.exp(yb[3])
        waR = yb[4] / alR
        ratio = np.sqrt(np.maximum(pot.mu - waR, 0.) / (pot.mu + waR))
        return [ya[1] - A,
                ya[0] - (pot.mu_eff(S0c) - wa0) * A * r0 / 3.,
                ya[2],
                yb[3] - 0.5 * np.log(1. - 2. * yb[2] / R),
                yb[0] + ratio * yb[1]]

    # --- solves ---
    def solve_fixed_A(self, A, r, y, verbose=0):
        """solve at fixed central amplitude A, from the initial guess y on the mesh r"""
        bc = lambda ya, yb: np.array(self._bc_common(ya, yb, A))
        sol = solve_bvp(self.rhs, bc, r, y, tol=self.tol, max_nodes=self.max_nodes, verbose=verbose)
        return sol, A

    def solve_arclength(self, r, y, A_guess, ref, tangent, ds, verbose=0):
        """pseudo-arclength step: unknown A (parameter), constraint
            tA (A - A_ref) / sA + tw (w - w_ref) / sw + tM (ln M - ln M_ref) = ds,
        with ref = (A_ref, w_ref, lnM_ref, sA, sw) and the unit tangent (tA, tw, tM); M = m(R). Including ln M keeps the
        steps small where the mass changes quickly at almost fixed A and omega (thin-wall stars)"""
        A_ref, w_ref, lnM_ref, sA, sw = ref
        tA, tw, tM = tangent
        def bc(ya, yb, p):
            A = p[0]
            res = self._bc_common(ya, yb, A)
            res.append(tA * (A - A_ref) / sA + tw * (ya[4] - w_ref) / sw + tM * (np.log(max(yb[2], 1e-300)) - lnM_ref) - ds)
            return np.array(res)
        rhs = lambda rr, yy, p: self.rhs(rr, yy)
        sol = solve_bvp(rhs, bc, r, y, p=[A_guess], tol=self.tol, max_nodes=self.max_nodes, verbose=verbose)
        return sol, (sol.p[0] if sol.p is not None else A_guess)

    # --- diagnostics ---
    def diagnostics(self, sol, A):
        r = np.concatenate([[0.], sol.x])
        y = sol.y
        f, g, m, phi, w = y
        X, wa, S, rho, pr, dphi = self._terms(sol.x, y)
        M = m[-1]
        N = np.trapezoid(2 * sol.x**2 * X * (f * f + g * g), sol.x)
        # radius containing 99% of the mass (linear interpolation)
        j = np.searchsorted(m, 0.99 * M)
        r99 = np.interp(0.99 * M, m[max(j - 1, 0):j + 1], sol.x[max(j - 1, 0):j + 1]) if 0 < j < len(m) else np.nan
        cm = m / sol.x
        # C4 (Collodel et al., as BosonStar::get_c4): M / r4, r4 the outermost radius where X^2 deviates from its
        # Schwarzschild value 1 / (1 - 2M/r) by more than 1e-3 (relative), on a fine uniform grid
        rf = np.linspace(sol.x[0], sol.x[-1], 200001)[1:-1]
        mf = sol.sol(rf)[2]
        dev = np.abs((1. - 2. * M / rf) / (1. - 2. * mf / rf) - 1.)
        out = np.where(dev > 1e-3)[0]
        if len(out) and out[-1] + 1 < len(rf):
            j4 = out[-1]
            r4 = rf[j4] + (rf[j4 + 1] - rf[j4]) * (1e-3 - dev[j4]) / (dev[j4 + 1] - dev[j4])
        else:
            r4 = np.nan
        # light rings: r phi' = 1 (stationary points of alpha / r)
        h = sol.x * dphi - 1.
        lr = sol.x[1:][np.sign(h[1:]) != np.sign(h[:-1])]
        # energy conditions: weak (rho >= 0) and dominant (|p_r|, |p_t| <= rho), normalised by the maximum of rho
        pt = -f * g / (2 * np.pi * sol.x) + self.pot.P(S)
        rmax = rho.max()
        wec = rho.min() / rmax
        dec = (np.maximum(np.abs(pr), np.abs(pt)) - rho).max() / rmax
        return dict(A=A, omega=w[0], M=M, N=N, Eb=M - self.pot.mu * N, r99=r99, C=M / r99, Cmax=cm.max(),
                    r_Cmax=sol.x[cm.argmax()], phi0=phi[0], alpha0=np.exp(phi[0]), light_rings=lr,
                    nodes=len(sol.x), g_nodes=int(np.sum(np.sign(g[1:]) != np.sign(g[:-1]))), wec=wec, dec=dec, C4=M / r4, omega_c=w[0] / np.exp(phi[0]))


def seed_from_DSdata(path, r_mesh):
    """initial guess from DiracStar's DSdata.dat (r F G X phi m f g) and its frequency"""
    d = np.genfromtxt(path)
    d = d[np.isfinite(d).all(axis=1)]
    r = d[:, 0]
    y = np.vstack([np.interp(r_mesh, r, d[:, 6]), np.interp(r_mesh, r, d[:, 7]), np.interp(r_mesh, r, d[:, 5]),
                   np.interp(r_mesh, r, d[:, 4]), np.zeros_like(r_mesh)])
    return y


def remesh(sol, R_new=None, n_min=400):
    """mesh and guess for the next solve: the previous mesh (optionally extended to R_new with the exponential tail)"""
    x, y = sol.x, sol.y.copy()
    if len(x) < n_min:
        xn = np.unique(np.concatenate([x, np.linspace(x[0], x[-1], n_min)]))
        y = sol.sol(xn)
        x = xn
    return x, y


def equidistribute(sol, n=3000):
    """a fresh mesh of n nodes equidistributing 1 + |d profiles/dr| (arc length of g, f, m, phi, each normalised by its
    maximum, plus a uniform part), so that thin walls and the core stay resolved while the mesh does not grow from step to
    step; returns the mesh and the solution interpolated onto it"""
    x = sol.x
    y = sol.sol(x)
    sc = [max(np.abs(y[i]).max(), 1e-30) for i in range(4)]
    d = np.sqrt(sum(np.gradient(y[i] / sc[i], x)**2 for i in range(4)))
    mon = 1. + d * (x[-1] - x[0]) / 20.
    cum = np.concatenate([[0.], np.cumsum(0.5 * (mon[1:] + mon[:-1]) * np.diff(x))])
    xn = np.interp(np.linspace(0, cum[-1], n), cum, x)
    xn[0], xn[-1] = x[0], x[-1]
    xn = np.unique(xn)
    return xn, sol.sol(xn)


def extend(sol, R_new, s, n_tail=400):
    """mesh and guess on [r0, R_new] from a solution on [r0, R_old]: truncated, or extended with the exponential tail
    (f, g ~ exp(-k r), m = M, phi from Schwarzschild)"""
    x, y = sol.x, sol.y
    if R_new <= x[-1]:
        keep = x < R_new
        xn = np.concatenate([x[keep], [R_new]])
        return xn, sol.sol(xn)
    f, g, m, phi, w = y[:, -1]
    M = m
    al = np.exp(phi)
    k = np.sqrt(max(s.pot.mu**2 - (w / al)**2, 1e-12)) / np.sqrt(1 - 2 * M / x[-1])
    xt = np.linspace(x[-1], R_new, n_tail)[1:]
    dec = np.exp(-k * (xt - x[-1]))
    yt = np.vstack([f * dec, g * dec, np.full_like(xt, M), 0.5 * np.log(1 - 2 * M / xt), np.full_like(xt, w)])
    return np.concatenate([x, xt]), np.hstack([y, yt])


def continue_family(s, sol, A, n_steps=200, ds=0.02, direction=+1, scales=None, R_rule=None, stop=None, verbose=True,
                    ds_min=1e-6, ds_max=0.2, n_mesh=3000):
    """pseudo-arclength continuation of a family in the (A / sA, omega / sw, ln M) space, starting from sol at amplitude A.
    scales = (sA, sw) normalise A and omega; the first tangent is +-(1, 0) in A (direction) unless a second point is
    known. R_rule(diag) -> outer radius to use for the next step. stop(diag) -> True ends the continuation.
    Returns the list of diagnostics (with 'sol' for each point)."""
    sA, sw = scales if scales is not None else (A, 1.)
    out = []
    d = s.diagnostics(sol, A); d['sol'] = sol; out.append(d)
    prev = None
    ds_prev = ds
    tangent = np.array([direction, 0., 0.])
    while len(out) < n_steps:
        cur = out[-1]
        if prev is not None:
            t = np.array([(cur['A'] - prev['A']) / sA, (cur['omega'] - prev['omega']) / sw, np.log(cur['M'] / prev['M'])])
            if np.linalg.norm(t) > 0:
                tangent = t / np.linalg.norm(t)
        R_new = R_rule(cur) if R_rule is not None else s.R
        x, y = extend(cur['sol'], R_new, s)
        s.R = R_new
        if prev is not None:
            # secant predictor in the radius scaled by r_99: y(rho) = y_cur(rho) + q (y_cur(rho) - y_prev(rho)), rho = r / r_99,
            # with q the ratio of this step to the previous one (thin walls move with the stellar radius)
            q = ds / max(ds_prev, 1e-12)
            r99p = cur['r99'] + q * (cur['r99'] - prev['r99'])
            if r99p > 0:
                rho = x / r99p
                yc = cur['sol'].sol(np.clip(rho * cur['r99'], x[0], cur['sol'].x[-1]))
                yp = prev['sol'].sol(np.clip(rho * prev['r99'], x[0], prev['sol'].x[-1]))
                y = yc + q * (yc - yp)
                y[4] = cur['omega'] + q * (cur['omega'] - prev['omega'])
        if True:   # always re-mesh, so that the mesh follows the solution along the family without growing
            class _Sol: pass
            tmp = _Sol(); tmp.x = x
            tmp.sol = lambda xx, xs=x, ys=y: np.vstack([np.interp(xx, xs, yy) for yy in ys])
            x, y = equidistribute(tmp, n_mesh)
        A_pred = cur['A'] + ds * tangent[0] * sA
        y = y.copy()
        if prev is None:
            # first step: scale g, f by the amplitude ratio
            y[4] = cur['omega'] + ds * tangent[1] * sw
            y[:2] *= A_pred / cur['A']
        sol_new, A_new = s.solve_arclength(x, y, A_pred, (cur['A'], cur['omega'], np.log(cur['M']), sA, sw), tangent, ds)
        ok = sol_new.status == 0
        if ok:
            dn = s.diagnostics(sol_new, A_new)
            # reject jumps to other branches or excited states
            if (dn['g_nodes'] != 0 or abs(dn['A'] - cur['A']) / sA > 3 * ds or abs(dn['omega'] - cur['omega']) / sw > 3 * ds
                    or abs(np.log(dn['M'] / cur['M'])) > 3 * ds):
                ok = False
        if not ok:
            ds *= 0.5
            if verbose:
                print('   step failed (%s), ds -> %.2e' % (sol_new.message[:40], ds))
            if ds < ds_min:
                print('   ds below minimum, stopping')
                break
            continue
        dn['sol'] = sol_new
        prev = cur
        ds_prev = ds
        out.append(dn)
        if verbose:
            print('A %.6f  omega %.6f  M %.5f  N %.5f  C %.4f  Cmax %.4f  LR %d  nodes %d  ds %.3g' % (
                dn['A'], dn['omega'], dn['M'], dn['N'], dn['C'], dn['Cmax'], len(dn['light_rings']), dn['nodes'], ds))
        if stop is not None and stop(dn):
            break
        if sol_new.niter <= 4:
            ds = min(ds * 1.5, ds_max)
    return out
