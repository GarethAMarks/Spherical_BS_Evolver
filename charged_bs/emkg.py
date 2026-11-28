"""
Full Einstein-Maxwell-Klein-Gordon charged solitonic boson stars (static, spherical).

Units G = c = mu = 1, Gaussian Maxwell; action as in the notes:
  S = int sqrt(-g)[R/16pi - 1/2(|D phi|^2 + V) - F^2/16pi],  D = nabla - i q A,
  V = |phi|^2 (1 - 2|phi|^2/sigma0^2)^2.
Ansatz: ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2, phi = A(r) e^{-i omega t}, A_mu dx^mu = Phi dt.
Unknowns per node (scaled, all O(1)):
  a = A/v, b = A'/v (v = sigma0/sqrt2), mt = sigma0^2 m (Misner-Sharp, incl. EM),
  L = ln alpha (alpha(0)=1, renormalised afterwards), wt = (omega + q Phi)/sigma0 in the
  unnormalised time, Qt = sigma0^2 Q.
Discretisation: trapezoidal rule on a uniform grid r in [r0, rmax]; global Newton (sparse LU).
BCs: a(r0) = lam, b(r0) = 0, mt(r0) = 0, Qt(r0) = 0, L(r0) = 0;  b = -X kappa a at rmax.
The frequency is an output (w(r0) is free and fixed by the outer BC).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

PI = np.pi
NV = 6


class Model:
    def __init__(self, sigma0, q, rmax, h, r0=1e-3, mode='a0', grid=None):
        self.mode = mode
        self.s0 = sigma0
        self.q = q
        self.v = sigma0 / np.sqrt(2)
        self.r = np.concatenate([[r0], np.arange(h, rmax + 1e-12, h)]) if grid is None else grid
        self.N = len(self.r)
        self.h = np.diff(self.r)

    # ---------------- right-hand side ----------------
    def F(self, r, Y):
        s0, q, v = self.s0, self.q, self.v
        a, b, mt, L, wt, Qt = Y.T
        m = mt / s0**2
        Q = Qt / s0**2
        X2 = 1.0 / (1.0 - 2.0 * m / r)
        X = np.sqrt(np.abs(X2))
        A, Ap = v * a, v * b
        Om2 = s0**2 * wt**2 * np.exp(-2 * L)
        V = v**2 * a**2 * (1 - a**2)**2
        dV = (1 - a**2) * (1 - 3 * a**2)
        kin = Om2 * A**2 + Ap**2 / X2
        rho = 0.5 * (kin + V)
        pr = 0.5 * (kin - V)
        E2 = Q**2 / r**4
        mp = 4 * PI * r**2 * rho + E2 * r**2 / 2
        Lp = X2 * (m + 4 * PI * r**3 * pr - Q**2 / (2 * r)) / r**2
        Qp = 4 * PI * r**2 * X * q * np.sqrt(Om2) * np.sign(wt) * A**2
        wp = q * np.exp(L) * X * Q / r**2
        XpX = X2 * (mp / r - m / r**2)
        App = X2 * (dV * A - Om2 * A) - Ap * (2 / r + Lp - XpX)
        return np.stack([Ap / v * 0 + b, App / v, s0**2 * mp, Lp, wp / s0, s0**2 * Qp], axis=1)

    # ---------------- residual ----------------
    def residual(self, Y, lam):
        r, h = self.r, self.h
        Fv = self.F(r, Y)
        R_int = (Y[1:] - Y[:-1]) / h[:, None] - 0.5 * (Fv[1:] + Fv[:-1])
        a, b, mt, L, wt, Qt = Y[0]
        g = (a - lam) if self.mode == 'a0' else ((wt - lam) if self.mode == 'w0' else self.pin(Y, lam))
        bc0 = np.array([g, b, mt, L, Qt])
        # outer Robin
        aN, bN, mN, LN, wN, QN = Y[-1]
        rN = r[-1]
        X2 = 1 / (1 - 2 * mN / self.s0**2 / rN)
        Om2 = self.s0**2 * wN**2 * np.exp(-2 * LN)
        kap = np.sqrt(max(1 - Om2, 1e-12))
        bc1 = np.array([bN + np.sqrt(X2) * kap * aN + aN / rN])
        return np.concatenate([bc0, R_int.ravel(), bc1])

    def pin(self, Y, Rw):
        k = int(np.searchsorted(self.r, Rw)) - 1
        th = (Rw - self.r[k]) / (self.r[k + 1] - self.r[k])
        return (1 - th) * Y[k, 0]**2 + th * Y[k + 1, 0]**2 - 0.5

    def pin_grad(self, Y, Rw):
        k = int(np.searchsorted(self.r, Rw)) - 1
        hk = self.r[k + 1] - self.r[k]
        th = (Rw - self.r[k]) / hk
        return k, 2 * (1 - th) * Y[k, 0], 2 * th * Y[k + 1, 0], (Y[k + 1, 0]**2 - Y[k, 0]**2) / hk

    # ---------------- sparse Jacobian by colouring ----------------
    def jacobian(self, Y, lam, R0=None):
        N = self.N
        if R0 is None:
            R0 = self.residual(Y, lam)
        rows, cols, vals = [], [], []
        nb0 = 5
        idx_int = np.arange(N - 1)
        for c in range(2):
            nodes = np.arange(c, N, 2)
            for k in range(NV):
                if k in (2, 5):
                    eps = 1e-7 * (np.abs(Y[nodes, k]) + 1e-4 * self.s0**2 * self.r[nodes]**3 + 1e-14)
                else:
                    eps = 1e-7 * (1 + np.abs(Y[nodes, k]))
                Yp = Y.copy()
                Yp[nodes, k] += eps
                d = self.residual(Yp, lam) - R0
                epsfull = np.zeros(N)
                epsfull[nodes] = eps
                # interval rows
                j = np.where(idx_int % 2 == c, idx_int, idx_int + 1)
                for e in range(NV):
                    rr = nb0 + idx_int * NV + e
                    rows.append(rr)
                    cols.append(j * NV + k)
                    vals.append(d[rr] / epsfull[j])

                if (N - 1) % 2 == c:
                    rr = np.array([len(R0) - 1])
                    rows.append(rr); cols.append(np.array([(N - 1) * NV + k])); vals.append(d[rr] / eps[-1])
        if self.mode == 'a0':
            k, g0, g1, gl, kv = 0, 1.0, 0.0, -1.0, 0
        elif self.mode == 'w0':
            k, g0, g1, gl, kv = 0, 1.0, 0.0, -1.0, 4
        else:
            (k, g0, g1, gl), kv = self.pin_grad(Y, lam), 0
        rows += [np.array([0, 0, 1, 2, 3, 4])]
        cols += [np.array([k * NV + kv, (k + 1) * NV + kv, 1, 2, 3, 5])]
        vals += [np.array([g0, g1, 1.0, 1.0, 1.0, 1.0])]
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(len(R0), N * NV))
        self._dlam = np.zeros(len(R0)); self._dlam[0] = gl
        return J

    # ---------------- Newton at fixed lam ----------------
    def safe(self, Y):
        return np.all(Y[:, 2] / self.s0**2 < 0.499 * self.r) and np.all(np.isfinite(Y))

    def solve(self, Y, lam, tol=1e-9, maxit=100, verbose=False, cap=None):
        for it in range(maxit):
            R0 = self.residual(Y, lam)
            nr = np.abs(R0).max()
            if verbose and it % 5 == 0:
                print(it, nr)
            if not np.isfinite(nr):
                return Y, False, it
            if nr < tol:
                return Y, True, it
            J = self.jacobian(Y, lam, R0)
            dY = spla.spsolve(J, -R0).reshape(Y.shape)
            n0 = np.linalg.norm(R0); t = 1.0
            while t > 1e-5:
                Yt = Y + t * dY
                if self.safe(Yt):
                    n1 = np.linalg.norm(self.residual(Yt, lam))
                    if np.isfinite(n1) and n1 < n0 * (1 - 1e-4 * t):
                        break
                t *= 0.5
            if t <= 1e-5:
                return Y, False, it
            Y = Yt
        return Y, False, maxit

    # ---------------- observables ----------------
    def observables(self, Y):
        s0, q, r = self.s0, self.q, self.r
        a, b, mt, L, wt, Qt = Y.T
        m, Q = mt / s0**2, Qt / s0**2
        rN = r[-1]
        M = m[-1] + Q[-1]**2 / (2 * rN)
        f = 1 - 2 * M / rN + Q[-1]**2 / rN**2
        ainf = np.exp(L[-1]) / np.sqrt(f)
        omega = s0 * wt[-1] / ainf + q * Q[-1] / rN
        # wall radius: where a^2 = 1/2 (outermost crossing)
        i = np.where(a**2 > 0.5)[0]
        if len(i) == 0:
            R = np.nan
        else:
            i = i[-1]
            R = r[i] + (0.5 - a[i]**2) / (a[i + 1]**2 - a[i]**2) * (r[i + 1] - r[i])
        # Noether charge
        X = 1 / np.sqrt(1 - 2 * m / r)
        Om = s0 * wt * np.exp(-L)
        n = 4 * PI * r**2 * X * Om * (self.v * a)**2
        Nn = np.sum(0.5 * (n[1:] + n[:-1]) * np.diff(r))
        # lapse at wall and centre (normalised)
        alpha = np.exp(L) / ainf
        alpha_c = alpha[0]
        fR = 1 - 2 * np.interp(R, r, m) / R if np.isfinite(R) else np.nan
        Rph = np.nan
        # outermost light ring: extremum of alpha^2/r^2
        H = alpha**2 / r**2
        dH = np.diff(H)
        ext = np.where(dH[:-1] * dH[1:] < 0)[0]
        nLR = len(ext)
        return dict(M=M, Q=Q[-1], N=Nn, omega=omega, R=R, C=M / R, QM=Q[-1] / M if M > 0 else np.nan,
                    alpha_c=alpha_c, a0=a[0], nLR=nLR, maxcomp=np.max(2 * m / r),
                    Mhat=s0**2 * M, Rhat=s0**2 * R, wedge=a[-1], Wc=wt[0] * np.exp(-L[0]),
                    Ws=np.interp(R, r, wt * np.exp(-L)) if np.isfinite(R) else np.nan,
                    alpha_R=np.interp(R, r, alpha) if np.isfinite(R) else np.nan,
                    fR_out=1 - 2 * M / R + Q[-1]**2 / R**2 if np.isfinite(R) else np.nan)


# ---------------- initial guess from the leading-order solution ----------------
def guess_from_LO(model, Wc, qh, pick='max', xmin=0.0):
    """Leading-order star at (W_c, qh) mapped onto the full grid."""
    import cbs_lo
    s0 = model.s0
    sols = cbs_lo.stars_for(Wc, qh, ngrid=1500)
    # choose the larger (compact) root
    sols = [p for p in sols if p['x'] > xmin]
    P = max(sols, key=lambda p: p['x']) if pick == 'max' else min(sols, key=lambda p: p['x'])
    sol = cbs_lo.core(Wc, qh)
    x_s = P['x']
    r = model.r
    x = s0**2 * r
    Rw = x_s / s0**2
    Y = np.zeros((model.N, NV))
    inside = x < x_s
    xi = np.clip(x, sol.t[0], x_s)
    mm, QQ, WW, LL = sol.sol(xi)
    # condensate amplitude incl. centrifugal shift: (1-s)(1-3s) = Omega^2
    Om2 = s0**2 * WW**2
    sA = (4 + np.sqrt(16 - 12 * (1 - Om2))) / 6
    prof = 1 / (1 + np.exp(2 * (r - Rw)))           # BPS wall: A^2 = v^2/(1+e^{2(r-R)})
    a = np.sqrt(np.where(inside, sA, 1.0) * prof)
    Y[:, 0] = a
    Y[:, 1] = np.gradient(a, r)
    # metric / charge: interior LO profile, exterior RN
    M, Q = P['M'], P['Q']
    m_ext = M - Q**2 / (2 * np.maximum(x, 1e-12))
    Y[:, 2] = np.where(inside, mm, m_ext)
    Y[:, 5] = np.where(inside, QQ, Q)
    # lapse: interior from LO core (alpha(0)=1 there already), exterior sqrt(f) rescaled to match
    up = P['up']
    Lin = LL - LL[0]
    L_s = (sol.sol(x_s)[3] - LL[0])
    f_ext = 1 - 2 * M / np.maximum(x, 1e-12) + Q**2 / np.maximum(x, 1e-12)**2
    Lext = 0.5 * np.log(np.maximum(f_ext, 1e-30)) - np.log(up) + L_s
    Y[:, 3] = np.where(inside, Lin, Lext)
    # w (unnormalised time): interior W*alpha_unnorm ; exterior from omega - q Q/r in normalised units
    ainf = np.exp(L_s) / up          # alpha_unnorm(inf)
    om = P['omega']
    w_in = WW * np.exp(Lin)
    w_ext = ainf * (om - qh * Q / np.maximum(x, 1e-12))
    Y[:, 4] = np.where(inside, w_in, w_ext)
    Y[0, 2] = 0; Y[0, 5] = 0; Y[0, 3] = 0
    return Y, P


def make_grid(Rw, hf, rmax, hc=0.1, kappa=0.04, r0=1e-3):
    """Graded grid: spacing hf at the wall Rw, growing linearly (rate kappa) to hc."""
    pts = [Rw]
    r = Rw
    while r < rmax:
        r += min(hc, hf + kappa * abs(r - Rw)); pts.append(r)
    r = Rw
    while True:
        r -= min(hc, hf + kappa * abs(r - Rw))
        if r <= 2 * hc: break
        pts.append(r)
    pts = np.sort(np.array(pts))
    first = np.arange(hc, pts[0] - 0.5 * hc, hc)
    return np.concatenate([[r0], first, pts])


def remesh(model_old, Y, grid):
    Ynew = np.empty((len(grid), NV))
    for k in range(NV):
        Ynew[:, k] = np.interp(grid, model_old.r, Y[:, k])
    # exterior beyond old rmax: extend with exponential tail / constants
    out = grid > model_old.r[-1]
    if np.any(out):
        Ynew[out, 0] = Y[-1, 0] * np.exp(-(grid[out] - model_old.r[-1]))
        Ynew[out, 1] = -Ynew[out, 0]
    return Ynew
