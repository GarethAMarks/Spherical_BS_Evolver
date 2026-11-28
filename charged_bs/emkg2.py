"""
Moving-mesh version of the EMKG solver (see emkg.py for equations/conventions).

Computational grid xi_i with the wall pinned at node kw (a^2 = 1/2 there).
Physical radius r_i = xi_i + D * phi(xi_i), with phi a smooth bump (=1 near the wall,
0 near the origin and far away) and D an extra unknown: wall motion = change of D.
Parameter lam = central local frequency (scaled), BC wt(r0) = lam.
Unknown vector: U = (Y.ravel(), D);  equations: 6N residuals + pin.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from emkg import Model, NV, PI


def bump(xi, xw, a1=3.0, a2=6.0):
    d = np.abs(xi - xw)
    phi = np.where(d <= a1, 1.0, np.where(d >= a1 + a2, 0.0, 0.5 * (1 + np.cos(PI * (d - a1) / a2))))
    return phi


def make_xi(xw, hf, rmax, hc=0.1, kappa=0.04, r0=1e-3):
    pts = [xw]; r = xw
    while r < rmax:
        r += min(hc, hf + kappa * abs(r - xw)); pts.append(r)
    r = xw
    while True:
        r -= min(hc, hf + kappa * abs(r - xw))
        if r <= 2 * hc: break
        pts.append(r)
    pts = np.sort(np.array(pts))
    first = np.arange(hc, pts[0] - 0.5 * hc, hc)
    xi = np.concatenate([[r0], first, pts])
    kw = int(np.argmin(np.abs(xi - xw)))
    return xi, kw


class MM:
    def __init__(self, s0, q, xi, kw, a1=3.0, a2=6.0):
        self.base = Model(s0, q, xi[-1], 0.1, grid=xi.copy())
        self.s0, self.q = s0, q
        self.xi, self.kw = xi, kw
        self.phi = bump(xi, xi[kw], a1, a2)
        self.phi[0] = 0.0
        self.N = len(xi)
        self.dphi_max = PI / (2 * a2)

    def r_of(self, D):
        return self.xi + D * self.phi

    def _set(self, D):
        r = self.r_of(D)
        self.base.r = r
        self.base.h = np.diff(r)
        return r

    def residual(self, Y, D, lam):
        self._set(D)
        b = self.base
        Fv = b.F(b.r, Y)
        R_int = (Y[1:] - Y[:-1]) - 0.5 * b.h[:, None] * (Fv[1:] + Fv[:-1])
        a, bb, mt, L, wt, Qt = Y[0]
        bc0 = np.array([wt - lam, bb, mt, L, Qt])
        aN, bN, mN, LN, wN, QN = Y[-1]
        rN = b.r[-1]
        X2 = 1 / (1 - 2 * mN / self.s0**2 / rN)
        Om2 = self.s0**2 * wN**2 * np.exp(-2 * LN)
        kap = np.sqrt(max(1 - Om2, 1e-12))
        bc1 = np.array([bN + np.sqrt(X2) * kap * aN + aN / rN])
        pin = np.array([Y[self.kw, 0]**2 - 0.5])
        return np.concatenate([bc0, R_int.ravel(), bc1, pin])

    def jacobian(self, Y, D, lam, R0=None):
        N = self.N
        if R0 is None:
            R0 = self.residual(Y, D, lam)
        nrow = len(R0)
        rows, cols, vals = [], [], []
        nb0 = 5
        idx = np.arange(N - 1)
        r = self._set(D)
        for c in range(2):
            nodes = np.arange(c, N, 2)
            for k in range(NV):
                if k in (2, 5):
                    eps = 1e-7 * (np.abs(Y[nodes, k]) + 1e-4 * self.s0**2 * r[nodes]**3 + 1e-14)
                else:
                    eps = 1e-7 * (1 + np.abs(Y[nodes, k]))
                Yp = Y.copy(); Yp[nodes, k] += eps
                d = self.residual(Yp, D, lam) - R0
                ef = np.zeros(N); ef[nodes] = eps
                j = np.where(idx % 2 == c, idx, idx + 1)
                for e in range(NV):
                    rr = nb0 + idx * NV + e
                    rows.append(rr); cols.append(j * NV + k); vals.append(d[rr] / ef[j])
                if (N - 1) % 2 == c:
                    rr = np.array([nrow - 2]); rows.append(rr); cols.append(np.array([(N - 1) * NV + k])); vals.append(d[rr] / eps[-1])
        # BC rows (analytic)
        rows.append(np.array([0, 1, 2, 3, 4])); cols.append(np.array([4, 1, 2, 3, 5])); vals.append(np.ones(5))
        # pin row
        rows.append(np.array([nrow - 1])); cols.append(np.array([self.kw * NV])); vals.append(np.array([2 * Y[self.kw, 0]]))
        # D column
        eD = 1e-6
        dD = (self.residual(Y, D + eD, lam) - R0) / eD
        nz = np.nonzero(dD)[0]
        rows.append(nz); cols.append(np.full(len(nz), N * NV)); vals.append(dD[nz])
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nrow, N * NV + 1))
        self._set(D)
        dl = np.zeros(nrow); dl[0] = -1.0
        return J, dl

    def safe(self, Y, D):
        r = self.r_of(D)
        return (np.all(np.isfinite(Y)) and np.all(Y[:, 2] / self.s0**2 < 0.5 * (1 - 1e-9) * r)
                and 1 - abs(D) * self.dphi_max > 0.2)

    def pack(self, Y, D): return np.concatenate([Y.ravel(), [D]])
    def unpack(self, U): return U[:-1].reshape(self.N, NV), U[-1]

    def solve(self, Y, D, lam, tol=1e-9, maxit=80):
        for it in range(maxit):
            R0 = self.residual(Y, D, lam)
            nr = np.abs(R0).max()
            if not np.isfinite(nr): return Y, D, False, it
            if nr < tol: return Y, D, True, it
            J, _ = self.jacobian(Y, D, lam, R0)
            dU = spla.spsolve(J, -R0)
            dY, dD = self.unpack(dU)
            n0 = np.linalg.norm(R0); t = 1.0
            while t > 1e-5:
                Yt, Dt = Y + t * dY, D + t * dD
                if self.safe(Yt, Dt):
                    n1 = np.linalg.norm(self.residual(Yt, Dt, lam))
                    if np.isfinite(n1) and n1 < n0 * (1 - 1e-4 * t): break
                t *= 0.5
            if t <= 1e-5: return Y, D, False, it
            Y, D = Yt, Dt
        return Y, D, False, maxit

    def observables(self, Y, D):
        self._set(D)
        o = self.base.observables(Y)
        o['D'] = D
        return o


def to_mm(model_from, Y, rfrom, s0, q, xw, hf, rmax, hc=0.1):
    """Interpolate a solution given on physical grid rfrom onto a new moving mesh centred at xw (D=0)."""
    xi, kw = make_xi(xw, hf, rmax, hc=hc)
    mm = MM(s0, q, xi, kw)
    Yn = np.empty((len(xi), NV))
    for k in range(NV):
        Yn[:, k] = np.interp(xi, rfrom, Y[:, k])
    return mm, Yn
