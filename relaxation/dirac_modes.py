"""
Radial pulsation modes of spherical Dirac stars, from the linearised equations of derivations/dirac_pulsation.py
(coefficients in dirac_pulsation_coeffs.py), polar-areal gauge, G = 1.

Perturbations ~ e^{lambda t} in the frame rotating with the star's frequency omega:
    lambda y_k = sum_j (A_kj + B_kj d/dr) y_j   for k = uR, uI, vR, vI   (Dirac),
    dX'   = sum_j (H_j + H1_j d/dr) y_j,   dphi' = sum_j (Q_j + Q1_j d/dr) y_j    (Hamiltonian constraint, slicing),
which is solved as the generalised eigenvalue problem  A Y = lambda B Y.

Discretisation: Chebyshev collocation on [-R, R] folded by parity (u, dX, dphi even, v odd, as in the cartoon evolution),
so only r > 0 points are used and the 1/r terms are regular. Optionally the Chebyshev coordinate is mapped to cluster the
points around a radius r_c (thin walls). Boundary conditions: the outermost Dirac rows are replaced by u = v = 0 (decay),
dX(0) = 0 (by interpolation to r = 0, replacing the dX equation at r = R), and dphi(R) + dX(R)/X(R) = 0 (alpha X -> 1,
replacing the dphi equation at r = R).

Spurious (unresolved) eigenvalues are removed by solving at two resolutions and keeping those that agree.
lambda = -i sigma: stable modes have real sigma (the field components oscillate at omega +- sigma). Re(lambda) > 0 is an
instability, oscillatory when Im(lambda) != 0. The global phase rotation gives lambda = 0.
"""
import sys, os
import numpy as np
from scipy.linalg import eig
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dirac_pulsation_coeffs as pc
from dirac_relax import Potential, DiracRelax


def cheb(N):
    """Chebyshev points x_j = cos(pi j / N), j = 0..N, and the differentiation matrix (Trefethen)"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    Xm = np.tile(x, (N + 1, 1)).T
    dX = Xm - Xm.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return x, D


class ModeSolver:
    def __init__(self, pot, sol, N=200, R=None, r_c=None, width=None, tortoise=0.):
        """pot: Potential; sol: the background (scipy solve_bvp solution of dirac_relax, y = f, g, m, phi, omega);
        N: number of collocation points with r > 0; r_c, width: optional clustering of points around r_c;
        tortoise in [0, 1): weight of the tortoise coordinate r* = int X/alpha dr in the collocation map, which resolves
        the short local wavelengths ~ alpha / (omega +- sigma) in the redshifted cores of ultracompact stars."""
        self.pot, self.sol = pot, sol
        R = R if R is not None else sol.x[-1]
        x, D = cheb(2 * N - 1)                 # 2N points on [-1, 1], none at 0
        if tortoise > 0.:
            # q(r) = (1 - t) r / R + t r*(r) / r*(R), an odd monotonic map; r(q) by interpolation, dr/dq from dq/dr
            rr = np.linspace(0., R, 200001)
            yb = sol.sol(np.clip(rr, sol.x[0], sol.x[-1]))
            Xb = 1 / np.sqrt(1 - 2 * np.where(rr > 0, yb[2] / np.maximum(rr, 1e-300), 0.))
            dstar = Xb * np.exp(-yb[3])
            rstar = np.concatenate([[0.], np.cumsum(0.5 * (dstar[1:] + dstar[:-1]) * np.diff(rr))])
            qq = (1 - tortoise) * rr / R + tortoise * rstar / rstar[-1]
            dqdr = (1 - tortoise) / R + tortoise * dstar / rstar[-1]
            ra = np.interp(np.abs(x), qq, rr) * np.sign(x)
            r_all, drdx = ra, 1. / np.interp(np.abs(x), qq, dqdr)
        elif r_c is None:
            r_all, drdx = R * x, R * np.ones_like(x)
        else:
            # r = R x + a [atan((x - xc)/w) + atan((x + xc)/w) - normalisation]: an odd map clustering points at +-r_c
            xc, wd = r_c / R, width / R
            s = lambda xx: np.arctan((xx - xc) / wd) + np.arctan((xx + xc) / wd)
            ds = lambda xx: (1 / wd) / (1 + ((xx - xc) / wd) ** 2) + (1 / wd) / (1 + ((xx + xc) / wd) ** 2)
            # x(q) with q uniform-Chebyshev: invert the map q = (s(x) / s(1)) mixed with the identity
            q = x
            mix = 0.5
            f_map = lambda xx: (1 - mix) * xx + mix * s(xx) / s(1.)
            xs = np.linspace(-1, 1, 200001); qs = f_map(xs)
            xx = np.interp(q, qs, xs)
            dqdx = (1 - mix) + mix * ds(xx) / s(1.)
            r_all, drdx = R * xx, R / dqdx       # dr/dq
        self.R = R
        pos = np.where(r_all > 0)[0]
        mirror = np.array([np.argmin(np.abs(r_all + r_all[j])) for j in pos])
        Dr = D / drdx[:, None]                   # d/dr = (dq/dr) d/dq, rows scaled
        self.r = r_all[pos]
        self.De = Dr[np.ix_(pos, pos)] + Dr[np.ix_(pos, mirror)]
        self.Do = Dr[np.ix_(pos, pos)] - Dr[np.ix_(pos, mirror)]
        # interpolation to r = 0 of an even function: barycentric weights of the 2N Chebyshev points at q = 0
        n = len(x)
        wts = (-1.) ** np.arange(n) * np.hstack([0.5, np.ones(n - 2), 0.5])
        q0 = 0.                                  # the map is odd, so r = 0 is q = 0
        lw = wts / (q0 - (x if r_c is None else x))
        lw = lw / lw.sum()
        self.interp0_even = lw[pos] + lw[mirror]
        self.N = len(pos)
        self._background()

    def _background(self):
        r = self.r
        s = DiracRelax(self.pot, R=self.sol.x[-1])
        yb = self.sol.sol(np.clip(r, self.sol.x[0], self.sol.x[-1]))
        out = r > self.sol.x[-1]
        # beyond the background domain: Schwarzschild exterior, field decaying
        if out.any():
            f, g, m, ph, w = self.sol.y[:, -1]
            k = np.sqrt(max(self.pot.mu**2 - (w / np.exp(ph))**2, 1e-12))
            yb[:, out] = np.vstack([f * np.exp(-k * (r[out] - self.sol.x[-1])), g * np.exp(-k * (r[out] - self.sol.x[-1])),
                                    np.full(out.sum(), m), 0.5 * np.log(1 - 2 * m / r[out]), np.full(out.sum(), w)])
        rhs = s.rhs(r, yb)
        f, g, m, ph, w = yb
        self.omega = w[0]
        X0 = 1 / np.sqrt(1 - 2 * m / r)
        S0 = (g * g - f * f) / (2 * np.pi)
        self.bg = dict(r=r, omega=self.omega, mu=self.pot.mu, g=g, f=f, dg=rhs[1], df=rhs[0], X0=X0, phi0=ph, dphi0=rhs[3],
                       U=self.pot.U(S0), dU=self.pot.dU(S0), ddU=self.pot.ddU(S0))

    def matrices(self):
        N, b = self.N, self.bg
        args = (b['r'], b['omega'], b['mu'], b['g'], b['f'], b['dg'], b['df'], b['X0'], b['phi0'], b['dphi0'], b['U'], b['dU'], b['ddU'])
        Dpar = [self.De, self.De, self.Do, self.Do, self.De, self.De]
        Id = np.eye(N)
        A = np.zeros((6 * N, 6 * N)); B = np.zeros((6 * N, 6 * N))
        def block(Acoef, Bcoef):
            return [np.diag(Acoef[j]) + np.diag(Bcoef[j]) @ Dpar[j] for j in range(6)]
        for k, fn in enumerate([pc.uR, pc.uI, pc.vR, pc.vI]):
            Ac, Bc = fn(*args)
            bl = block(Ac, Bc)
            for j in range(6):
                A[k * N:(k + 1) * N, j * N:(j + 1) * N] = bl[j]
            B[k * N:(k + 1) * N, k * N:(k + 1) * N] = Id
        for k, fn in ((4, pc.ham), (5, pc.sli)):
            Ac, Bc = fn(*args)
            bl = block(Ac, Bc)
            for j in range(6):
                A[k * N:(k + 1) * N, j * N:(j + 1) * N] = -bl[j]
            A[k * N:(k + 1) * N, k * N:(k + 1) * N] += Dpar[k]
        # boundary conditions at the outermost point (index of max r)
        jR = int(np.argmax(self.r))
        for k in range(4):
            row = k * N + jR
            A[row, :] = 0.; B[row, :] = 0.; A[row, row] = 1.
        row = 4 * N + jR                                    # dX(0) = 0
        A[row, :] = 0.; A[row, 4 * N:5 * N] = self.interp0_even
        row = 5 * N + jR                                    # dphi(R) + dX(R)/X(R) = 0
        A[row, :] = 0.; A[row, 5 * N + jR] = 1.; A[row, 4 * N + jR] = 1. / b['X0'][jR]
        return A, B

    def eigenvalues(self, return_vectors=False):
        A, B = self.matrices()
        lam, V = eig(A, B)
        ok = np.isfinite(lam)
        if return_vectors:
            return lam[ok], V[:, ok]
        return lam[ok]


def robust_modes(pot, sol, N1=160, N2=220, tol=1e-4, lam_max=5., **kw):
    """eigenvalues that agree (relative/absolute tol) between two resolutions, sorted by |lambda|"""
    l1 = ModeSolver(pot, sol, N=N1, **kw).eigenvalues()
    l2 = ModeSolver(pot, sol, N=N2, **kw).eigenvalues()
    l1 = l1[np.abs(l1) < lam_max]; l2 = l2[np.abs(l2) < lam_max]
    keep = []
    for z in l1:
        d = np.abs(l2 - z)
        if d.min() < tol * max(1., abs(z)):
            keep.append(0.5 * (z + l2[d.argmin()]))
    keep = np.array(keep)
    return keep[np.argsort(np.abs(keep))] if len(keep) else keep
