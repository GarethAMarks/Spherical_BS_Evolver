"""
Non-radial (polar, axisymmetric M = 0, multipole L >= 1) linear modes of free spherical Proca stars, from the projected
equations of derivations/proca_nonradial.py (coefficients in proca_nonradial_coeffs.py), Regge-Wheeler gauge,
signature (-,+,+,+), G = 1. Same formulation as dirac_nonradial_modes.py (see its docstring for the reasoning):

Fields (~ e^{lambda t} in the frame rotating with the star):
    matter  aR0 aI0 (dA_t = e^{-i w t} a0 Y), aR1 aI1 (dA_r = e^{-i w t} a1 Y), aR2 aI2 (dA_theta = e^{-i w t} a2 Y'),
    metric  H0, H1, H2, K  (K = 0 at L = 1, where the RW gauge leaves it free).
Equations (square system): the 6 projected Proca equations (t, r, theta components, real and imaginary parts) and the
Einstein components Eb, Erth, Etth, Etr (Eb dropped at L = 1, where it is not independent). The r and theta Proca
components contain lambda^2 (second time derivatives of a1, a2); they are linearised with auxiliary fields lambda * y.
Unused Einstein components (Ett, Err and the angular trace: Ea for L >= 2, Ea - Eb for L = 1, the only combination
present when Y'' = -Y) are evaluated on the eigenvectors as a check.
Centre: field = (r / sqrt(r^2 + r_s^2))^p x even function, p = L for a0, a2, H0, H2, K and L - 1 for a1, H1.
Outer boundary: matter fields vanish; K, H1, H0 obey F' + (lambda X / alpha + (L + 1) / r) F = 0.
Classification (classified()): agreement between two resolutions and residual of the unused components < 1e-2.
"""
import sys, os
import numpy as np
from scipy.linalg import eig
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proca_nonradial_coeffs as pc
from dirac_modes import cheb

MAT = ['aR0', 'aI0', 'aR1', 'aI1', 'aR2', 'aI2']
MET = ['H0', 'H1', 'H2', 'K']
PROCA_ROWS = ['Ptre', 'Ptim', 'Prre', 'Prim', 'Pthre', 'Pthim']


def background_on(sol, mu, r):
    """background (f, E, m, phi, omega) on r, with a Schwarzschild exterior and decaying field beyond the domain"""
    yb = sol.sol(np.clip(r, sol.x[0], sol.x[-1]))
    out = r > sol.x[-1]
    if out.any():
        f, E, m, ph, w = sol.y[:, -1]
        k = np.sqrt(max(mu**2 - (w / np.exp(ph))**2, 1e-12))
        dec = np.exp(-k * (r[out] - sol.x[-1]))
        yb[:, out] = np.vstack([f * dec, E * dec, np.full(out.sum(), m), 0.5 * np.log(1 - 2 * m / r[out]), np.full(out.sum(), w)])
    return yb


class ProcaModes:
    def __init__(self, sol, L, N=120, R=None, mu=1., r_s=2., fields=None, rows=None, bc=None, regular=True):
        self.sol, self.L, self.mu = sol, L, mu
        R = R if R is not None else sol.x[-1]
        self.R = R
        x, D = cheb(2 * N - 1)
        r_all = R * x
        pos = np.where(r_all > 0)[0]
        mirror = np.array([np.argmin(np.abs(r_all + r_all[j])) for j in pos])
        Dr = D / R; D2 = Dr @ Dr
        self.r = r_all[pos]; self.N = len(pos)
        self.D = {+1: Dr[np.ix_(pos, pos)] + Dr[np.ix_(pos, mirror)], -1: Dr[np.ix_(pos, pos)] - Dr[np.ix_(pos, mirror)]}
        self.D2 = {+1: D2[np.ix_(pos, pos)] + D2[np.ix_(pos, mirror)], -1: D2[np.ix_(pos, pos)] - D2[np.ix_(pos, mirror)]}
        if fields is None:
            fields = MAT + (['H0', 'H1', 'H2'] if L == 1 else MET)
        if rows is None:
            rows = PROCA_ROWS + (['Erth', 'Etth', 'Etr'] if L == 1 else ['Eb', 'Erth', 'Etth', 'Etr'])
        self.fields, self.rowlist = list(fields), list(rows)
        self.idx = {n: k for k, n in enumerate(self.fields)}
        self.bcmap = bc if bc is not None else dict(Ptre='aR0', Ptim='aI0', Prre='aR1', Prim='aI1', Pthre='aR2', Pthim='aI2',
                                                     Erth='H0', Etth='H1', Etr=('H2' if L == 1 else 'K'))
        # leading powers at the centre
        self.pw = {n: (L - 1 if n in ('aR1', 'aI1', 'H1') else L) for n in MAT + MET}
        if L == 0:
            self.pw.update(aR1=1, aI1=1)
        self.regular = regular
        rr = self.r; q = rr / np.sqrt(rr**2 + r_s**2); dq = r_s**2 / (rr**2 + r_s**2)**1.5
        ddq = -3 * rr * r_s**2 / (rr**2 + r_s**2)**2.5
        self.wts = {}
        for n, p in self.pw.items():
            w = q**p
            dw = p * q**(p - 1) * dq if p > 0 else 0 * rr
            ddw = (p * (p - 1) * q**(p - 2) * dq**2 if p > 1 else 0 * rr) + (p * q**(p - 1) * ddq if p > 0 else 0 * rr)
            self.wts[n] = (w, dw, ddw)
        f, E, m, ph, w = background_on(sol, mu, self.r)
        self.omega = w[0]
        self.bgargs = (self.r, self.omega, mu, float(L * (L + 1)), f, E, 1 / np.sqrt(1 - 2 * m / self.r), np.exp(ph))

    def _ops(self, n):
        w, dw, ddw = self.wts[n]
        I0 = np.diag(w); I1 = np.diag(w) @ self.D[1] + np.diag(dw)
        I2 = np.diag(w) @ self.D2[1] + 2 * np.diag(dw) @ self.D[1] + np.diag(ddw)
        return I0, I1, I2

    def _op(self, coeffs):
        N, nf = self.N, len(self.fields)
        A = np.zeros((N, nf * N)); B = np.zeros((N, nf * N)); C = np.zeros((N, nf * N))
        for (n, d), c in coeffs.items():
            if n not in self.idx:
                continue
            j = self.idx[n]; sl = slice(j * N, (j + 1) * N)
            I0, I1, I2 = self._ops(n)
            c = np.diag(c)
            if d == '': A[:, sl] += c @ I0
            elif d == '_r': A[:, sl] += c @ I1
            elif d == '_rr': A[:, sl] += c @ I2
            elif d == '_t': B[:, sl] += c @ I0
            elif d == '_tr': B[:, sl] += c @ I1
            elif d == '_tt': C[:, sl] += c @ I0
            else: raise ValueError(d)
        return A, B, C

    def matrices(self):
        """linear pencil (A, B), A z = lambda B z, z = (fields, lambda * fields with second time derivatives)"""
        N, nf = self.N, len(self.fields)
        jR = int(np.argmax(self.r))
        A = np.zeros((nf * N, nf * N)); B = np.zeros_like(A); C = np.zeros_like(A)
        for k, name in enumerate(self.rowlist):
            Ar, Br, Cr = self._op(pc.ROWS[name](*self.bgargs))
            sl = slice(k * N, (k + 1) * N)
            A[sl] = Ar; B[sl] = Br; C[sl] = Cr
            fb = self.bcmap.get(name)
            if fb is not None:
                row = k * N + jR
                A[row] = 0.; B[row] = 0.; C[row] = 0.
                j = self.idx[fb]; s2 = slice(j * N, (j + 1) * N)
                if fb in MAT:
                    A[row, j * N + jR] = 1.
                else:
                    I0, I1, _ = self._ops(fb)
                    Xb, alb = self.bgargs[6][jR], self.bgargs[7][jR]
                    A[row, s2] = I1[jR] + (self.L + 1) / self.r[jR] * I0[jR]
                    B[row, s2] = (Xb / alb) * I0[jR]
        # companion linearisation of lambda^2 C for the field blocks that need it
        aux = [j for j in range(nf) if C[:, j * N:(j + 1) * N].any()]
        self.aux = aux
        na = len(aux)
        Z = nf * N + na * N
        AA = np.zeros((Z, Z)); BB = np.zeros((Z, Z))
        AA[:nf * N, :nf * N] = A; BB[:nf * N, :nf * N] = -B
        for q, j in enumerate(aux):
            cols = slice(nf * N + q * N, nf * N + (q + 1) * N)
            BB[:nf * N, cols] = -C[:, j * N:(j + 1) * N]
            AA[cols, cols] = np.eye(N)
            BB[cols, j * N:(j + 1) * N] = np.eye(N)
        return AA, BB

    def eigen(self, vectors=False):
        A, B = self.matrices()
        lam, V = eig(A, B)
        ok = np.isfinite(lam) & (np.abs(lam) < 1e8)
        return (lam[ok], V[:, ok]) if vectors else lam[ok]

    def split(self, v):
        N = self.N
        return {n: v[k * N:(k + 1) * N] for n, k in self.idx.items()}

    def unused(self):
        if self.L == 1:
            return ['Ett', 'Err', 'Ea-Eb']
        return [n for n in ('Ett', 'Err', 'Ea', 'Etr', 'Eb', 'Erth', 'Etth') if n not in self.rowlist]

    def residuals(self, lam, v, names=None):
        v = v[:len(self.fields) * self.N]
        out = {}
        for name in (names or self.unused()):
            parts = [(1., 'Ea'), (-1., 'Eb')] if name == 'Ea-Eb' else [(1., name)]
            res = 0.; scale = 0.
            for sgn, nm in parts:
                Ar, Br, Cr = self._op(pc.ROWS[nm](*self.bgargs))
                M = Ar + lam * Br + lam**2 * Cr
                res = res + sgn * (M @ v)
                scale = scale + (np.abs(Ar) + abs(lam) * np.abs(Br) + abs(lam)**2 * np.abs(Cr)) @ np.abs(v)
            inner = self.r < 0.8 * self.R
            out[name] = np.max(np.abs(res[inner])) / max(np.max(scale[inner]), 1e-300)
        return out


def classified(sol, L, N1=80, N2=110, R=None, lam_max=2.5, tol=2e-3, atol=3e-4, **kw):
    ms = [ProcaModes(sol, L, N=N, R=R, **kw) for N in (N1, N2)]
    l1 = ms[0].eigen(); l1 = l1[(np.abs(l1) < lam_max) & (np.abs(l1) > 1e-8)]
    l2, V2 = ms[1].eigen(vectors=True)
    s = (np.abs(l2) < lam_max) & (np.abs(l2) > 1e-8); l2 = l2[s]; V2 = V2[:, s]
    out = []
    for k, z in enumerate(l2):
        conv = len(l1) > 0 and np.min(np.abs(l1 - z)) < tol * abs(z) + atol
        res = ms[1].residuals(z, V2[:, k])
        out.append((z, conv, max(res.values())))
    out.sort(key=lambda t: -t[0].real)
    return out
