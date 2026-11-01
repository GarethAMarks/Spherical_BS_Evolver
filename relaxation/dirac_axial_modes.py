"""
Axial-parity (parity (-1)^(L+1)), M = 0 linear modes of free spherical Dirac stars, from derivations/dirac_axial.py
(coefficients in dirac_axial_coeffs.py), Regge-Wheeler gauge, signature (-,+,+,+), G = 1. Together with the polar sector
(dirac_nonradial_modes.py) this covers all non-radial modes: by the combined rotation x flavour symmetry, M != 0 modes are
rotations of M = 0 modes with the same L and parity.

Fields (~ e^{lambda t} in the rotating frame): cR1 cI1 cR2 cI2 (upper components, (sigma.n)(Y chi_up, e^{i phi} Y' chi_dn)),
dR1 dI1 dR2 dI2 (lower components), h0 (h_tphi = h0 sin Y'), h1 (h_rphi = h1 sin Y'; h1 = 0 at L = 1, a gauge choice).
Equations: the 8 projected Dirac equations and Etph, Erph (L >= 2) or Etph (L = 1). By the phi component of the Bianchi
identity, Ethph = 0 then follows for L >= 2, and at L = 1 Erph obeys a first-order equation whose solution is singular at
the centre; the unused component is evaluated on every eigenvector.
Centre: field = (r / sqrt(r^2 + r_s^2))^p x even function, p = L - 1 (c), L (d), L + 1 (h0, h1).
Outer boundary: matter vanishes; h0, h1 obey F' + (lambda X / alpha + L / r) F = 0 (static exterior h0 ~ r^-L).
"""
import sys, os
import numpy as np
from scipy.linalg import eig
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dirac_axial_coeffs as ac
from dirac_modes import cheb
from dirac_nonradial_modes import background_on

MAT = ['cR1', 'cI1', 'cR2', 'cI2', 'dR1', 'dI1', 'dR2', 'dI2']
DROWS = ['A1re', 'A1im', 'A2re', 'A2im', 'C1re', 'C1im', 'C2re', 'C2im']


class AxialModes:
    def __init__(self, sol, L, N=120, R=None, mu=1., r_s=2., cowling=False):
        self.sol, self.L, self.mu = sol, L, mu
        R = R if R is not None else sol.x[-1]
        self.R = R
        x, D = cheb(2 * N - 1)
        r_all = R * x
        pos = np.where(r_all > 0)[0]
        mirror = np.array([np.argmin(np.abs(r_all + r_all[j])) for j in pos])
        Dr = D / R; D2 = Dr @ Dr
        self.r = r_all[pos]; self.N = len(pos)
        self.D = Dr[np.ix_(pos, pos)] + Dr[np.ix_(pos, mirror)]          # even functions
        self.D2 = D2[np.ix_(pos, pos)] + D2[np.ix_(pos, mirror)]
        if cowling:
            self.fields, self.rowlist = list(MAT), list(DROWS)
        elif L == 1:
            self.fields, self.rowlist = MAT + ['h0'], DROWS + ['Etph']
        else:
            self.fields, self.rowlist = MAT + ['h0', 'h1'], DROWS + ['Etph', 'Erph']
        self.idx = {n: k for k, n in enumerate(self.fields)}
        self.bcmap = {n: f for n, f in zip(DROWS, MAT)}
        self.bcmap.update(Etph='h0', Erph='h1')
        self.unused = ['Erph'] if L == 1 else ['Ethph']
        rr = self.r; q = rr / np.sqrt(rr**2 + r_s**2); dq = r_s**2 / (rr**2 + r_s**2)**1.5
        ddq = -3 * rr * r_s**2 / (rr**2 + r_s**2)**2.5
        self.wts = {}
        for n in MAT + ['h0', 'h1']:
            p = L - 1 if n[0] == 'c' else (L if n[0] == 'd' else L + 1)
            w = q**p
            dw = p * q**(p - 1) * dq if p > 0 else 0 * rr
            ddw = (p * (p - 1) * q**(p - 2) * dq**2 if p > 1 else 0 * rr) + (p * q**(p - 1) * ddq if p > 0 else 0 * rr)
            self.wts[n] = (w, dw, ddw)
        f, g, m, ph, w = background_on(sol, mu, self.r)
        self.omega = w[0]
        self.bgargs = (self.r, self.omega, mu, float(L * (L + 1)), g, f, 1 / np.sqrt(1 - 2 * m / self.r), np.exp(ph))

    def _ops(self, n):
        w, dw, ddw = self.wts[n]
        return (np.diag(w), np.diag(w) @ self.D + np.diag(dw),
                np.diag(w) @ self.D2 + 2 * np.diag(dw) @ self.D + np.diag(ddw))

    def _op(self, coeffs):
        N, nf = self.N, len(self.fields)
        A = np.zeros((N, nf * N)); B = np.zeros((N, nf * N)); C = np.zeros((N, nf * N))
        for (n, d), c in coeffs.items():
            if n not in self.idx:
                continue
            j = self.idx[n]; sl = slice(j * N, (j + 1) * N)
            I0, I1, I2 = self._ops(n); c = np.diag(c)
            if d == '': A[:, sl] += c @ I0
            elif d == '_r': A[:, sl] += c @ I1
            elif d == '_rr': A[:, sl] += c @ I2
            elif d == '_t': B[:, sl] += c @ I0
            elif d == '_tr': B[:, sl] += c @ I1
            elif d == '_tt': C[:, sl] += c @ I0
            else: raise ValueError(d)
        return A, B, C

    def matrices(self):
        N, nf = self.N, len(self.fields)
        jR = int(np.argmax(self.r))
        A = np.zeros((nf * N, nf * N)); B = np.zeros_like(A); C = np.zeros_like(A)
        for k, name in enumerate(self.rowlist):
            Ar, Br, Cr = self._op(ac.ROWS[name](*self.bgargs))
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
                    A[row, s2] = I1[jR] + self.L / self.r[jR] * I0[jR]
                    B[row, s2] = (Xb / alb) * I0[jR]
        aux = [j for j in range(nf) if C[:, j * N:(j + 1) * N].any()]
        na = len(aux); Z = nf * N + na * N
        AA = np.zeros((Z, Z)); BB = np.zeros((Z, Z))
        AA[:nf * N, :nf * N] = A; BB[:nf * N, :nf * N] = -B
        for q, j in enumerate(aux):
            cols = slice(nf * N + q * N, nf * N + (q + 1) * N)
            BB[:nf * N, cols] = -C[:, j * N:(j + 1) * N]
            AA[cols, cols] = np.eye(N); BB[cols, j * N:(j + 1) * N] = np.eye(N)
        return AA, BB

    def eigen(self, vectors=False):
        A, B = self.matrices()
        lam, V = eig(A, B)
        ok = np.isfinite(lam) & (np.abs(lam) < 1e8)
        return (lam[ok], V[:, ok]) if vectors else lam[ok]

    def split(self, v):
        return {n: v[k * self.N:(k + 1) * self.N] for n, k in self.idx.items()}

    def residuals(self, lam, v, names=None):
        v = v[:len(self.fields) * self.N]
        out = {}
        for name in (names or self.unused):
            Ar, Br, Cr = self._op(ac.ROWS[name](*self.bgargs))
            res = (Ar + lam * Br + lam**2 * Cr) @ v
            sc = (np.abs(Ar) + abs(lam) * np.abs(Br) + abs(lam)**2 * np.abs(Cr)) @ np.abs(v)
            inner = self.r < 0.8 * self.R
            out[name] = np.max(np.abs(res[inner])) / max(np.max(sc[inner]), 1e-300)
        return out


def classified(sol, L, N1=80, N2=110, R=None, lam_max=2.5, tol=2e-3, atol=3e-4, **kw):
    ms = [AxialModes(sol, L, N=N, R=R, **kw) for N in (N1, N2)]
    l1 = ms[0].eigen(); l1 = l1[(np.abs(l1) < lam_max) & (np.abs(l1) > 1e-8)]
    l2, V2 = ms[1].eigen(vectors=True)
    s = (np.abs(l2) < lam_max) & (np.abs(l2) > 1e-8); l2 = l2[s]; V2 = V2[:, s]
    out = []
    for k, z in enumerate(l2):
        conv = len(l1) > 0 and np.min(np.abs(l1 - z)) < tol * abs(z) + atol
        res = ms[1].residuals(z, V2[:, k]) if len(ms[1].fields) > 8 else {'none': 0.}
        out.append((z, conv, max(res.values())))
    out.sort(key=lambda t: -t[0].real)
    return out
