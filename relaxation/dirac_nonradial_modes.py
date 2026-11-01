"""
Non-radial (polar, axisymmetric M = 0, total angular momentum L >= 1) linear modes of free spherical Dirac stars, from the
projected equations of derivations/dirac_nonradial.py (coefficients in dirac_nonradial_coeffs.py), Regge-Wheeler gauge,
signature (-,+,+,+), G = 1.

Fields (all ~ e^{lambda t} in the frame rotating with the star):
    matter  aR1 aI1 aR2 aI2 (upper components: Y chi_up and spin-flipped e^{i phi} Y' chi_dn of the + spinor),
            bR1 bI1 bR2 bI2 (lower components, sigma.n times the same),
    metric  H0, H1, H2, K  (h_tt = alpha^2 H0 Y, h_tr = H1 Y, h_rr = X^2 H2 Y, h_AB = r^2 K Y Omega_AB).
Equations (default, square system): the 8 projected Dirac equations and the Einstein components Eb (traceless angular
part, algebraic in H0 - H2; note it is sourced by the spin-flip amplitudes aR2, bR2, so H0 != H2 in general), Erth (r-theta),
Etth (t-theta momentum constraint) and Etr (r momentum constraint). For lambda != 0 the projected Bianchi identities
then imply Ett = Ea = 0, and Err obeys a first-order equation whose solution is singular at the centre, so it vanishes for
regular modes. The unused components (Ett, Err, Ea) are evaluated on every eigenvector (residuals()) as a check.
At lambda = 0 exactly the square system admits static Hamiltonian-violating solutions; these are discarded (|lambda| <
1e-8), so zero modes (bifurcations) are located by following eigenvalues through zero along a family.
For L = 1 the RW gauge leaves residual gauge freedom (K = 0 is imposed); the L = 1 spectrum is not reliable with the
square system (rect() / eigen_projected() with full=False recovers the translation zero mode).

Discretisation: Chebyshev collocation on [-R, R] folded by parity, with the leading power at the centre factored out
(field = w(r) x even function, w = (r / sqrt(r^2 + r_s^2))^p, p = L for a, H0, H2, K and L - 1 for b, H1), so that
singular branches cannot be represented. Outer boundary: matter fields vanish at R; the metric functions K, H1, H0 obey the
decaying/outgoing Robin condition F' + (lambda X / alpha + (L + 1) / r) F = 0 (linear in lambda). Dirichlet conditions on
the metric are NOT used: they are inconsistent with the slowly decaying exterior metric, and the boundary residual spreads
through the spectral derivative into O(1) constraint violations in the interior.
lambda = -i sigma: Re lambda > 0 is an instability. Radiating modes acquire Re lambda < 0 through the Robin condition.

Classification (classified()): an eigenvalue is physical if it agrees between two resolutions AND the un-imposed Einstein
components vanish on its eigenvector (relative residual < 1e-2). Spurious growing modes (a high-frequency band near
|Im lambda| ~ mu + omega, and weak continuum box modes with Re lambda < 1e-3 whose frequencies move with R) fail one or both.
"""
import sys, os
import numpy as np
from scipy.linalg import eig
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dirac_nonradial_coeffs as nc
from dirac_modes import cheb

MAT = ['aR1', 'aI1', 'aR2', 'aI2', 'bR1', 'bI1', 'bR2', 'bI2']
MET = ['H0', 'H1', 'H2', 'K']
DIRAC_ROWS = ['A1re', 'A1im', 'A2re', 'A2im', 'C1re', 'C1im', 'C2re', 'C2im']


def background_on(sol, pot_mu, r):
    """background (f, g, m, phi, omega) on r, with a Schwarzschild exterior and decaying field beyond the solution domain"""
    yb = sol.sol(np.clip(r, sol.x[0], sol.x[-1]))
    out = r > sol.x[-1]
    if out.any():
        f, g, m, ph, w = sol.y[:, -1]
        k = np.sqrt(max(pot_mu**2 - (w / np.exp(ph))**2, 1e-12))
        dec = np.exp(-k * (r[out] - sol.x[-1]))
        yb[:, out] = np.vstack([f * dec, g * dec, np.full(out.sum(), m), 0.5 * np.log(1 - 2 * m / r[out]), np.full(out.sum(), w)])
    return yb


class NonRadialModes:
    def __init__(self, sol, L, N=120, R=None, mu=1., tortoise=0., cowling=False, einstein=None, bc=None, regular=True,
                 r_s=2., outer='robin'):
        """sol: dirac_relax solve_bvp solution (y = f, g, m, phi, omega); L >= 1; N collocation points with r > 0;
        cowling: freeze the metric (H = 0), for comparison"""
        self.sol, self.L, self.mu, self.cowling = sol, L, mu, cowling
        self.einstein = einstein if einstein is not None else (['Erth', 'Etth', 'Etr'] if L == 1 else ['Eb', 'Erth', 'Etth', 'Etr'])
        self.bc_metric = bc if bc is not None else ({'Erth': 'H0', 'Etth': 'H1', 'Etr': 'H2'} if L == 1 else
                                                   {'Erth': 'H1', 'Etth': 'H2', 'Etr': 'K'})
        R = R if R is not None else sol.x[-1]
        self.R = R
        x, D = cheb(2 * N - 1)
        if tortoise > 0.:
            rr = np.linspace(0., R, 200001)
            yb = background_on(sol, mu, np.maximum(rr, 1e-8))
            Xb = 1 / np.sqrt(1 - 2 * yb[2] / np.maximum(rr, 1e-8))
            dstar = Xb * np.exp(-yb[3])
            rstar = np.concatenate([[0.], np.cumsum(0.5 * (dstar[1:] + dstar[:-1]) * np.diff(rr))])
            qq = (1 - tortoise) * rr / R + tortoise * rstar / rstar[-1]
            dqdr = (1 - tortoise) / R + tortoise * dstar / rstar[-1]
            r_all = np.interp(np.abs(x), qq, rr) * np.sign(x)
            drdx = 1. / np.interp(np.abs(x), qq, dqdr)
        else:
            r_all, drdx = R * x, R * np.ones_like(x)
        pos = np.where(r_all > 0)[0]
        mirror = np.array([np.argmin(np.abs(r_all + r_all[j])) for j in pos])
        Dr = D / drdx[:, None]
        D2 = Dr @ Dr
        self.r = r_all[pos]
        self.N = len(pos)
        self.D = {+1: Dr[np.ix_(pos, pos)] + Dr[np.ix_(pos, mirror)], -1: Dr[np.ix_(pos, pos)] - Dr[np.ix_(pos, mirror)]}
        self.D2 = {+1: D2[np.ix_(pos, pos)] + D2[np.ix_(pos, mirror)], -1: D2[np.ix_(pos, pos)] - D2[np.ix_(pos, mirror)]}
        pL = (-1) ** L
        self.parity = {n: (pL if n[0] == 'a' else -pL) for n in MAT}
        self.parity.update(H0=pL, H2=pL, K=pL, H1=-pL)
        # regular factorisation: field = w(r) * (even smooth function), w = (r / sqrt(r^2 + r_s^2))^p, with p the
        # leading power at the centre: L for a, H0, H2, K; L - 1 for b and H1. Singular (e.g. constraint-violating)
        # branches then cannot be represented.
        self.regular = regular
        self.outer = outer
        if regular:
            rr = self.r; q = rr / np.sqrt(rr**2 + r_s**2); dq = r_s**2 / (rr**2 + r_s**2)**1.5
            ddq = -3 * rr * r_s**2 / (rr**2 + r_s**2)**2.5
            self.wts = {}
            for n in MAT + MET:
                p = L if (n[0] == 'a' or n in ('H0', 'H2', 'K')) else L - 1
                w = q**p
                dw = p * q**(p - 1) * dq if p > 0 else 0 * rr
                ddw = (p * (p - 1) * q**(p - 2) * dq**2 if p > 1 else 0 * rr) + (p * q**(p - 1) * ddq if p > 0 else 0 * rr)
                self.wts[n] = (w, dw, ddw)
        if cowling:
            self.fields = list(MAT)
        else:
            self.fields = MAT + (['H0', 'H1', 'H2'] if L == 1 else MET)
        self.idx = {n: k for k, n in enumerate(self.fields)}
        yb = background_on(sol, mu, self.r)
        f, g, m, ph, w = yb
        self.omega = w[0]
        self.bgargs = (self.r, self.omega, mu, float(L * (L + 1)), g, f, 1 / np.sqrt(1 - 2 * m / self.r), np.exp(ph))

    def rows(self):
        if self.cowling:
            return DIRAC_ROWS
        return DIRAC_ROWS + list(self.einstein)

    def _op(self, coeffs):
        """(A_row, B_row, C_row): row operator = A + lambda B + lambda^2 C acting on the stacked fields"""
        N, nf = self.N, len(self.fields)
        A = np.zeros((N, nf * N)); B = np.zeros((N, nf * N)); C = np.zeros((N, nf * N))
        for (n, d), c in coeffs.items():
            if n not in self.idx:
                continue          # field frozen (K for L = 1, metric in the Cowling approximation)
            j = self.idx[n]; sl = slice(j * N, (j + 1) * N)
            if self.regular:
                w, dw, ddw = self.wts[n]
                I0 = np.diag(w); I1 = np.diag(w) @ self.D[1] + np.diag(dw)
                I2 = np.diag(w) @ self.D2[1] + 2 * np.diag(dw) @ self.D[1] + np.diag(ddw)
            else:
                p = self.parity[n]
                I0 = np.eye(N); I1 = self.D[p]; I2 = self.D2[p]
            c = np.diag(c)
            if d == '':
                A[:, sl] += c @ I0
            elif d == '_r':
                A[:, sl] += c @ I1
            elif d == '_rr':
                A[:, sl] += c @ I2
            elif d == '_t':
                B[:, sl] += c @ I0
            elif d == '_tr':
                B[:, sl] += c @ I1
            elif d == '_tt':
                C[:, sl] += c @ I0
            else:
                raise ValueError(d)
        return A, B, C

    def matrices(self):
        rows = self.rows()
        N, nf = self.N, len(self.fields)
        A = np.zeros((nf * N, nf * N)); B = np.zeros_like(A)
        jR = int(np.argmax(self.r))
        for k, name in enumerate(rows):
            Ar, Br, Cr = self._op(nc.ROWS[name](*self.bgargs))
            assert not Cr.any(), name
            A[k * N:(k + 1) * N] = Ar; B[k * N:(k + 1) * N] = -Br
        self.nrows = len(rows)
        # outer boundary: Dirichlet on the matter fields (replacing the Dirac rows) and on H1, K, H2 (replacing Erth,
        # Ett, Etth respectively); for L = 1 H2 replaces Ett
        bc = {name: fld for name, fld in zip(DIRAC_ROWS, MAT)}
        if not self.cowling:
            bc.update(self.bc_metric)
        for k, name in enumerate(rows):
            if name in bc:
                row = k * N + jR
                A[row] = 0.; B[row] = 0.
                fld = bc[name]
                j = self.idx[fld]; sl = slice(j * N, (j + 1) * N)
                if fld in MAT or self.outer == 'dirichlet':
                    A[row, j * N + jR] = 1.
                else:
                    # outgoing / decaying exterior: F' + (lambda X / alpha + (L + 1) / r) F = 0 at R
                    I0, I1 = self._ops1(fld)
                    Xb, alb = self.bgargs[6][jR], self.bgargs[7][jR]
                    A[row, sl] = I1[jR] + (self.L + 1) / self.r[jR] * I0[jR]
                    B[row, sl] = -(Xb / alb) * I0[jR]
        return A, B

    def _ops1(self, n):
        N = self.N
        if self.regular:
            w, dw, ddw = self.wts[n]
            return np.diag(w), np.diag(w) @ self.D[1] + np.diag(dw)
        p = self.parity[n]
        return np.eye(N), self.D[p]

    def rect(self):
        """rectangular pencil (A, B) with all first-order Einstein components Ett, Etr, Etth, Erth, Eb (13 row blocks for
        12 fields; 12 for 11 at L = 1), rows normalised; outer boundary rows replaced by Dirichlet conditions"""
        names = DIRAC_ROWS + (['Ett', 'Etr', 'Etth', 'Erth'] if self.L == 1 else ['Ett', 'Etr', 'Etth', 'Erth', 'Eb'])
        bcm = {'Ett': 'K', 'Etr': 'H2', 'Etth': 'H1', 'Erth': 'H0'}
        if self.L == 1:
            bcm = {'Ett': 'H2', 'Etr': None, 'Etth': 'H1', 'Erth': 'H0'}
        bcm.update({n: f for n, f in zip(DIRAC_ROWS, MAT)})
        N, nf = self.N, len(self.fields)
        jR = int(np.argmax(self.r))
        As, Bs = [], []
        for name in names:
            Ar, Br, Cr = self._op(nc.ROWS[name](*self.bgargs))
            assert not Cr.any()
            Br = -Br
            fb = bcm.get(name)
            if fb is not None:
                Ar[jR] = 0.; Br[jR] = 0.; Ar[jR, self.idx[fb] * N + jR] = 1.
            As.append(Ar); Bs.append(Br)
        A = np.vstack(As); B = np.vstack(Bs)
        nrm = np.sqrt((A**2).sum(1) + (B**2).sum(1))
        return A / nrm[:, None], B / nrm[:, None]

    def rect_full(self):
        """rectangular pencil with ALL seven Einstein components (L >= 2): the lambda^2 terms (K_tt, H2_tt in Err, Ea) are
        linearised with auxiliary fields Kd = lambda K, H2d = lambda H2.  Fields: 12 + 2 blocks; rows: 8 Dirac + 7
        Einstein + 2 auxiliary.  The three extra rows are the Bianchi identities, so true modes are exact null vectors
        (up to discretisation error)."""
        assert self.L >= 2 and not self.cowling
        N = self.N; nf = len(self.fields) + 2
        jR = int(np.argmax(self.r))
        iKd, iH2d = len(self.fields), len(self.fields) + 1
        names = DIRAC_ROWS + ['Ett', 'Etr', 'Err', 'Etth', 'Erth', 'Ea', 'Eb']
        bcm = {n: f for n, f in zip(DIRAC_ROWS, MAT)}
        bcm.update({'Ett': 'K', 'Etr': 'H2', 'Etth': 'H1', 'Erth': 'H0'})
        As, Bs = [], []
        for name in names:
            Ar, Br, Cr = self._op(nc.ROWS[name](*self.bgargs))
            A = np.zeros((N, nf * N)); B = np.zeros((N, nf * N))
            A[:, :Ar.shape[1]] = Ar; B[:, :Br.shape[1]] = -Br
            # lambda^2 C y: C acts only on K and H2 columns -> lambda (C on Kd, H2d)
            for fld, iaux in (('K', iKd), ('H2', iH2d)):
                j = self.idx[fld]
                blk = Cr[:, j * N:(j + 1) * N]
                if blk.any():
                    B[:, iaux * N:(iaux + 1) * N] -= blk
            other = Cr.copy()
            for fld in ('K', 'H2'):
                j = self.idx[fld]; other[:, j * N:(j + 1) * N] = 0.
            assert not other.any(), name
            fb = bcm.get(name)
            if fb is not None:
                A[jR] = 0.; B[jR] = 0.; A[jR, self.idx[fb] * N + jR] = 1.
            As.append(A); Bs.append(B)
        Id = np.eye(N)
        for fld, iaux in (('K', iKd), ('H2', iH2d)):
            A = np.zeros((N, nf * N)); B = np.zeros((N, nf * N))
            A[:, iaux * N:(iaux + 1) * N] = Id; B[:, self.idx[fld] * N:(self.idx[fld] + 1) * N] = Id
            As.append(A); Bs.append(B)
        A = np.vstack(As); B = np.vstack(Bs)
        nrm = np.sqrt((A**2).sum(1) + (B**2).sum(1))
        return A / nrm[:, None], B / nrm[:, None]

    def eigen_projected(self, seeds=(1, 2), tol=1e-6, vectors=False, full=True):
        """eigenvalues of the rectangular pencil: those of W^T (A - lambda B) common to two random projections W"""
        A, B = self.rect_full() if (full and self.L >= 2) else self.rect()
        res = []
        for sd in seeds:
            rng = np.random.default_rng(sd)
            W = rng.standard_normal((A.shape[0], A.shape[1]))
            l, V = eig(W.T @ A, W.T @ B)
            ok = np.isfinite(l) & (np.abs(l) < 1e8)
            res.append((l[ok], V[:, ok]))
        (l1, V1), (l2, _) = res
        keep = [k for k, z in enumerate(l1) if np.min(np.abs(l2 - z)) < tol * max(1., abs(z))]
        if vectors:
            return l1[keep], V1[:, keep]
        return l1[keep]

    def smin(self, lam, full=True):
        """smallest singular value of the rectangular A - lambda B (relative)"""
        A, B = self.rect_full() if (full and self.L >= 2) else self.rect()
        M = A - lam * B
        sv = np.linalg.svd(M, compute_uv=False)
        return sv[-1] / sv[0]

    def eigen(self, vectors=False):
        A, B = self.matrices()
        lam, V = eig(A, B)
        ok = np.isfinite(lam) & (np.abs(lam) < 1e8)
        return (lam[ok], V[:, ok]) if vectors else lam[ok]

    def split(self, v):
        N = self.N
        return {n: v[k * N:(k + 1) * N] for n, k in self.idx.items()}

    def residuals(self, lam, v):
        """relative residuals of the Einstein components not used (Etr, Err, Ea) on an eigenvector"""
        out = {}
        v = v[:len(self.fields) * self.N]
        for name in ('Ett', 'Etr', 'Err', 'Ea', 'Eb', 'Erth', 'Etth'):
            Ar, Br, Cr = self._op(nc.ROWS[name](*self.bgargs))
            res = (Ar + lam * Br + lam**2 * Cr) @ v
            scale = (np.abs(Ar) + abs(lam) * np.abs(Br) + abs(lam)**2 * np.abs(Cr)) @ np.abs(v)
            inner = self.r < 0.8 * self.R
            out[name] = np.max(np.abs(res[inner])) / max(np.max(scale[inner]), 1e-300)
        return out


def robust(sol, L, N1=100, N2=140, R=None, tol=2e-3, lam_max=3., **kw):
    """eigenvalues agreeing between two resolutions, sorted by decreasing Re"""
    ls = []
    for N in (N1, N2):
        l = NonRadialModes(sol, L, N=N, R=R, **kw).eigen()
        ls.append(l[np.abs(l) < lam_max])
    keep = []
    for z in ls[0]:
        d = np.abs(ls[1] - z)
        if len(d) and d.min() < tol * max(1., abs(z)):
            keep.append(0.5 * (z + ls[1][d.argmin()]))
    keep = np.array(keep)
    return keep[np.argsort(-keep.real)] if len(keep) else keep


def classified(sol, L, N1=80, N2=110, R=None, lam_max=2.5, tol=2e-3, atol=3e-4, res_max=1e-2, **kw):
    """eigenvalues of the square regular system (Dirac + Eb, Erth, Etth, Etr; Ett, Etr, Etth, Erth at L = 1), each
    flagged by (i) agreement between resolutions N1, N2 and (ii) the residual of the Einstein components NOT imposed
    (relative, r < 0.8 R) on the N2 eigenvector. Returns a structured list of (lambda, converged, residual)."""
    ms = [NonRadialModes(sol, L, N=N, R=R, **kw) for N in (N1, N2)]
    # |lambda| < 1e-8: the static (lambda = 0) Hamiltonian-violating solutions of the square system, excluded
    l1 = ms[0].eigen(); l1 = l1[(np.abs(l1) < lam_max) & (np.abs(l1) > 1e-8)]
    l2, V2 = ms[1].eigen(vectors=True)
    s = (np.abs(l2) < lam_max) & (np.abs(l2) > 1e-8); l2 = l2[s]; V2 = V2[:, s]
    unused = [n for n in ('Ett', 'Err', 'Ea') if n not in ms[1].rows()]
    out = []
    for k, z in enumerate(l2):
        conv = len(l1) > 0 and np.min(np.abs(l1 - z)) < tol * abs(z) + atol
        res = ms[1].residuals(z, V2[:, k])
        out.append((z, conv, max(res[n] for n in unused)))
    out.sort(key=lambda t: -t[0].real)
    return out
