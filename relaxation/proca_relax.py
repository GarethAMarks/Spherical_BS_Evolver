"""
Relaxation (collocation) solver for free spherical Proca stars, in the conventions of ProcaStar (signature (-,+,+,+),
G = 1, D = 4):

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,  alpha = exp(phi),  X^2 = 1 / (1 - 2 m / r),
    A = e^{-i w t} (f dt + i g dr),  E = f' - w g,  g = -w E / (alpha^2 mu^2),
    rho = [E^2 / (alpha^2 X^2) + mu^2 (f^2 / alpha^2 + g^2 / X^2)] / 2,
    p_r = [-E^2 / (alpha^2 X^2) + mu^2 (f^2 / alpha^2 + g^2 / X^2)] / 2,
    m' = 4 pi r^2 rho,  phi' = (X^2 - 1) / (2 r) + 4 pi r X^2 p_r,
    f' = E + w g,  E' = mu^2 f X^2 - 2 E / r + E (phi' + X'/X)        (Gauss law, the t component of the field equation),
    w' = 0 (frequency as a grid variable).
Boundary conditions: f(r0) / alpha(r0) = A (A_central of ProcaStar), E(r0) = mu^2 f(r0) r0 / 3, m(r0) = 0;
phi(R) = ln(1 - 2 m(R)/R) / 2 (alpha -> 1), and decay E(R) + mu^2 f(R) / k = 0, k = sqrt(mu^2 - w^2 / alpha(R)^2).
seed_from_PSdata reads ProcaStar's PSdata.dat (r f g X phi E m).
"""
import os, re, subprocess
import numpy as np
from scipy.integrate import solve_bvp


class ProcaRelax:
    def __init__(self, mu=1., r0=1e-4, R=60., tol=1e-8, max_nodes=300000):
        self.mu, self.r0, self.R, self.tol, self.max_nodes = mu, r0, R, tol, max_nodes

    def terms(self, r, y):
        f, E, m, phi, w = y
        mu = self.mu
        al = np.exp(phi)
        X = 1. / np.sqrt(np.maximum(1. - 2. * m / r, 1e-300))
        g = -w * E / (al**2 * mu**2)
        Et = E**2 / (al**2 * X**2)
        V = mu**2 * (f**2 / al**2 + g**2 / X**2)
        rho = 0.5 * (Et + V); pr = 0.5 * (-Et + V)
        dX = X * (1 - X**2) / (2 * r) + 4 * np.pi * r * X**3 * rho
        dphi = (X**2 - 1) / (2 * r) + 4 * np.pi * r * X**2 * pr
        return al, X, g, rho, pr, dX, dphi

    def rhs(self, r, y):
        f, E, m, phi, w = y
        al, X, g, rho, pr, dX, dphi = self.terms(r, y)
        df = E + w * g
        dE = self.mu**2 * f * X**2 - 2 * E / r + E * (dphi + dX / X)
        return np.vstack([df, dE, 4 * np.pi * r**2 * rho, dphi, 0 * r])

    def bc(self, ya, yb, A):
        mu, r0, R = self.mu, self.r0, self.R
        alR = np.exp(yb[3]); k = np.sqrt(max(mu**2 - (yb[4] / alR)**2, 1e-14))
        return np.array([ya[0] - A * np.exp(ya[3]), ya[1] - mu**2 * ya[0] * r0 / 3., ya[2],
                         yb[3] - 0.5 * np.log(1. - 2. * yb[2] / R), yb[1] + mu**2 * yb[0] / k])

    def solve(self, A, r, y):
        return solve_bvp(self.rhs, lambda a, b: self.bc(a, b, A), r, y, tol=self.tol, max_nodes=self.max_nodes)

    def diagnostics(self, sol, A):
        f, E, m, phi, w = sol.y
        al, X, g, rho, pr, dX, dphi = self.terms(sol.x, sol.y)
        M = m[-1]
        j = np.searchsorted(m, 0.99 * M)
        r99 = np.interp(0.99 * M, m[max(j - 1, 0):j + 1], sol.x[max(j - 1, 0):j + 1])
        # Noether charge: Q = int sqrt(-g) j^t, j^a = (i/2)(Fbar^{ab} A_b - F^{ab} Abar_b) -> 4 pi int r^2 X E (-g) ... use
        # the standard form for this ansatz: Q = 4 pi int r^2 (w f / alpha - ... ) is not needed here
        return dict(A=A, omega=w[0], M=M, r99=r99, C=M / r99, alpha0=np.exp(phi[0]),
                    f_nodes=int(np.sum(np.sign(f[1:]) != np.sign(f[:-1]))), E_nodes=int(np.sum(np.sign(E[1:-5]) != np.sign(E[:-6]))))


def seed_from_PSdata(path, r_mesh):
    d = np.genfromtxt(path)
    d = d[np.isfinite(d).all(axis=1)]
    rr = d[:, 0]
    return np.vstack([np.interp(r_mesh, rr, d[:, 1]), np.interp(r_mesh, rr, d[:, 5]), np.interp(r_mesh, rr, d[:, 6]),
                      np.interp(r_mesh, rr, d[:, 4]), np.zeros_like(r_mesh)])


def shoot(binary, par_template, work, A, eigen=0, R=60., n=12001, fguess=0.9):
    os.makedirs(work, exist_ok=True)
    par = open(par_template).read()
    sub = {'dirac': 0, 'proca': 1, 'lambda': 0, 'solitonic': 0, 'A_central': A, 'cycle_only': 0, 'stop_time': 0, 'eigen': eigen,
           'R': R, 'n_gridpoints': n, 'frequency_guess': fguess, 'perturb': 0, 'isotropic': 0, 'BS_resolution_factor': 1}
    for k, v in sub.items():
        par = re.sub(r'^%s = [^\n#]*' % k, '%s = %s ' % (k, v), par, flags=re.M)
    open(os.path.join(work, 'BSParams.par'), 'w').write(par)
    out = subprocess.run([os.path.abspath(binary)], cwd=work, capture_output=True, text=True).stdout
    w = float(re.search(r'Final frequency after enforcing lapse condition is ([0-9.eE+-]+)', out).group(1))
    return os.path.join(work, 'PSdata.dat'), w


def background(binary, par_template, A, R=60., work='proca_shoot', eigen=0):
    path, w = shoot(binary, par_template, work, A, eigen=eigen, R=R)
    s = ProcaRelax(R=R)
    r = np.linspace(s.r0, R, 6001)
    y = seed_from_PSdata(path, r); y[4] = w
    sol = s.solve(A, r, y)
    assert sol.status == 0, sol.message
    return s, sol, s.diagnostics(sol, A)
