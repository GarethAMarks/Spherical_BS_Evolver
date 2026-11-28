"""
Leading-order (sigma_0 -> 0) charged thin-shell solitonic boson stars.

Action (G = c = 1, Gaussian Maxwell units, mu = 1):
  S = int sqrt(-g) [ R/16pi - (1/2)(|D phi|^2 + V) - F^2/16pi ],
  D = nabla - i q A,  V = mu^2 |phi|^2 (1 - 2|phi|^2/sigma_0^2)^2.
Static ansatz: ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,
  phi = A(r) e^{-i omega t},  A_mu dx^mu = Phi(r) dt.
Double-scaling limit: sigma_0 -> 0 with
  x = sigma_0^2 r, m = sigma_0^2 m_phys, Q = sigma_0^2 Q_phys,
  W = Omega/sigma_0 (Omega = (omega + q Phi)/alpha local frequency),
  qh = q / sigma_0  (fixed).
Core (A pinned at v = sigma_0/sqrt2): stiff charged condensate,
  rho = p = W^2/4,  rho_e = qh W / 2   (units of sigma_0^4).
ODEs (' = d/dx):
  m' = pi x^2 W^2 + Q^2/(2 x^2)
  Q' = 2 pi qh x^2 X W
  W' = qh X Q / x^2 - W F,   F = alpha'/alpha = X^2 (m + pi x^3 W^2 - Q^2/(2x)) / x^2
Wall: pure-tension membrane, Sigma = sigma_0^2/8 (mu=1). Israel:
  (J1) u_- - u_+ = pi x_s / 2
  (J2) [ (1-u_+^2)/2 - Q^2/(2x^2) ] / u_+ - [ m_-/x + pi x^2 W^2 - Q^2/(2x^2) ] / u_- + pi x/2 = 0
Exterior: Reissner-Nordstrom, f = 1 - 2M/r + Q^2/r^2.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PI = np.pi


def rhs(x, y, qh):
    m, Q, W, lna = y
    X2 = 1.0 / (1.0 - 2.0 * m / x)
    X = np.sqrt(X2)
    F = X2 * (m + PI * x**3 * W**2 - Q**2 / (2 * x)) / x**2
    return [PI * x**2 * W**2 + Q**2 / (2 * x**2),
            2 * PI * qh * x**2 * X * W,
            qh * X * Q / x**2 - W * F,
            F]


def core(Wc, qh, xmax=5.0):
    """Integrate the core from the centre until J1 can no longer be met."""
    x0 = 1e-6 / max(Wc, 1e-6)
    y0 = [PI / 3 * Wc**2 * x0**3, 2 * PI / 3 * qh * Wc * x0**3, Wc, 0.0]

    def ev(x, y, qh):
        m = y[0]
        a = 1 - 2 * m / x
        return (np.sqrt(a) if a > 0 else -1.0) - PI * x / 2
    ev.terminal = True
    ev.direction = -1
    sol = solve_ivp(rhs, (x0, xmax), y0, args=(qh,), rtol=1e-11, atol=1e-14,
                    dense_output=True, events=ev, method='DOP853')
    return sol


def J2(x, sol, qh):
    m, Q, W, lna = sol.sol(x)
    um = np.sqrt(1 - 2 * m / x)
    up = um - PI * x / 2
    if up <= 0:
        return np.nan
    return (((1 - up**2) / 2 - Q**2 / (2 * x**2)) / up
            - (m / x + PI * x**2 * W**2 - Q**2 / (2 * x**2)) / um + PI * x / 2)


def stars_for(Wc, qh, ngrid=1500):
    """All equilibrium wall positions for a given central W_c and charge qh."""
    sol = core(Wc, qh)
    xe = sol.t[-1]
    xs = np.linspace(sol.t[0] * 10, xe * (1 - 1e-9), ngrid)
    g = np.array([J2(x, sol, qh) for x in xs])
    out = []
    for i in range(len(xs) - 1):
        if np.isfinite(g[i]) and np.isfinite(g[i + 1]) and g[i] * g[i + 1] < 0:
            xr = brentq(J2, xs[i], xs[i + 1], args=(sol, qh), xtol=1e-14)
            out.append(props(xr, sol, qh, Wc))
    return out


def props(x, sol, qh, Wc):
    m, Q, W, lna = sol.sol(x)
    um = np.sqrt(1 - 2 * m / x)
    up = um - PI * x / 2
    M = (x * (1 - up**2) + Q**2 / x) / 2
    # lapse normalisation: alpha continuous at wall, alpha_+ = up
    lna_shift = np.log(up) - lna
    omega = W * up + qh * Q / x
    # Noether charge N (scaled: sigma0^3 N) = Q/qh, or direct integral when qh = 0
    if qh > 0:
        N = Q / qh
    else:
        xx = np.linspace(sol.t[0], x, 4000)
        mm, QQ, WW, _ = sol.sol(xx)
        N = np.trapezoid(2 * PI * xx**2 * WW / np.sqrt(1 - 2 * mm / xx), xx)
    # light rings: H = alpha^2/r^2. interior extrema where x F = 1
    xx = np.linspace(sol.t[0] * 10, x, 4000)
    yy = sol.sol(xx)
    mm, QQ, WW, _ = yy
    X2 = 1 / (1 - 2 * mm / xx)
    F = X2 * (mm + PI * xx**3 * WW**2 - QQ**2 / (2 * xx)) / xx**2
    s = xx * F - 1
    inner = [xx[i] for i in range(len(xx) - 1) if s[i] * s[i + 1] < 0]
    # exterior photon sphere of RN
    disc = 9 * M**2 - 8 * Q**2
    rph = (3 * M + np.sqrt(disc)) / 2 if disc >= 0 else np.nan
    # slopes of H just inside/outside the wall (kink)
    Fin = (m + PI * x**3 * W**2 - Q**2 / (2 * x)) / (x**2 * um**2)
    Fout = (M - Q**2 / x) / (x**2 * up**2)
    alpha0 = np.exp(lna_shift)  # central lapse
    return dict(qh=qh, Wc=Wc, x=x, M=M, Q=Q, N=N, omega=omega, Ws=W, um=um, up=up,
                C=M / x, QM=Q / M if M > 0 else np.nan, rph=rph, inner_LR=inner,
                Fin=Fin, Fout=Fout, alpha_c=alpha0, ys=Wc * x,
                kappa=qh / Wc)
