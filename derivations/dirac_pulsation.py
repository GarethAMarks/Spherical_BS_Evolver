"""
Linearised radial pulsation equations of spherical Dirac stars (two spinors, any self-interaction U(S)), in the polar-areal
gauge and the (-,+,+,+) signature, G = 1:

    ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2,    alpha = exp(phi),   alpha(t, r), X(t, r),
    psi_m = ( u(t,r) Omega_{-1,m},  i v(t,r) Omega_{+1,m} )   (u, v the complex orthonormal-frame amplitudes of ComplexSpinorField).

1. Field equations. With shift zero and g_T = 1, the extrinsic curvature has K_T = 0 and K_rr = -X X_t / alpha, so
   K = -X_t / (alpha X). The 3+1 Dirac equations (check_dirac_3p1.py, compact form, mu -> mu_eff = mu + U'(S) by
   check_dirac_selfint_3p1.py) become
       u_t = -X_t/(2X) u - i alpha mu_eff u + i [ (alpha/X)(v' + v/r) + alpha v/r + alpha' v/(2X) ]
       v_t = -X_t/(2X) v + i alpha mu_eff v - i [ (alpha/X)(u' + u/r) - alpha u/r + alpha' u/(2X) ]
   with S = (|u|^2 - |v|^2)/(2 pi), and the orthonormal-frame sources (P = S U' - U)
       rho = mu S + U + Srr_k + 2 Stt_k,  p_r = Srr_k + P,  j = Im(ubar u' + vbar v')/(2 pi X),
       Srr_k = Re(vbar u' - u vbar')/(2 pi X),  Stt_k = -Re(u vbar)/(2 pi r).
   The Einstein tensor of the metric is computed here, and G_tt = 8 pi alpha^2 rho, G_tr = 8 pi T_tr = -8 pi alpha X j,
   G_rr = 8 pi X^2 p_r give
       X'   = X (1 - X^2)/(2 r) + 4 pi r X^3 rho                 (Hamiltonian constraint)
       X_t  = -4 pi r alpha X^2 j                               (momentum constraint)
       phi' = (X^2 - 1)/(2 r) + 4 pi r X^2 p_r                   (polar slicing)
   The script checks the Einstein-tensor identities behind these.

2. Linearisation about a stationary star u = e^{-i w t} g(r), v = e^{-i w t} f(r) (f, g real, the physical amplitudes of
   DiracStar), X0(r), phi0(r). In the rotating frame
       u = e^{-i w t} (g + uR + i uI),  v = e^{-i w t} (f + vR + i vI),  X = X0 + dX,  phi = phi0 + dphi,
   with all perturbations ~ e^{lambda t}. The phase factor cancels from S, rho, p_r and j, and X_t is replaced by the
   momentum constraint. The result is
       lambda y_k = sum_j ( A_kj y_j + B_kj y_j' ),   k = uR, uI, vR, vI      (Dirac),
       dX'   = sum_j ( H_j y_j + H1_j y_j' ),     dphi' = sum_j ( Q_j y_j + Q1_j y_j' ),
   for y = (uR, uI, vR, vI, dX, dphi). The coefficients depend on r, w, the background (g, f, g', f', X0, phi0, phi0'),
   mu, and the potential through U, U', U'' at the background S. They are written to dirac_pulsation_coeffs.py, which
   relaxation/dirac_modes.py uses.

   The global phase rotation (uI = g, vI = f, the rest zero) is a zero mode, and the check on it is run here.
   Stability: lambda = -i sigma, with sigma the mode frequency in the rotating frame (the field oscillates at w +- sigma).
   Real sigma: stable oscillation. Re(lambda) > 0: instability, oscillatory if Im(lambda) != 0.
Run:  python3 derivations/dirac_pulsation.py   (a few seconds; writes relaxation/dirac_pulsation_coeffs.py)
"""
import os
import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi_', positive=True)
coords = [t, r, th, ph]
al = sp.Function('alpha')(t, r)
X = sp.Function('X')(t, r)

# ---------------- Einstein tensor of the polar-areal metric ----------------
g = sp.diag(-al**2, X**2, r**2, r**2 * sp.sin(th)**2)
gi = g.inv()
Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d]))
             for d in range(4)) / 2 for c in range(4)] for b in range(4)] for a in range(4)]
def Ric(b, c):
    return sp.simplify(sum(sp.diff(Gam[a][b][c], coords[a]) - sp.diff(Gam[a][b][a], coords[c])
                           + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(4))
                           for a in range(4)))
Rc = sp.Matrix(4, 4, lambda b, c: Ric(b, c) if b <= c else 0)
for b in range(4):
    for c in range(b):
        Rc[b, c] = Rc[c, b]
Rs = sp.simplify(sum(gi[a, b] * Rc[a, b] for a in range(4) for b in range(4)))
Gt = sp.simplify(Rc - Rs * g / 2)

# the three equations quoted above, for symbolic orthonormal-frame sources
rho, j, pr = sp.symbols('rho j p_r', real=True)
Xr = X * (1 - X**2) / (2 * r) + 4 * sp.pi * r * X**3 * rho
Xt = -4 * sp.pi * r * al * X**2 * j
phr = (X**2 - 1) / (2 * r) + 4 * sp.pi * r * X**2 * pr
chk_tt = sp.simplify((Gt[0, 0] - 8 * sp.pi * al**2 * rho).subs(sp.Derivative(X, r), Xr))
chk_tr = sp.simplify((Gt[0, 1] - 8 * sp.pi * (-al * X * j)).subs(sp.Derivative(X, t), Xt))
chk_rr = sp.simplify((Gt[1, 1] - 8 * sp.pi * X**2 * pr).subs(sp.Derivative(al, r), al * phr))
assert chk_tt == 0, chk_tt
assert chk_tr == 0, chk_tr
assert chk_rr == 0, chk_rr
print("1. Einstein equations: X' (Hamiltonian), X_t = -4 pi r alpha X^2 j (momentum), phi' (slicing): OK")

# ---------------- linearisation ----------------
eps, w, mu = sp.symbols('epsilon omega mu', real=True)
gb, fb, dgb, dfb, X0, ph0, dph0 = sp.symbols('g f dg df X0 phi0 dphi0', real=True)
Ub, dUb, ddUb = sp.symbols('U dU ddU', real=True)      # U, U', U'' at the background S
names = ['uR', 'uI', 'vR', 'vI', 'dX', 'dphi']
y = sp.symbols(' '.join(names), real=True)
yp = sp.symbols(' '.join(n + '_r' for n in names), real=True)
uR, uI, vR, vI, dX, dph = y
uRp, uIp, vRp, vIp, dXp, dphp = yp

I = sp.I
U = gb + eps * (uR + I * uI); V = fb + eps * (vR + I * vI)
Up = dgb + eps * (uRp + I * uIp); Vp = dfb + eps * (vRp + I * vIp)
Xs = X0 + eps * dX
phs = ph0 + eps * dph
phps = dph0 + eps * dphp
als = sp.exp(phs)
cj = sp.conjugate
Re = lambda e: (e + cj(e)) / 2
Im = lambda e: (e - cj(e)) / (2 * I)

S0 = (gb**2 - fb**2) / (2 * sp.pi)
S = sp.expand(Re(U * cj(U) - V * cj(V)) / (2 * sp.pi))
dS = sp.diff(S, eps).subs(eps, 0)
# potential expanded to first order about the background S
Ueps = Ub + eps * dUb * dS
dUeps = dUb + eps * ddUb * dS
mu_eff = mu + dUeps
P = S0 * dUb - Ub + eps * S0 * ddUb * dS          # P' = S U''
Srr = Re(cj(V) * Up - U * cj(Vp)) / (2 * sp.pi * Xs)
Stt = -Re(U * cj(V)) / (2 * sp.pi * r)
rho_e = mu * S + Ueps + Srr + 2 * Stt
pr_e = Srr + P
j_e = Im(cj(U) * Up + cj(V) * Vp) / (2 * sp.pi * Xs)
Xt_e = -4 * sp.pi * r * als * Xs**2 * j_e

# rotating frame: U_t = i w U + u_t e^{i w t}
RHSu = I * w * U - Xt_e / (2 * Xs) * U - I * als * mu_eff * U + I * ((als / Xs) * (Vp + V / r) + als * V / r + als * phps * V / (2 * Xs))
RHSv = I * w * V - Xt_e / (2 * Xs) * V + I * als * mu_eff * V - I * ((als / Xs) * (Up + U / r) - als * U / r + als * phps * U / (2 * Xs))

def lin(e):
    return sp.expand(sp.diff(e, eps).subs(eps, 0))
def zeroth(e):
    return sp.expand(e.subs(eps, 0))

dirac = [lin(Re(RHSu)), lin(Im(RHSu)), lin(Re(RHSv)), lin(Im(RHSv))]
ham = lin(Xs * (1 - Xs**2) / (2 * r) + 4 * sp.pi * r * Xs**3 * rho_e)
sli = lin((Xs**2 - 1) / (2 * r) + 4 * sp.pi * r * Xs**2 * pr_e)
# background equations (zeroth order) for the check of the zero mode
bg = [zeroth(Re(RHSu)), zeroth(Im(RHSu)), zeroth(Re(RHSv)), zeroth(Im(RHSv))]

# the background stationary equations: Im parts give g', f' (Re parts vanish identically for real f, g)
sol_bg = sp.solve([bg[1], bg[3]], [dgb, dfb], dict=True)[0]
assert all(sp.simplify(b) == 0 for b in (bg[0], bg[2]))
# zero mode: global phase, uI = g, vI = f, everything else zero
zm = {uR: 0, uI: gb, vR: 0, vI: fb, dX: 0, dph: 0, uRp: 0, uIp: dgb, vRp: 0, vIp: dfb, dXp: 0, dphp: 0}
for k, e in enumerate(dirac):
    assert sp.simplify(e.subs(zm).subs(sol_bg)) == 0, ('zero mode', k)
assert sp.simplify(ham.subs(zm)) == 0 and sp.simplify(sli.subs(zm)) == 0
print('2. linearised equations derived; the global phase rotation is a zero mode: OK')

# check against the DiracStar radial equations: G' = X[G/r + (mu_eff + w/alpha) F] with G = sqrt(alpha) r g
dg_ds = sp.simplify(sol_bg[dgb]); df_ds = sp.simplify(sol_bg[dfb])
a0 = sp.exp(ph0)
dg_expect = (X0 - 1) * gb / r - dph0 * gb / 2 + X0 * (mu + dUb + w / a0) * fb
df_expect = -(1 + X0) * fb / r - dph0 * fb / 2 + X0 * (mu + dUb - w / a0) * gb
assert sp.simplify(dg_ds - dg_expect) == 0 and sp.simplify(df_ds - df_expect) == 0
print('   the background equations are those of DiracStar / dirac_relax.py: OK')

# ---------------- coefficients ----------------
def coeffs(e):
    A = [sp.simplify(e.coeff(v)) for v in y]
    B = [sp.simplify(e.coeff(v)) for v in yp]
    rest = sp.simplify(e - sum(a * v for a, v in zip(A, y)) - sum(b * v for b, v in zip(B, yp)))
    assert rest == 0, rest
    return A, B
rows = {'uR': dirac[0], 'uI': dirac[1], 'vR': dirac[2], 'vI': dirac[3], 'ham': ham, 'sli': sli}
args = (r, w, mu, gb, fb, dgb, dfb, X0, ph0, dph0, Ub, dUb, ddUb)
out = ['"""generated by derivations/dirac_pulsation.py: coefficients of the linearised Dirac-star pulsation equations',
       '   row = sum_j A[j] y_j + B[j] y_j\',  y = (uR, uI, vR, vI, dX, dphi);  rows uR..vI give lambda y_k, ham gives dX\', sli gives dphi\'',
       '   arguments: r, w, mu, g, f, dg, df, X0, phi0, dphi0, U, dU, ddU (potential at the background S)"""',
       'import numpy as np', 'from numpy import exp, sqrt, pi', '']
for name, e in rows.items():
    A, B = coeffs(e)
    out.append('def %s(r, omega, mu, g, f, dg, df, X0, phi0, dphi0, U, dU, ddU):' % name)
    out.append('    z = np.zeros_like(r)')
    out.append('    A = [%s]' % ', '.join('(%s) + z' % sp.pycode(sp.simplify(a)).replace('math.', '') for a in A))
    out.append('    B = [%s]' % ', '.join('(%s) + z' % sp.pycode(sp.simplify(b)).replace('math.', '') for b in B))
    out.append('    return A, B')
    out.append('')
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'relaxation', 'dirac_pulsation_coeffs.py')
src = '\n'.join(out)
src = src.replace('epsilon', 'eps').replace('phi0', 'phi0')
open(path, 'w').write(src)
print('3. coefficients written to relaxation/dirac_pulsation_coeffs.py')
for name, e in rows.items():
    A, B = coeffs(e)
    print('   %-4s: depends on %s' % (name, ', '.join(n for n, a, b in zip(names, A, B) if a != 0 or b != 0)))
