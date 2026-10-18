"""
Brute-force check of the modified-cartoon reductions used by ComplexScalarField and ComplexVectorField.

The evolution code works on the line w^a = 0 of a (D-1)-dimensional spatial slice with SO(D-1) symmetry,
using the conformal variables chi, h_zz, h_ww (gamma_ij = h_ij / chi). Here we build an explicit SO(D-1)
symmetric configuration in full Cartesian coordinates (z, w^1, ..., w^{D-2}), compute every covariant
quantity by brute force, evaluate on w = 0, and compare with the reduced expressions implemented in C++:

  * the scalar-field terms of ComplexScalarField::rhs_CSF (as implemented in the repo), as a check of the method;
  * the Proca terms of ComplexVectorField::rhs_CVF: divergences D_i X^i, D_i E^i, Lie derivatives along the
    shift of a covector and a vector, X^i d_i alpha and D^i Z.

Smooth test profiles with the correct parity are used: scalars are even in r, radial components of vectors /
covectors are odd. Run:  python3 derivations/check_cartoon_reduction.py
"""
import sympy as sp

def check(D, z0=sp.Rational(7, 10), tol=1e-20):
    N = D - 1                       # spatial dimension
    n = D - 2                       # number of extra (w) directions
    xs = sp.symbols('x0:%d' % N, real=True)     # x0 = z, x1.. = w^a
    r = sp.sqrt(sum(x**2 for x in xs))

    # smooth SO(D-1)-symmetric profiles (generic coefficients); even functions for scalars, odd for radial components
    chi_f   = lambda s: 1 + sp.Rational(3, 10)*s**2 - sp.Rational(1, 20)*s**4
    hzz_f   = lambda s: 1 + sp.Rational(1, 5)*s**2 + sp.Rational(1, 10)*s**4
    hww_f   = lambda s: 1 + sp.Rational(1, 10)*s**2 - sp.Rational(1, 30)*s**4
    alpha_f = lambda s: sp.Rational(4, 5) + sp.Rational(1, 10)*s**2
    beta_f  = lambda s: sp.Rational(1, 20)*s - sp.Rational(1, 50)*s**3      # beta^i = beta(r) x^i / r
    phi_f   = lambda s: sp.Rational(1, 5) + sp.Rational(1, 10)*s**2 - sp.Rational(1, 40)*s**4
    P_f     = lambda s: sp.Rational(3, 10) - sp.Rational(1, 5)*s**2          # X_phi
    Xr_f    = lambda s: sp.Rational(2, 5)*s + sp.Rational(1, 10)*s**3        # X_i = Xr(r) x_i / r  (covector)
    Er_f    = lambda s: -sp.Rational(3, 10)*s + sp.Rational(1, 20)*s**3      # E^i = Er(r) x^i / r  (vector)
    Z_f     = lambda s: sp.Rational(1, 10)*s**2 + sp.Rational(1, 7)

    chi, hzz, hww = chi_f(r), hzz_f(r), hww_f(r)
    xhat = [x / r for x in xs]

    # physical metric gamma_ij = [h_ww delta_ij + (h_zz - h_ww) xhat_i xhat_j] / chi, and its inverse
    g = sp.Matrix(N, N, lambda i, j: (hww*(1 if i == j else 0) + (hzz - hww)*xhat[i]*xhat[j]) / chi)
    gi = sp.Matrix(N, N, lambda i, j: chi*((1 if i == j else 0)/hww + (1/hzz - 1/hww)*xhat[i]*xhat[j]))

    pt = {xs[0]: z0}
    pt.update({x: 0 for x in xs[1:]})
    ev = lambda e: sp.N(sp.sympify(e).subs(pt), 30)

    # numeric values at the point
    G = sp.Matrix(N, N, lambda i, j: ev(g[i, j]))
    GI = sp.Matrix(N, N, lambda i, j: ev(gi[i, j]))
    inv_err = (G*GI - sp.eye(N)).norm(); assert inv_err < 1e-20, ("inverse metric", inv_err)
    dg = [[[ev(sp.diff(g[i, j], xs[k])) for k in range(N)] for j in range(N)] for i in range(N)]
    Gam = [[[sum(sp.Rational(1, 2)*GI[k, l]*(dg[j][l][i] + dg[i][l][j] - dg[i][j][l]) for l in range(N))
             for j in range(N)] for i in range(N)] for k in range(N)]       # Gam[k][i][j] = Gamma^k_ij

    alpha, phi, P, Z = alpha_f(r), phi_f(r), P_f(r), Z_f(r)
    beta_up = [beta_f(r)*xh for xh in xhat]
    X_dn = [Xr_f(r)*xh for xh in xhat]
    E_up = [Er_f(r)*xh for xh in xhat]
    X_up = [sum(gi[i, j]*X_dn[j] for j in range(N)) for i in range(N)]

    # ---------- brute force ----------
    d = lambda f, i: ev(sp.diff(f, xs[i]))
    dd = lambda f, i, j: ev(sp.diff(f, xs[i], xs[j]))

    lap_phi = sum(GI[i, j]*(dd(phi, i, j) - sum(Gam[k][i][j]*d(phi, k) for k in range(N))) for i in range(N) for j in range(N))
    grad_a_grad_phi = sum(GI[i, j]*d(alpha, i)*d(phi, j) for i in range(N) for j in range(N))
    divX = sum(d(X_up[i], i) for i in range(N)) + sum(Gam[i][i][k]*ev(X_up[k]) for i in range(N) for k in range(N))
    divE = sum(d(E_up[i], i) for i in range(N)) + sum(Gam[i][i][k]*ev(E_up[k]) for i in range(N) for k in range(N))
    Xup_dalpha = sum(ev(X_up[i])*d(alpha, i) for i in range(N))
    DZ_up_z = sum(GI[0, j]*d(Z, j) for j in range(N))
    LieX_z = sum(ev(beta_up[j])*d(X_dn[0], j) + ev(X_dn[j])*d(beta_up[j], 0) for j in range(N))
    LieE_z = sum(ev(beta_up[j])*d(E_up[0], j) - ev(E_up[j])*d(beta_up[0], j) for j in range(N))
    # transverse components that the reduction assumes vanish
    transverse = [ev(X_dn[1]) if N > 1 else 0, ev(E_up[1]) if N > 1 else 0,
                  (sum(ev(beta_up[j])*d(X_dn[1], j) + ev(X_dn[j])*d(beta_up[j], 1) for j in range(N)) if N > 1 else 0)]

    # ---------- reduced (1D) expressions, exactly as in the C++ ----------
    s = sp.Symbol('s', positive=True)
    f1 = lambda F: (F(s).subs(s, z0), sp.diff(F(s), s).subs(s, z0), sp.diff(F(s), s, 2).subs(s, z0))
    (c, dc, _), (Hz, dHz, _), (Hw, dHw, _) = f1(chi_f), f1(hzz_f), f1(hww_f)
    (a, da, _), (b, db, _), (ph, dph, ddph) = f1(alpha_f), f1(beta_f), f1(phi_f)
    (Pv, dP, _), (Xv, dX, _), (Ev, dE, _), (Zv, dZ, _) = f1(P_f), f1(Xr_f), f1(Er_f), f1(Z_f)
    z = z0
    hZZ, hWW = 1/Hz, 1/Hw
    chris_Zzz = hZZ*sp.Rational(1, 2)*dHz                        # Spacetime::auxiliary_quantities_at_point
    chris_Zww = hZZ*(-sp.Rational(1, 2)*dHw + (Hz - Hw)/z)

    # ComplexScalarField::rhs_CSF: the gradient terms of the K_phi equation (multiplying -1/2) as implemented
    cD_zz = ddph - chris_Zzz*dph
    cD_ww = dph/z - chris_Zww*dph
    csf_terms = -sp.Rational(1, 2)*c*(hZZ*da*dph + a*(hZZ*cD_zz + n*hWW*cD_ww)) + sp.Rational(1, 4)*(D - 3)*a*hZZ*dc*dph
    csf_brute = -sp.Rational(1, 2)*(a*lap_phi + grad_a_grad_phi)

    # ComplexVectorField reductions
    gZZ = c*hZZ                                              # gamma^zz
    dgZZ = dc*hZZ - c*hZZ**2*dHz
    dlnsqrtg = sp.Rational(1, 2)*hZZ*dHz + sp.Rational(1, 2)*n*hWW*dHw - sp.Rational(1, 2)*(D - 1)*dc/c
    Xz_up, dXz_up = gZZ*Xv, dgZZ*Xv + gZZ*dX
    divX_red = dXz_up + n*Xz_up/z + Xz_up*dlnsqrtg
    divE_red = dE + n*Ev/z + Ev*dlnsqrtg
    results = {
        'scalar field Laplacian terms (as implemented)': (csf_terms, csf_brute),
        'D_i X^i': (divX_red, divX),
        'D_i E^i': (divE_red, divE),
        'X^i d_i alpha': (Xz_up*da, Xup_dalpha),
        'D^z Z': (gZZ*dZ, DZ_up_z),
        'Lie_beta X_z': (b*dX + Xv*db, LieX_z),
        'Lie_beta E^z': (b*dE - Ev*db, LieE_z),
    }
    ok = True
    print("D = %d (spatial dim %d), z = %.2f" % (D, N, z0))
    for k, (red, bf) in results.items():
        err = abs(sp.N(red, 30) - bf)
        ok &= err < tol
        print("  %-48s reduced = % .15f   brute = % .15f   |diff| = %.1e" % (k, float(red), float(bf), float(err)))
    terr = max(abs(float(t)) for t in transverse)
    ok &= terr < tol
    print("  %-48s max = %.1e" % ('transverse components (must vanish)', terr))
    return ok

if __name__ == '__main__':
    all_ok = all([check(4), check(5), check(6, z0=sp.Rational(13, 10))])
    print("\nALL CHECKS PASSED" if all_ok else "\nSOME CHECKS FAILED")
