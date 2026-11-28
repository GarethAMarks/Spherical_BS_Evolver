import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, time, pickle, sys
from emkg2 import MM, to_mm, make_xi, NV
from nat2 import start

def Z(mm, Y, D, lam, wl, wD):
    return np.concatenate([Y.ravel()/np.sqrt(mm.N), [D*wD, lam*wl]])

def corr(mm, Y, D, lam, Zp, tz, ds, wl, wD, tol=1e-9, maxit=15):
    sY = 1/np.sqrt(mm.N)
    def arcres(Y, D, lam): return tz @ (Z(mm, Y, D, lam, wl, wD) - Zp) - ds
    for it in range(maxit):
        R0 = mm.residual(Y, D, lam); g = arcres(Y, D, lam)
        nr = max(np.abs(R0).max(), abs(g))
        if not np.isfinite(nr): return Y, D, lam, False, it
        if nr < tol: return Y, D, lam, True, it
        J, dl = mm.jacobian(Y, D, lam, R0)
        row = np.concatenate([tz[:-2]*sY, [tz[-2]*wD]]); rl = tz[-1]*wl
        lu = spla.splu(J.tocsc())
        x1 = lu.solve(-R0); x2 = lu.solve(-dl)
        den = row @ x2 + rl
        dlam = (-g - row @ x1)/den
        dU = x1 + dlam*x2
        dY, dD = mm.unpack(dU)
        n0 = np.linalg.norm(np.concatenate([R0, [g]])); t = 1.0
        while t > 1e-5:
            Yt, Dt, lt = Y + t*dY, D + t*dD, lam + t*dlam
            if mm.safe(Yt, Dt):
                n1 = np.linalg.norm(np.concatenate([mm.residual(Yt, Dt, lt), [arcres(Yt, Dt, lt)]]))
                if np.isfinite(n1) and n1 < n0*(1 - 1e-4*t): break
            t *= 0.5
        if t <= 1e-5: return Y, D, lam, False, it
        Y, D, lam = Yt, Dt, lt
    return Y, D, lam, False, maxit

def interp_phys(Y, rfrom, rto):
    return np.stack([np.interp(rto, rfrom, Y[:, k]) for k in range(NV)], axis=1)

def arc(s0, qh, rmax, W2a=25.0, nsteps=4000, wl=5.0, wD=0.2, dsmax=0.05, fac=0.08, hmax=0.03, hc=0.1, stop=None, tag='', verbose=True, init=None, save=True, direction=1, snap=None):
    # snap(o, mm, (Y, D, lam)): optional callback after every accepted step (used by make_charged_id.py to keep profiles)
    if init is None:
        mm, Y, D, ok = start(s0, qh, rmax, W2a, hc=hc)
        lam = np.sqrt(W2a)
    else:
        mm, Y, D, lam = init
    # second point: slightly smaller Wc (towards compact) for direction=1
    l1 = lam*np.exp(-2e-3*direction)
    Y1, D1, ok, _ = mm.solve(Y, D, l1)
    pts = [(Y, D, lam), (Y1, D1, l1)]
    out = [mm.observables(Y, D), mm.observables(Y1, D1)]
    for o, p in zip(out, pts): o['Wc2'] = p[2]**2
    ds = 0.005; t0 = time.time(); hfc = mm.xi[mm.kw+1]-mm.xi[mm.kw]
    for step in range(nsteps):
        (Ya, Da, la), (Yb, Db, lb) = pts[-2], pts[-1]
        za, zb = Z(mm, Ya, Da, la, wl, wD), Z(mm, Yb, Db, lb, wl, wD)
        tz = zb - za; tz /= np.linalg.norm(tz)
        zp = zb + ds*tz
        Yp = zp[:-2].reshape(Yb.shape)*np.sqrt(mm.N); Dp = zp[-2]/wD; lp = zp[-1]/wl
        Yn, Dn, ln, ok, it = corr(mm, Yp, Dp, lp, zb, tz, ds, wl, wD)
        if not ok:
            ds *= 0.5
            if ds < 1e-7: print("  ds underflow"); break
            continue
        pts.append((Yn, Dn, ln)); pts = pts[-2:]
        o = mm.observables(Yn, Dn); o['Wc2'] = ln**2; out.append(o)
        u = np.sqrt(max(o['fR_out'], 1e-12))
        if abs(Dn) > 1.5 or min(hmax, fac*u) < 0.7*hfc:
            # recentre both points onto a new mesh (physical-space interpolation), re-correct the newest
            hf = min(hmax, fac*u); xi, kw = make_xi(o['R'], hf, rmax, hc=hc)
            mn = MM(s0, qh*s0, xi, kw)
            ra, rb = mm.r_of(pts[0][1]), mm.r_of(pts[1][1])
            Yan = interp_phys(pts[0][0], ra, xi); Ybn = interp_phys(pts[1][0], rb, xi)
            za, zb = Z(mn, Yan, 0.0, pts[0][2], wl, wD), Z(mn, Ybn, 0.0, pts[1][2], wl, wD)
            tz = zb - za; tz /= np.linalg.norm(tz)
            Yc, Dc, lc, okc, _ = corr(mn, Ybn, 0.0, pts[1][2], zb, tz, 0.0, wl, wD, maxit=30)
            if not okc:
                Yc, Dc, okc, _ = mn.solve(Ybn, 0.0, pts[1][2]); lc = pts[1][2]
            if not okc: print("  recentre failed"); break
            shiftY = Yc - Ybn
            pts = [(Yan + shiftY, Dc, pts[0][2] + (lc - pts[1][2])), (Yc, Dc, lc)]
            mm = mn; hfc = hf
        ds = min(dsmax, ds*(1.4 if it <= 3 else (1.0 if it <= 6 else 0.6)))
        if verbose and len(out) % 50 == 0:
            print("  %4d Wc2=%.5f Mhat=%.5f C=%.4f Q/M=%.5f ac=%.4f uR=%.4f N=%d ds=%.3g (%.0fs)"%(len(out), o['Wc2'], o['Mhat'], o['C'], o['QM'], o['alpha_c'], u, mm.N, ds, time.time()-t0), flush=True)
        if snap: snap(o, mm, pts[-1])
        if o['R'] > rmax-25: print("  rmax reached"); break
        if stop and stop(o, out): break
    Yl, Dl, ll = pts[-1]
    if not save: return out, (mm, Yl, Dl, ll)
    pickle.dump(out, open('arc2_%s_%s%s.pkl'%(s0, qh, tag), 'wb'))
    pickle.dump(dict(Y=Yl, r=mm.r_of(Dl), lam=ll, s0=s0, qh=qh), open('arc2end_%s_%s%s.pkl'%(s0, qh, tag), 'wb'))
    return out

def fmin_ext(o):
    M, Q, R = o['M'], o['Q'], o['R']
    rs = Q**2/M
    return (1 - 2*M/R + Q**2/R**2) if rs <= R else (1 - M**2/Q**2)

def classify(s0, qh, rmax, hc=0.1, fac=0.08, hmax=0.03, tag='', verbose=False, fend=1e-4):
    st = dict(Mmax=0)
    def stop(o, out):
        o['fmin'] = fmin_ext(o)
        if o['Mhat'] > st['Mmax']: st['Mmax'] = o['Mhat']
        if o['Mhat'] < 0.97*st['Mmax']: return True
        return o['fmin'] < fend
    t0 = time.time()
    out = arc(s0, qh, rmax, stop=stop, hc=hc, fac=fac, hmax=hmax, tag=tag, verbose=verbose)
    M = np.array([o['Mhat'] for o in out]); i = int(np.argmax(M)); o = out[i]; L = out[-1]
    status = 'EXTREMAL-END' if L.get('fmin', 1) < fend and i >= len(out)-3 else ('FOLD' if L['Mhat'] < 0.97*M.max() else 'INCOMPLETE')
    print("s0=%.3f qh=%.3f %-12s Mmax: Mhat=%.5f C=%.4f Q/M=%.5f ac=%.4f fmin=%.2e Wc2=%.4g | last Q/M=%.5f C=%.4f (n=%d, %.0fs)"%(
        s0, qh, status, o['Mhat'], o['C'], o['QM'], o['alpha_c'], fmin_ext(o), o['Wc2'], L['QM'], L['C'], len(out), time.time()-t0), flush=True)
    return status, out

if __name__ == "__main__":
    for a in sys.argv[1:]:
        p = a.split(':'); fine = len(p) > 3 and p[3] == 'fine'
        classify(float(p[0]), float(p[1]), float(p[2]), **(dict(hc=0.05, fac=0.04, hmax=0.02, tag='_fine') if fine else {}))
