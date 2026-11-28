import numpy as np, sys, pickle, time
from emkg import Model, guess_from_LO
from emkg2 import MM, to_mm
def start(s0, qh, rmax, W2, hf=0.03, hc=0.1):
    q=qh*s0
    m=Model(s0,q,rmax=rmax,h=0.1,mode='w0')
    Y,P=guess_from_LO(m,np.sqrt(W2),qh,pick='min',xmin=6*s0**2)
    Y,ok,it=m.solve(Y,Y[0,4],maxit=300)
    o=m.observables(Y)
    mm,Yn=to_mm(m,Y,m.r,s0,q,o['R'],hf,rmax,hc)
    Yn,D,ok,it=mm.solve(Yn,0.0,np.sqrt(W2))
    return mm,Yn,D,ok
def recentre(mm,Y,D,lam,fac,hmax,hc,rmax):
    o=mm.observables(Y,D); u=np.sqrt(max(o['fR_out'],1e-12)); hf=min(hmax,fac*u)
    r=mm.r_of(D)
    mn,Yn=to_mm(None,Y,r,mm.s0,mm.q,o['R'],hf,rmax,hc)
    Yn,Dn,ok,it=mn.solve(Yn,0.0,lam)
    return mn,Yn,Dn,ok,hf
def natural(s0,qh,rmax,W2a=25.0,W2b=0.005,dlog=0.01,fac=0.08,hmax=0.03,hc=0.1,tag='',verbose=True):
    mm,Y,D,ok=start(s0,qh,rmax,W2a,hc=hc)
    assert ok
    lam=np.sqrt(W2a); rows=[]; t0=time.time(); hfc=mm.xi[mm.kw+1]-mm.xi[mm.kw]
    step=dlog; Yprev=None
    while lam**2>W2b:
        lt=lam*np.exp(-step)
        Yp=Y if Yprev is None else Y+(Y-Yprev[0])*(step/Yprev[2])
        Dp=D if Yprev is None else D+(D-Yprev[1])*(step/Yprev[2])
        Yn,Dn,ok,it=mm.solve(Yp,Dp,lt)
        if not ok:
            Yn,Dn,ok,it=mm.solve(Y,D,lt)
        if not ok:
            step*=0.5
            if step<1e-5: print("  stuck at Wc2=%.5g"%lam**2); break
            continue
        Yprev=(Y,D,step); Y,D,lam=Yn,Dn,lt
        o=mm.observables(Y,D); o['Wc2']=lam**2; rows.append(o)
        u=np.sqrt(max(o['fR_out'],1e-12))
        if abs(D)>1.0 or min(hmax,fac*u)<0.7*hfc:
            mm,Y,D,ok,hfc=recentre(mm,Y,D,lam,fac,hmax,hc,rmax); Yprev=None
            if not ok: print("  recentre failed at Wc2=%.5g"%lam**2); break
        if o['R']>rmax-25: print("  rmax reached"); break
        step=min(0.05,step*(1.3 if it<=4 else 0.8))
        if verbose and len(rows)%40==0: print("   Wc2=%8.4g Mhat=%.5f C=%.4f Q/M=%.5f ac=%.4f uR=%.4f N=%d (%.0fs)"%(o['Wc2'],o['Mhat'],o['C'],o['QM'],o['alpha_c'],u,mm.N,time.time()-t0),flush=True)
    pickle.dump(rows,open('nat2_%s_%s%s.pkl'%(s0,qh,tag),'wb'))
    print("s0=%.3f qh=%.2f: %d pts down to Wc2=%.4g (%.0fs)"%(s0,qh,len(rows),lam**2,time.time()-t0),flush=True)
    return rows
if __name__=="__main__":
    a=sys.argv
    natural(float(a[1]),float(a[2]),float(a[3]))
