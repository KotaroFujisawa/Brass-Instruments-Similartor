"""Causal two-way nonlinear waveguide coupled to Adachi-Sato 2D lips.

Each time sample scatters waves at EVERY area junction, propagates both ways,
applies causal wall losses, and applies a weakly nonlinear Burgers correction.
The actual incoming wave is used in the lip/flow boundary equation. No stored
excitation or fixed reflection impulse is used in the coupled simulation.
"""
from trumpet import *
from nonlinear_propagation import SOUND_DIFFUSIVITY


def loss_fit(terms=8):
    """Positive diffusive fit sqrt(s) ~ sum w_j*s/(s+theta_j)."""
    theta=2*np.pi*np.geomspace(3,160000,terms)
    f=np.geomspace(40,25000,300);s=2j*np.pi*f
    mat=s[:,None]/(s[:,None]+theta)/np.sqrt(s[:,None])
    weights,_=optimize.nnls(np.r_[mat.real,mat.imag],np.r_[np.ones(len(f)),np.zeros(len(f))])
    err=mat@weights-1
    return theta,weights,float(np.max(abs(err)))


def setup(valves=(),fs=192000,depth=None,losses=True):
    if depth is None:depth=json.loads((ROOT/'data/model.json').read_text())['cup_depth_m']
    _,_,x,r=geometry(tuple(valves),depth)
    dx=C/fs;n=int(round(x[-1]/dx))
    centers=(np.arange(n)+.5)*dx
    rad=np.interp(centers,x,r);area=np.pi*rad**2
    reflection=(area[:-1]-area[1:])/(area[:-1]+area[1:])
    theta,weights,error=loss_fit()
    dt=1/fs
    a=(2-theta*dt)/(2+theta*dt);b=theta*dt/(2+theta*dt)
    # A short passive section: H=1/(1+K*sum w*s/(s+theta)).
    # K = dx*(2*sqrt(pi)*1.045)/(c*632.8*r).
    k=dt*(2*np.sqrt(np.pi)*1.045/632.8)/rad
    if not losses:k*=0
    den=1+k*np.sum(weights*(1-b))
    return dict(fs=fs,dx=dx,radius=rad,area=area,reflection=reflection,
        theta=theta,weights=weights,a=a,b=b,loss_k=k,loss_den=den,
        z0=RHO*C/area[0],zbell=RHO*C/area[-1],length=n*dx,
        loss_fit_max_relative_error=error,cup_depth_m=depth,
        viscosity=SOUND_DIFFUSIVITY/2 if losses else 0.)


@njit(cache=True)
def _mc(q,sl):
    n=len(q);sl[0]=0.;sl[n-1]=0.
    for i in range(1,n-1):
        a=q[i]-q[i-1];b=q[i+1]-q[i]
        if a*b<=0:sl[i]=0.
        else:sl[i]=np.sign(a)*min(2*abs(a),2*abs(b),abs(a+b)/2)


@njit(cache=True)
def _nl_rhs(q,rad,sign,sl,flux,rhs):
    _mc(q,sl)
    for i in range(len(q)-1):
        left=q[i]+.5*sl[i];right=q[i+1]-.5*sl[i+1]
        speed=max(abs(left),abs(right))
        flux[i]=.25*sign*(left*left+right*right)-.5*speed*(right-left)
    rhs[0]=0.;rhs[-1]=0.
    for i in range(1,len(q)-1):rhs[i]=-(flux[i]-flux[i-1])/rad[i]


@njit(cache=True)
def nonlinear_correction(p,mean,rad,sign,beta,ratio,nu_ratio,q,tmp,sl,flux,rhs):
    """Exact linear travel is split from nonlinear speed correction.

    q=a*p_ac removes leading geometric spreading. Flux is +/-q^2/2,
    coefficient beta/(rho*c*a). SSP-RK2 + MUSCL-MC entropy flux.
    Mean flow is carried by the linear network, excluded from acoustic distortion.
    """
    n=len(p)
    for i in range(n):q[i]=rad[i]*(p[i]-mean[i])
    coeff=beta/(RHO*C)*ratio
    _nl_rhs(q,rad,sign,sl,flux,rhs)
    for i in range(n):tmp[i]=q[i]+coeff*rhs[i]
    _nl_rhs(tmp,rad,sign,sl,flux,rhs)
    for i in range(1,n-1):
        updated=.5*(q[i]+tmp[i]+coeff*rhs[i])
        updated+=nu_ratio*(q[i+1]-2*q[i]+q[i-1])
        p[i]+=(updated-q[i])/rad[i]


@njit(cache=True)
def volume_loss(p,mean,rad,nu_ratio,work):
    for i in range(len(p)):work[i]=rad[i]*(p[i]-mean[i])
    for i in range(1,len(p)-1):
        p[i]+=nu_ratio*(work[i+1]-2*work[i]+work[i-1])/rad[i]


@njit(cache=True)
def run_core(rad,refl,weights,a,b,loss_k,den,fs,seconds,flip,pressure,yeq,beta,
             hold,viscosity,source,feedback=True,matched_end=False,mean_hz=20.):
    """source empty => self-excited lips; otherwise source specifies outgoing p+.

    feedback=False is a diagnostic ablation: the return travels normally, but
    is removed ONLY from the lip pressure boundary. It is never used for delivery.
    """
    cells=len(rad);steps=int(seconds*fs);terms=len(weights);dt=1/fs;dx=C/fs
    f=np.zeros(cells);g=np.zeros(cells);fn=np.zeros(cells);gn=np.zeros(cells)
    mf=np.zeros(cells);mg=np.zeros(cells)
    vf=np.zeros((cells,terms));vg=np.zeros((cells,terms))
    yf=np.zeros(cells);yg=np.zeros(cells)
    q=np.zeros(cells);tmp=q.copy();sl=q.copy();flux=q.copy();rhs=q.copy()
    z0=RHO*C/(np.pi*rad[0]**2);zb=RHO*C/(np.pi*rad[-1]**2)
    rb=.61**2/.25;tau=.61*rad[-1]/C
    rate=rb/((rb+1)*tau);ri=(rb-1)/(rb+1)
    ra=(2-rate*dt)/(2+rate*dt);rbf=rate*dt/(2+rate*dt)
    bellstate=0.;bellprev=0.
    k=1.5*flip;m=1.5/(4*np.pi*np.pi*flip);omega=2*np.pi*flip
    xx=.001;yy=yeq;xp=xx;yp=yy;p=0.;plip=0.;ua=0.
    output=np.zeros((steps,12))
    mean_rate=1-np.exp(-2*np.pi*mean_hz*dt)
    maxmach=0.
    for it in range(steps):
        # Waves arriving at boundaries, including all prior nonlinear return trips.
        incoming=g[0]
        h=2*incoming if feedback else 0.
        t=it/fs
        env=1.
        if t<.045:env=.5-.5*np.cos(np.pi*t/.045)
        elif t>hold:env=.5+.5*np.cos(np.pi*min(1.,(t-hold)/.065))
        pm=pressure*env
        if len(source)>0:
            outgoing=source[it];p=outgoing+incoming
            u=(outgoing-incoming)/z0;ua=u;ul=0.
        else:
            quality=3. if yy>0 else .5;damp=omega/quality
            fx=2*B/m*(pm-p)*(.004-yy)
            fy=2*B/m*((pm-p)*xx+D*plip)
            if yy<0:fy-=6*k/m*yy
            xn=(2*xx-(1-damp*dt/2)*xp+dt*dt*(fx-omega*omega*(xx-.001)))/(1+damp*dt/2)
            yn=(2*yy-(1-damp*dt/2)*yp+dt*dt*(fy-omega*omega*(yy-yeq)))/(1+damp*dt/2)
            vx=(xn-xp)/(2*dt);vy=(yn-yp)/(2*dt)
            ul=B*(xx*vy-(yy-.004)*vx)
            xp,yp=xx,yy;xx,yy=xn,yn
            area=max(2*B*yy,0.)
            if area>1e-12:
                aa=RHO*(.5/area**2-1/(SCUP*area)+1/SCUP**2)
                inert=RHO*D/(area*dt);bb=inert+z0
                cc=pm-h-z0*ul+inert*ua
                ua=2*cc/(bb+np.sqrt(bb*bb+4*aa*abs(cc)))
                p=z0*(ua+ul)+h
                plip=p-RHO*ua*abs(ua)*(1/(SCUP*area)-1/SCUP**2)
            else:
                ua=0.;p=z0*ul+h;plip=p
            u=ua+ul;outgoing=z0*u+(incoming if feedback else 0.)
        if not np.isfinite(p) or abs(p)>1e6 or abs(xx)>.05 or abs(yy)>.02:
            raise ValueError('Coupled boundary diverged; no clipping applied')
        # Passive frequency-dependent bell reflection. Tustin realization.
        bellstate=ra*bellstate+rbf*(f[-1]+bellprev)
        bellprev=f[-1]
        returned=ri*f[-1]-(1+ri)*bellstate
        if matched_end:returned=0.
        bellflow=(f[-1]-returned)/zb
        bellp=f[-1]+returned
        fn[0]=outgoing;gn[-1]=returned
        for j in range(cells-1):
            rr=refl[j]
            fn[j+1]=(1+rr)*f[j]-rr*g[j+1]
            gn[j]=rr*f[j]+(1-rr)*g[j+1]
        # Positive-real per-cell wall filters, fixed in space, both directions.
        for j in range(cells):
            sf=0.;sg=0.
            for ell in range(terms):
                sf+=weights[ell]*(a[ell]*vf[j,ell]+b[ell]*yf[j])
                sg+=weights[ell]*(a[ell]*vg[j,ell]+b[ell]*yg[j])
            ff=(fn[j]+loss_k[j]*sf)/den[j]
            gg=(gn[j]+loss_k[j]*sg)/den[j]
            for ell in range(terms):
                vf[j,ell]=a[ell]*vf[j,ell]+b[ell]*(ff+yf[j])
                vg[j,ell]=a[ell]*vg[j,ell]+b[ell]*(gg+yg[j])
            yf[j]=ff;yg[j]=gg
            fn[j]=ff;gn[j]=gg
            mf[j]+=mean_rate*(ff-mf[j]);mg[j]+=mean_rate*(gg-mg[j])
        if beta!=0:
            nonlinear_correction(fn,mf,rad,1.,beta,1/C,0.,q,tmp,sl,flux,rhs)
            nonlinear_correction(gn,mg,rad,-1.,beta,1/C,0.,q,tmp,sl,flux,rhs)
        if viscosity>0:
            volume_loss(fn,mf,rad,viscosity*dt/dx**2,q)
            volume_loss(gn,mg,rad,viscosity*dt/dx**2,q)
        f,fn=fn,f;g,gn=gn,g
        if it%128==0:
            for j in range(cells):
                mach=max(abs(f[j]-mf[j]),abs(g[j]-mg[j]))/(RHO*C*C)
                maxmach=max(maxmach,mach)
            if maxmach*abs(beta)>.35:raise ValueError('Nonlinear CFL/weak-wave bound exceeded')
        output[it,0]=p;output[it,1]=u;output[it,2]=xx;output[it,3]=yy
        output[it,4]=ua;output[it,5]=ul;output[it,6]=pm;output[it,7]=plip
        output[it,8]=incoming;output[it,9]=outgoing;output[it,10]=bellflow
        output[it,11]=bellp
    return output,maxmach


def simulate_coupled(cfg,seconds,flip,pressure,yeq,beta=1.2,hold=None,source=None,
                     feedback=True,matched_end=False,mean_hz=20.):
    if hold is None:hold=seconds-.14
    if source is None:source=np.zeros(0)
    return run_core(cfg['radius'],cfg['reflection'],cfg['weights'],cfg['a'],cfg['b'],
        cfg['loss_k'],cfg['loss_den'],cfg['fs'],seconds,flip,pressure,yeq,beta,
        hold,cfg['viscosity'],source,feedback,matched_end,mean_hz)


def radiated(trace,fs):
    pressure=RHO/(4*np.pi)*np.gradient(trace[:,10])*fs
    return signal.resample_poly(pressure,48000,fs)


def impedance(cfg,freq):
    """Analytical linear limit of the discrete lossy waveguide, for validation."""
    omega=2*np.pi*np.asarray(freq);zz=np.exp(-1j*omega/cfg['fs'])
    s=2*cfg['fs']*(1-zz)/(1+zz)
    frac=np.sum(cfg['weights'][:,None]*s[None,:]/(s[None,:]+cfg['theta'][:,None]),axis=0)
    radius=cfg['radius'];zbell=cfg['zbell'];rr=.61**2/.25
    inert=s*.61*radius[-1]/C
    z=zbell*rr*inert/(rr+inert)
    for j in range(len(radius)-1,-1,-1):
        h=1/(1+cfg['loss_k'][j]*frac)
        prop=1j*omega/cfg['fs']-np.log(h)
        zc=RHO*C/(np.pi*radius[j]**2)
        z=(z*np.cosh(prop)+zc*np.sinh(prop))/(z/zc*np.sinh(prop)+np.cosh(prop))
    return z


if __name__=='__main__':
    import time
    cfg=setup();print('cells',len(cfg['radius']),'loss fit',cfg['loss_fit_max_relative_error'],flush=True)
    settings=json.loads((ROOT/'data/performance.json').read_text())
    s=next(s for s in settings if s['note']=='Bb4')
    for beta in [0.,1.2]:
        t=time.time()
        tr,mach=simulate_coupled(cfg,.45,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'],beta=beta,hold=.45)
        print('beta',beta,'elapsed',time.time()-t,'mach',mach,measure(tr,fs=cfg['fs'],start=.25),flush=True)
        np.savez_compressed(ROOT/'tmp'/f'coupled_probe_{beta}.npz',trace=tr,fs=cfg['fs'])
