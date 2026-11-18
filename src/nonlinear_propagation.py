"""Prescribed-incident-wave, variable-area generalised Burgers propagation.

The Science of Brass Instruments, sections 6.2.2-6.2.3, Eqs.6.30,6.33,6.34.
This is a one-way propagation approximation, NOT a bidirectionally coupled
nonlinear lip/bore solver. It does not modify the lip excitation by returning
the distorted reflected wave. The legacy linear solver remains available.
"""
from trumpet import *
from scipy.fft import rfft, irfft, rfftfreq, next_fast_len

GAMMA=1.4
BETA=(GAMMA+1)/2
ETA=1.846e-5                 # Pa s; same 300 K reference as book section 4.7
THERMAL_CONDUCTIVITY=.0262   # W/(m K), assumed representative dry air
CP=1005.                    # J/(kg K)
BULK_VISCOSITY=0.           # relaxation/bulk viscosity not modelled
SOUND_DIFFUSIVITY=(4*ETA/3+BULK_VISCOSITY+THERMAL_CONDUCTIVITY*(GAMMA-1)/CP)/RHO


@njit(cache=True)
def _slope(q):
    n=len(q);sl=np.empty(n)
    for i in range(n):
        dm=q[i]-q[(i-1)%n];dp=q[(i+1)%n]-q[i]
        if dm*dp<=0:sl[i]=0.
        else:sl[i]=np.sign(dm)*min(2*abs(dm),2*abs(dp),abs(dm+dp)/2)
    return sl


@njit(cache=True)
def _rhs(q,coefficient,dt):
    """MUSCL-MC reconstruction, local Lax-Friedrichs entropy flux."""
    n=len(q);sl=_slope(q);flux=np.empty(n);rhs=np.empty(n)
    for i in range(n):
        j=(i+1)%n
        left=q[i]+sl[i]/2;right=q[j]-sl[j]/2
        speed=coefficient*max(abs(left),abs(right))
        flux[i]=-.25*coefficient*(left*left+right*right)-.5*speed*(right-left)
    for i in range(n):rhs[i]=-(flux[i]-flux[(i-1)%n])/dt
    return rhs


@njit(cache=True)
def burgers_step(q,coefficient,distance,dt,cfl=.4):
    """SSP-RK2 spatial marching, with CFL-controlled nonlinear substeps.

    q_x + d/dtau(-coefficient*q^2/2)=0. Here distance plays the role of
    evolution time and the retarded-time sample interval is the mesh spacing.
    Shock dissipation is entropy-consistent numerical flux, not clipping.
    """
    if coefficient==0:return q.copy(),0
    remaining=distance;steps=0
    out=q.copy()
    while remaining>1e-14:
        speed=coefficient*np.max(np.abs(out))
        dx=min(remaining,cfl*dt/max(speed,1e-30))
        v=out+dx*_rhs(out,coefficient,dt)
        out=.5*(out+v+dx*_rhs(v,coefficient,dt))
        remaining-=dx;steps+=1
    return out,steps


def propagation_geometry(valves,depth):
    """Select the main bore beyond the mouthpiece from the same TMM mesh."""
    lengths,radii,x,r=geometry(tuple(valves),depth,dx=.002)
    edges=np.r_[0.,np.cumsum(lengths)]
    split=int(np.searchsorted(edges,depth+.061))
    # split is a node, so the source and body acoustics use the exact same cut.
    body_lengths=lengths[split:]
    body_radii=radii[split:]
    return dict(lengths=lengths,radii=radii,split=split,
                body_lengths=body_lengths,body_radii=body_radii,
                body_x=np.cumsum(body_lengths)-body_lengths/2,
                length=float(sum(body_lengths)),start_m=float(edges[split]),
                inlet_diameter=float(2*body_radii[0]),outlet_diameter=float(2*body_radii[-1]))


def incident_transfer(geo,fs=192000,nfft=65536):
    """p+ at the mouthpipe / lip total volume flow, from the loaded linear TMM.

    p+=(p+rho*c*u)/2; the source contains reflected-wave contributions from the
    original linear resonator. This prevents using standing mouthpiece pressure
    as if it were a pure travelling wave. DC mean flow is not acoustic pressure.
    """
    freq=rfftfreq(nfft,1/fs)
    _,hfull=tmm(freq,geo['lengths'],geo['radii'])
    zb,hbody=tmm(freq,geo['body_lengths'],geo['body_radii'])
    zlocal=RHO*C/(np.pi*(geo['inlet_diameter']/2)**2)
    transfer=.5*(zb+zlocal)*hfull/hbody
    # Separate acoustic flow from the player's slowly varying mean airflow.
    # A 30 Hz analog Butterworth high-pass has negligible gain error at Bb3;
    # unlike just zeroing the FFT's DC bin it remains DC-free after truncation.
    s=2j*np.pi*freq;wc=2*np.pi*30.
    transfer*=s*s/(s*s+np.sqrt(2)*wc*s+wc*wc)
    transfer[0]=0.
    transfer[-1]=transfer[-1].real
    impulse=irfft(transfer,nfft)
    return impulse[:int(.12*fs)]


def propagate(drive,fs,geo,beta=BETA,dx=.01,wall=True,bulk=True,record=False):
    """Strang split Eq.6.34 + fractional wall losses + volume diffusion.

    q=p*D(x)/D(0) removes geometric spreading exactly. Then
        q_x - beta/(rho*c^3)*D(0)/D(x)*q*q_tau = linear losses.
    This follows Eq.6.34 directly, avoiding the Dmin normalization ambiguity
    between the printed Eqs.6.35-6.37. Diameter is interpolated from the bore.
    Retarded time removes the ordinary flight time; it is restored at output.
    """
    pad=int(.06*fs)
    n=next_fast_len(len(drive)+2*pad)
    q=np.zeros(n);q[pad:pad+len(drive)]=drive
    freq=rfftfreq(n,1/fs);omega=2*np.pi*freq
    length=geo['length'];count=int(np.ceil(length/dx));step=length/count
    locations=(np.arange(count)+.5)*step
    radius=np.interp(locations,geo['body_x'],geo['body_radii'])
    d0=geo['inlet_diameter'];dend=geo['outlet_diameter']
    # Leading Keefe loss + causal dispersion from the existing linear model.
    wall_integral=0.;total_steps=0;snapshots=[]
    if beta==0:
        wall_integral=float(np.sum(step/radius))
        attenuation=(1+1j)*2*np.pi*1.045/(C*632.8)*np.sqrt(freq)*wall_integral if wall else 0.
        attenuation+=SOUND_DIFFUSIVITY*omega**2*length/(2*C**3) if bulk else 0.
        q=irfft(rfft(q)*np.exp(-attenuation),n)
    else:
        for j,a in enumerate(radius):
            alpha=(1+1j)*2*np.pi*1.045/(C*632.8*a)*np.sqrt(freq) if wall else 0.
            alpha+=SOUND_DIFFUSIVITY*omega**2/(2*C**3) if bulk else 0.
            half_filter=np.exp(-alpha*step/2)
            q=irfft(rfft(q)*half_filter,n)
            coeff=beta/(RHO*C**3)*d0/(2*a)
            q,steps=burgers_step(q,coeff,step,1/fs)
            total_steps+=steps
            q=irfft(rfft(q)*half_filter,n)
            if record and j in [0,count//3,2*count//3,count-1]:
                snapshots.append((float((j+1)*step),q[pad:pad+len(drive)].copy()*d0/(2*a)))
    # Both versions have the same exact length, spreading and ordinary delay.
    spectrum=rfft(q)*np.exp(-1j*omega*length/C)
    out=irfft(spectrum,n)[pad:pad+len(drive)]*d0/dend
    info=dict(distance_m=length,equivalent_nonlinear_length_m=float(np.sum(step*d0/(2*radius))),
              nonlinear_substeps=total_steps,spatial_step_m=step,beta=beta,
              peak_input_mach=float(np.max(abs(drive))/(RHO*C*C)),
              inlet_diameter_m=d0,outlet_diameter_m=dend)
    return out,info,snapshots


def radiate_incident(incident,fs,geo,distance=1.):
    """Passive bell load then the same compact free-space radiation as V1.

    The reflected component is included in terminal U, but NOT propagated back
    toward the source. Nonlinear distortion outside the bell is not included.
    """
    n=next_fast_len(len(incident)+int(.12*fs))
    f=rfftfreq(n,1/fs);w=2*np.pi*f
    a=geo['outlet_diameter']/2
    z0=RHO*C/(np.pi*a*a)
    inert=1j*w*.61*a/C;resistance=.61**2/.25
    zrad=z0*resistance*inert/(resistance+inert)
    reflection=(zrad-z0)/(zrad+z0)
    u=rfft(incident,n)*(1-reflection)/z0
    radiated=irfft(RHO/(4*np.pi*distance)*1j*w*u,n)
    return radiated[:len(incident)]
