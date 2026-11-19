"""Trumpet: Adachi & Sato (1996) two-dimensional lips + geometry-derived bore.

SI units throughout. The acoustic resonator is a lossy cylindrical TMM chain,
including the cup, throat, backbore, valves, leadpipe and bell. No sample playback,
oscillator bank, pitch shifter or imposed audio-frequency lip excitation is used.
"""
from __future__ import annotations
import sys
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.deps'))
if os.environ.get('BRASS_DEPS'):
    sys.path.insert(0,os.environ['BRASS_DEPS'])
import json
import math
import numpy as np
from scipy import signal, optimize
from scipy.interpolate import PchipInterpolator
from numba import njit

RHO = 1.2
C = 343.0
SCUP = 2.3e-4
B = 7e-3
D = 2e-3
FS = 96000
NFFT = 32768


def body_profile():
    return np.loadtxt(ROOT / 'data/york_trumpet_fig7_52.csv', delimiter=',', skiprows=1)


def mouthpiece(depth=0.018):
    """Reconstructed cup; assumed throat/backbore, calibrated to book's 806 Hz.

    This is NOT a measured Bach 7C profile. Its calibration target is the 7C
    resonance in The Science of Brass Instruments, Fig. 4.40, printed p.147.
    """
    x = np.linspace(0, depth, 41)
    r0 = np.sqrt(SCUP/np.pi)
    throat = 0.0037/2
    r = throat + (r0-throat)*np.sqrt(np.maximum(0, 1-(x/depth)**2))
    x = np.r_[x, depth+0.006, depth+0.061]
    r = np.r_[r, throat, 0.0046]
    return x, r


def geometry(valves=(), depth=0.018, tuning=0.0, dx=0.002):
    body = body_profile()
    xb, rb = body[:,0].copy(), body[:,1]/2
    # Nominal length from book Table 7.1: c/(2*116.54 Hz).
    nominal = C/(2*116.54094038)
    extension = sum(nominal*(2**(s/12)-1) for s in valves)+tuning
    # Add straight valve tubing within the nearly cylindrical section.
    split = 0.55
    bi = np.searchsorted(xb, split)
    rs = np.interp(split, xb, rb)
    if extension < 0:
        # Remove the interval [split, right] from the cylindrical tuning section.
        right = split-extension
        if right > .65:
            raise ValueError('Shortening extends beyond the cylindrical tuning section.')
        tail = xb > right
        xb = np.r_[xb[:bi], split, xb[tail]+extension]
        rb = np.r_[rb[:bi], rs, rb[tail]]
    else:
        xb = np.r_[xb[:bi], split, split+extension, xb[bi:]+extension]
        rb = np.r_[rb[:bi], rs, rs, rb[bi:]]
    xm,rm=mouthpiece(depth)
    # Mouthpiece exit joins the digitised mouthpipe entrance. No cup double-count.
    x=np.r_[xm, xm[-1]+xb[1:]]
    r=np.r_[rm,rb[1:]]
    edges=np.linspace(0,x[-1],int(np.ceil(x[-1]/dx))+1)
    centers=(edges[1:]+edges[:-1])/2
    radius=np.interp(centers,x,r)
    return np.diff(edges),radius, x,r


def tmm(freq, lengths, radius):
    """Book Eqs. 4.108, 4.109 (leading causal term), 4.111, 4.112.

    Radiation is a passive first-order RL approximation matching Eq.4.121
    at low ka; its high-frequency asymptote is only approximate.
    """
    f=np.maximum(np.asarray(freq,float), 0.05)
    w=2*np.pi*f
    zbell=RHO*C/(np.pi*radius[-1]**2)
    # Z/Zc = r*(j*k*0.61*a)/(r+j*k*0.61*a), r=0.61^2/0.25.
    # Its LF expansion is (ka)^2/4 + j*0.61*ka.
    inert=1j*w*0.61*radius[-1]/C
    resistance=0.61**2/0.25
    p=zbell*resistance*inert/(resistance+inert)
    u=np.ones(f.shape,complex)
    for length,a in zip(lengths[::-1],radius[::-1]):
        rv=632.8*np.sqrt(f)*a
        # Equal real/imag sqrt(f) coefficients ensure causal wall-loss term.
        g=1j*w/C+(1+1j)*w/C*1.045/rv
        zc=RHO*C/(np.pi*a*a)*(1+(1-1j)*0.369/rv)
        co=np.cosh(g*length); si=np.sinh(g*length)
        p,u=co*p+zc*si*u, si/zc*p+co*u
    return p/u, 1/u


def calibrate_mouthpiece():
    def peak(depth):
        x,r=mouthpiece(depth)
        e=np.linspace(0,x[-1],200)
        a=np.interp((e[1:]+e[:-1])/2,x,r)
        def loss(f):return -abs(tmm(np.array([f]),np.diff(e),a)[0][0])
        return optimize.minimize_scalar(loss,bounds=(450,1400),method='bounded').x
    depth=optimize.brentq(lambda dep:peak(dep)-806, .004,.035,xtol=1e-9)
    return depth,peak(depth)


def resonator(valves=(), depth=.018, tuning=0., fs=FS, nfft=NFFT):
    lengths,rad,x,r=geometry(valves,depth,tuning)
    freq=np.fft.rfftfreq(nfft,1/fs)
    zin,transfer=tmm(freq,lengths,rad)
    z0=RHO*C/SCUP
    reflection=(zin-z0)/(zin+z0)
    reflection[0]=-1
    reflection[-1]=reflection[-1].real
    raw=np.fft.irfft(reflection,nfft)
    # Adachi-Sato Sec.I.B: twice the even component for t>0.
    # Retain a 100 ms causal response with a short smooth tail.
    n=min(int(.1*fs),nfft//2)
    impulse=raw[:n].copy()
    impulse[1:]+=raw[-1:-n:-1]
    impulse[-int(.01*fs):]*=np.cos(np.linspace(0,np.pi/2,int(.01*fs)))**2
    # Exact zero-frequency pressure release. A small DC residual would otherwise
    # generate an unphysical steady mouthpiece pressure.
    dc_adjust=-1-impulse.sum()
    weights=np.hanning(n); weights[:int(.015*fs)]=0
    impulse+=dc_adjust*weights/weights.sum()
    actual=np.fft.rfft(impulse,nfft)
    peaks,_=signal.find_peaks(abs(zin),prominence=z0)
    peak_f=freq[peaks][freq[peaks]<2500]
    return dict(r=impulse,z0=z0,freq=freq,zin=zin,transfer=transfer,
                lengths=lengths,radius=rad,x=x,radius_points=r,
                reflection=reflection,actual_reflection=actual,peaks=peak_f,
                acausal_energy_ratio=float(np.sum(raw[nfft//2:]**2)/np.sum(raw**2)),
                dc_adjust=float(dc_adjust),fs=fs)


@njit(cache=True)
def simulate(r,z0,fs,duration,flip,pressure,yeq,attack=.045,release=.065,hold=-1.):
    """Equations (1)-(8) of Adachi-Sato, coupled implicitly to Eq.(9).

    Damped centered lip difference; backward Euler air-column inertia and exact
    quadratic flow solve. Signed Q|Q| extension permits brief flow reversal.
    No displacement/pressure clipping: bad trajectories raise instead.
    """
    count=int(duration*fs)
    if hold<0:hold=duration-release-.06
    out=np.zeros((count,8))
    hist=np.zeros(len(r)*2)
    k=1.5*flip
    m=1.5/(4*np.pi*np.pi*flip)
    omega=2*np.pi*flip
    dt=1/fs
    xx=1e-3; yy=yeq
    xp=xx;yp=yy
    p=0.;plip=0.;ua=0.
    cursor=0
    # r[0] is inside the simultaneous nonlinear boundary solve.
    zeff=z0*(1+r[0])/(1-r[0])
    for n in range(count):
        t=n/fs
        env=1.
        if t<attack:env=.5-.5*np.cos(np.pi*t/attack)
        elif t>hold:
            env=.5+.5*np.cos(np.pi*min(1.,(t-hold)/release))
        pm=pressure*env
        q=3.0 if yy>0 else .5
        damp=omega/q
        fx=2*B/m*(pm-p)*(0.004-yy)
        fy=2*B/m*((pm-p)*xx+D*plip)
        if yy<0:fy-=6*k/m*yy
        xn=(2*xx-(1-damp*dt/2)*xp+dt*dt*(fx-omega*omega*(xx-.001)))/(1+damp*dt/2)
        yn=(2*yy-(1-damp*dt/2)*yp+dt*dt*(fy-omega*omega*(yy-yeq)))/(1+damp*dt/2)
        vx=(xn-xp)/(2*dt);vy=(yn-yp)/(2*dt)
        ul=B*(xx*vy-(yy-.004)*vx)
        xp,yp=xx,yy
        xx,yy=xn,yn
        area=max(2*B*yy,0.)
        # Recent outgoing p+Z0*U in a duplicated ring, contiguous for dot product.
        memory=0.
        for j in range(1,len(r)):
            memory+=r[j]*hist[cursor+j]
        h=memory/(1-r[0])
        if area>1e-12:
            # Sum of Eqs.7/8: A Ua|Ua| + (rho*d/S)/dt * Ua = pm-p + inert*Ua_old.
            aa=RHO*(.5/(area*area)-1/(SCUP*area)+1/(SCUP*SCUP))
            inert=RHO*D/(area*dt)
            bb=inert+zeff
            cc=pm-h-zeff*ul+inert*ua
            ua=2*cc/(bb+np.sqrt(bb*bb+4*aa*abs(cc)))
            p=zeff*(ua+ul)+h
            plip=p-RHO*ua*abs(ua)*(1/(SCUP*area)-1/(SCUP*SCUP))
        else:
            ua=0.
            p=zeff*ul+h
            plip=p
        u=ua+ul
        if not np.isfinite(p) or abs(p)>1e6 or abs(xx)>.05 or abs(yy)>.02:
            raise ValueError('Unstable/nonphysical lip trajectory; adjust integration or embouchure')
        hist[cursor]=p+z0*u
        hist[cursor+len(r)]=hist[cursor]
        cursor-=1
        if cursor<0:cursor=len(r)-1
        out[n,0]=p;out[n,1]=u;out[n,2]=xx;out[n,3]=yy
        out[n,4]=ua;out[n,5]=ul;out[n,6]=pm;out[n,7]=plip
    return out


def measure(trace,fs=FS,start=.18,end=None,*,min_hz=150,max_hz=1200):
    section=trace[int(start*fs):int(end*fs) if end else len(trace)]
    y=section[:,0]
    y=y-y.mean()
    if len(y)<64:return dict(frequency=0.,rms=0.,periodicity=0.)
    ac=signal.fftconvolve(y,y[::-1],mode='full')[len(y)-1:]
    ac/=np.maximum(1,np.arange(len(y),0,-1))
    lo=int(fs/max_hz);hi=min(int(fs/min_hz),len(ac)-2)
    peaks,_=signal.find_peaks(ac[lo:hi])
    if not len(peaks) or ac[0]<1e-9:return dict(frequency=0.,rms=float(np.std(y)),periodicity=0.)
    peaks=peaks+lo
    best=peaks[np.argmax(ac[peaks])]
    # First peak within 99% of maximum avoids octave/subharmonic selection.
    good=peaks[ac[peaks]>.99*ac[best]]
    best=good[0]
    delta=.5*(ac[best-1]-ac[best+1])/(ac[best-1]-2*ac[best]+ac[best+1])
    return dict(frequency=float(fs/(best+delta)),rms=float(np.std(y)),
                periodicity=float(ac[best]/ac[0]),
                closure_fraction=float(np.mean(section[:,3]<=0)))


def radiate(trace,res):
    """Compact-monopole approximation at 1 m, linear bore transfer.

    Ubell/Uin comes from the same full mouthpiece/bore TMM, avoiding the common
    mistake of delivering the internal mouthpiece pressure as the audible sound.
    """
    n=signal.fftconvolve(trace[:,1],np.fft.irfft(res['transfer'],NFFT)[:int(.12*res['fs'])])
    # rho/(4*pi*r)*dUbell/dt, free-space compact radiator. No room/reverb.
    y=RHO/(4*np.pi)*np.gradient(n)*res['fs']
    return y[:len(trace)]


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',action='store_true')
    args=parser.parse_args()
    dep,fr=calibrate_mouthpiece()
    print('Mouthpiece depth',dep,'resonance',fr,flush=True)
    res=resonator(depth=dep)
    print('Peaks',res['peaks'],'acausal',res['acausal_energy_ratio'],'dc',res['dc_adjust'],flush=True)
    np.savez(ROOT/'tmp/resonator_open.npz',**res)
    (ROOT/'data/model.json').write_text(json.dumps(dict(cup_depth_m=dep,isolated_mouthpiece_hz=fr),indent=2))
    for flip in [250,490,680,960]:
        for yeq in [0.,.00002,.00006,.00012]:
            try:
                tr=simulate(res['r'],res['z0'],FS,.42,flip,5000,yeq,hold=.42)
                print(flip,yeq,measure(tr,start=.24),flush=True)
            except ValueError as e: print(flip,yeq,str(e),flush=True)
