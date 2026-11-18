"""Check cylindrical geometry, a radial eigenmode, and shock conservation."""
import argparse
import os

from axisymmetric_cuda import Device, Run
from axisymmetric_mesh import annular_mesh, trumpet_mesh
from render_axisymmetric import OUT, analyze, settings, synthesize
from scipy.special import j0, jn_zeros
from trumpet import C, RHO, json, np, signal

P0=1/1.4


def pressure(q):
    return .4*(q[:,3]-.5*np.sum(q[:,1:3]**2,axis=1)/(1+q[:,0]))


def exact_sod(s):
    """Density and pressure in the Sod Riemann solution, x/t in c0 units."""
    g=1.4;left_sound=np.sqrt(g);star_p=.303130178;star_u=.92745262
    left_d=star_p**(1/g);right_d=.26557371
    tail=star_u-np.sqrt(g*star_p/left_d)
    shoulder=1.75215573
    fan=(2*left_sound-(g-1)*s)/(g+1)
    d=np.select([s<-left_sound,s<tail,s<star_u,s<shoulder],
                [1.,np.maximum(0.,fan/left_sound)**(2/(g-1)),left_d,right_d],default=.125)
    p=np.select([s<-left_sound,s<tail,s<shoulder],
                [1.,np.maximum(0.,fan/left_sound)**(2*g/(g-1)),star_p],default=.1)
    return d,p


def shock(device, nx):
    z=np.linspace(0,1,nx+1);mesh=annular_mesh(z,np.full_like(z,.03),dr=.006,loss=False)
    x=mesh['centers'][:,0];before=np.zeros((len(x),4),'f4')
    before[:,0]=np.where(x<.5,1.,.125)-1
    before[:,3]=(np.where(x<.5,1.,.1)-P0)/.4
    tau=.2
    with Run(device,mesh,tau/C,initial=before,mode=2,guard=2.45,trace_fs=960000) as sim:
        sim.until();after=sim.field();elapsed=sim.steps*sim.dt
        fallback=sim.fallbacks
    exact_d,exact_p=exact_sod((x-.5)/(C*elapsed));p=P0+pressure(after)
    vol=mesh['volumes'];diff=(after-before).astype('f8')
    xbar=mesh['centers'][::mesh['rings'][0],0]
    pbar=p[::mesh['rings'][0]];mask=(xbar>.78)&(xbar<.92)
    ids=np.flatnonzero(mask);jump=ids[np.argmin(np.diff(pbar)[ids])]
    shock_location=.5*(xbar[jump]+xbar[jump+1]);expected=.5+1.75215573*C*elapsed
    result=dict(nx=nx,pressure_mae=float(np.average(abs(p-exact_p),weights=vol)),
        density_mae=float(np.average(abs(1+after[:,0]-exact_d),weights=vol)),
        normalized_mass_change=float(np.sum(vol*diff[:,0])/vol.sum()),
        normalized_energy_change=float(np.sum(vol*diff[:,3])/vol.sum()),
        shock_location_m=float(shock_location),exact_shock_location_m=float(expected),
        max_radial_momentum=float(max(abs(after[:,2]))),fallbacks=int(fallback))
    return result,xbar,pbar,x,exact_p


def radial_mode(device):
    z=np.linspace(0,.032,9);r=.053
    mesh=annular_mesh(z,np.full_like(z,r),dr=r/32,loss=False)
    radius=mesh['centers'][:,1]
    k=jn_zeros(1,1)[0]/r;target=C*k/(2*np.pi);period=1/target
    profile=j0(k*radius);vol=mesh['volumes']
    profile-=np.average(profile,weights=vol)
    before=np.zeros((len(radius),4),'f4')
    before[:,0]=1e-4*profile;before[:,3]=before[:,0]/.4
    ts=np.arange(121)*period/20
    with Run(device,mesh,ts[-1]+1e-5,initial=before,mode=1,trace_fs=960000) as sim:
        amplitude=[]
        for t in ts:
            sim.until(t);amplitude.append(np.average(pressure(sim.field())*profile,weights=vol))
    amplitude=np.array(amplitude)/amplitude[0]
    peaks=signal.find_peaks(amplitude)[0]
    left,mid,right=amplitude[peaks-1],amplitude[peaks],amplitude[peaks+1]
    offset=.5*(left-right)/(left-2*mid+right)
    fitted=float(1/np.mean(np.diff(ts[peaks]+offset*period/20)))
    return dict(radius_m=r,rings=32,frequency_hz=fitted,analytic_frequency_hz=float(target),
                error_percent=float(100*(fitted/target-1)),last_peak_amplitude=float(amplitude[peaks[-1]])),ts,amplitude


def solver(device):
    os.makedirs(OUT,exist_ok=True)
    mesh=trumpet_mesh();q=np.zeros((len(mesh['cells']),4),'f4');q[:,0]=.003;q[:,3]=.003/.4
    with Run(device,mesh,.008,initial=q,mode=1) as sim:
        sim.until();after=sim.field()
    rest=float(np.max(abs(q-after)))
    coarse,x0,p0,_,_=shock(device,200)
    fine,x1,p1,xe,pe=shock(device,400)
    radial,t,a=radial_mode(device)
    passed=(rest<1e-6 and fine['pressure_mae']<.009 and fine['pressure_mae']<coarse['pressure_mae']
            and abs(fine['shock_location_m']-fine['exact_shock_location_m'])<.012
            and abs(fine['normalized_mass_change'])<1e-5 and abs(fine['normalized_energy_change'])<1e-5
            and abs(radial['error_percent'])<3.)
    report=dict(rest_max_state_error=rest,sod_coarse=coarse,sod_fine=fine,radial_mode=radial,passed=passed)
    (OUT/'solver_tests.json').write_text(json.dumps(report,indent=2))
    np.savez(OUT/'solver_curves.npz',x_coarse=x0,p_coarse=p0,x_fine=x1,p_fine=p1,
             x_exact=xe,p_exact=pe,radial_time_s=t,radial_amplitude=a,radial_f_ref=radial['analytic_frequency_hz'])
    print(json.dumps(report,indent=2),flush=True)
    if not passed:raise RuntimeError('Solver checks exceeded tolerance')


def coupling(device):
    s=next(s for s in settings() if s['note']=='Bb4')
    data=np.load(OUT/'Bb4_coupled.npz');trace=data['trace']
    no_feedback,_=synthesize(device,s,.45,hold=.45,feedback=False)
    fine,_=synthesize(device,s,.45,hold=.45,dz=.0026,dr=.0022)
    start,end=.25,.44
    a=analyze(s,trace,start,end);b=analyze(s,no_feedback,start,end);c=analyze(s,fine,start,end)
    convergence=(1200*np.log2(c['audio_frequency_hz']/a['audio_frequency_hz'])
                 if a['audio_frequency_hz'] and c['audio_frequency_hz'] else None)
    report=dict(coupled=a,without_feedback=b,refined=c,refined_pitch_difference_cents=convergence,
                passed=a['boundary_residual_max_pa']<.02 and a['mouthpiece_rms_pa']>500
                and a['mouthpiece_rms_pa']>2*b['mouthpiece_rms_pa'])
    (OUT/'coupling_tests.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    if not report['passed']:raise RuntimeError('Lip coupling check failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--coupling',action='store_true')
    args=parser.parse_args();gpu=Device()
    if args.coupling:coupling(gpu)
    else:solver(gpu)
