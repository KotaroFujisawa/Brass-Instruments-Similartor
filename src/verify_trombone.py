"""Check geometry, small-signal impedance and the coupled trombone."""
import argparse
import os

from trombone_axisymmetric import (OUT, MOUTHPIECE, PRESET, TRACE_FS, Device, Run, analyze,
                                   mesh, profile, slide_travel, synthesize)
from axisymmetric_mesh import annular_mesh
from trumpet import C, RHO, json, np
from scipy.special import jn_zeros
import verify_axisymmetric as fluid_tests


def response(device, grid, duration):
    with Run(device,grid,duration,mode=3,pressure=100.) as sim:
        sim.until()
        trace=sim.history()
    t=(np.arange(len(trace))+1)/TRACE_FS
    area=np.pi*grid['radius'][0]**2
    flow=area*(100/(RHO*C))*np.exp(-.5*((t-.00045)/.00007)**2)
    pressure=trace[:,8]+RHO*C/area*flow
    n=131072
    f=np.fft.rfftfreq(n,1/TRACE_FS)
    z=np.fft.rfft(pressure,n)/(np.fft.rfft(flow,n)+1e-30)
    valid=f<1200
    return f[valid],abs(z[valid])


def peak_in(f,z,low,high):
    indices=np.flatnonzero((f>=low)&(f<=high))
    i=indices[np.argmax(z[indices])]
    y=np.log(z[i-1:i+2]);frac=.5*(y[0]-y[2])/(y[0]-2*y[1]+y[2])
    return float(f[i]+frac*(f[1]-f[0]))


def geometry_and_solver(device):
    os.makedirs(OUT,exist_ok=True)
    fluid_tests.OUT=OUT
    fluid_tests.solver(device)
    xyz=np.loadtxt(MOUTHPIECE,delimiter=',',skiprows=1)
    z,r=xyz.T
    volume=np.pi*np.sum(np.diff(z)*(r[:-1]**2+r[1:]**2+r[:-1]*r[1:])/3)*1e6
    zs=np.linspace(0,z[-1],round(z[-1]/.0024)+1)
    detached=annular_mesh(zs,np.interp(zs,z,r),dr=.0045)
    fm,zm=response(device,detached,.13)
    mp=peak_in(fm,zm,360,700)
    f,z1=response(device,mesh(1),.31)
    _,z7=response(device,mesh(7),.31)
    peak_windows=[(25,55),(96,128),(151,209),(215,249),(269,305),(319,369)]
    peaks=[peak_in(f,z1,*band) for band in peak_windows]
    book=[38,111,172,232,292,347]
    longitudinal_length=float(profile(1)[0][-1]);seven=float(profile(7)[0][-1])
    result=dict(detached_volume_ml=float(volume),book_volume_ml=12.,detached_euler_resonance_hz=mp,
                book_detached_impedance_peak_hz=535.,first_position_euler_peaks_hz=peaks,
                book_first_position_peaks_approx_hz=book,first_position_length_m=longitudinal_length,
                seventh_position_length_m=seven,seventh_position_travel_m=slide_travel(7),
                bell_first_radial_cutoff_ideal_cylinder_hz=float(C*jn_zeros(1,1)[0]/(2*np.pi*.1065)))
    result['passed']=(abs(volume-12)<.02 and abs(mp-535)<90
                      and max(abs(np.array(peaks)-book))<19
                      and abs(seven-longitudinal_length-2*slide_travel(7))<1e-10)
    (OUT/'geometry_verification.json').write_text(json.dumps(result,indent=2))
    np.savez(OUT/'impedance_curves.npz',f_hz=f,pos1=z1,pos7=z7,mp_f_hz=fm,mp_z=zm)
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Geometry/impedance comparison is outside the verification range.')


def coupling(device):
    notes=json.loads(PRESET.read_text());s=next(n for n in notes if n['note']=='Bb3')
    with np.load(OUT/'Bb3_coupled.npz') as saved:
        trace=saved['trace']
    ablation,_=synthesize(device,s,.58,hold=.58,feedback=False)
    fine,_=synthesize(device,s,.58,hold=.58,dz=.0045,dr=.0035)
    a=analyze(s,trace,.34,.56);b=analyze(s,ablation,.34,.56);c=analyze(s,fine,.34,.56)
    cents=1200*np.log2(c['frequency_hz']/a['frequency_hz']) if min(a['frequency_hz'],c['frequency_hz'])>0 else None
    report=dict(coupled=a,without_return_to_lips=b,refined_grid=c,refined_change_cents=cents,
        passed=a['boundary_residual_pa']<.02 and a['mouthpiece_rms_pa']>500 and
               a['radiated_rms_pa']>10*max(b['radiated_rms_pa'],1e-9))
    (OUT/'coupling_tests.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    if not report['passed']:raise RuntimeError('Closed-loop lip coupling check failed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--coupling',action='store_true')
    args=parser.parse_args();device=Device()
    if args.coupling:coupling(device)
    else:geometry_and_solver(device)
