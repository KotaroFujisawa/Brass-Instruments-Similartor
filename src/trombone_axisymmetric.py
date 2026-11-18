"""Conn 8H tenor trombone on the axisymmetric Euler/lip solver."""
import argparse
import hashlib
import os
import time
from functools import cache

from axisymmetric_cuda import Device, Run
from axisymmetric_mesh import annular_mesh
from render_axisymmetric import radiation, TRACE_FS
from trumpet import C, RHO, ROOT, json, np, signal
from scipy.io import wavfile

OUT = ROOT/'output/trombone_axisymmetric'
PRESET = ROOT/'data/performance_trombone.json'
BORE = ROOT/'data/conn8h_fig4_22.csv'
MOUTHPIECE = ROOT/'data/wick5al_calibrated.csv'
LIPS = dict(width=.0105, thickness=.003, joint=.006, rest_x=.0015, mass_product=3.25,
            cup_area=np.pi*.012875**2, quality=3.)
SCALE = [('Bb2',46,1),('C3',48,6),('D3',50,4),('Eb3',51,3),('F3',53,1),
         ('G3',55,4),('A3',57,2),('Bb3',58,1),('C4',60,3),('D4',62,1),
         ('Eb4',63,3),('F4',65,1),('G4',67,4),('A4',69,2),('Bb4',70,1)]


def slide_travel(position):
    if position not in range(1,8):
        raise ValueError('Slide position must be from 1 through 7.')
    bb1 = 440*2**((34-69)/12)
    return C/(4*bb1)*(2**((position-1)/12)-1)


@cache
def profile(position, trim=0.):
    mouthpiece = np.loadtxt(MOUTHPIECE, delimiter=',', skiprows=1)
    bore = np.loadtxt(BORE, delimiter=',', skiprows=1)
    # The outer slide moves at both ends of the U. Its overlapping bore, rather
    # than the bell taper, is lengthened. A small separate correction in the
    # parallel legs represents main tuning and precise slide adjustment.
    shift = slide_travel(position)
    z = bore[:,0].copy()
    z[z>=.946] += 2*shift
    z += trim*np.clip((bore[:,0]-.260)/(.756-.260),0.,1.)
    return np.r_[mouthpiece[:,0],z[1:]], np.r_[mouthpiece[:,1],bore[1:,1]]


@cache
def mesh(position, trim=0., dz=.006, dr=.0045):
    x, r = profile(position,trim)
    shift=slide_travel(position)
    anchors = [0., .080, .256, .756+trim, .766+trim, .946+2*shift+trim,
               .956+2*shift+trim, x[-1]-.105, x[-1]-.022, x[-1]]
    # The straight pipe can be discretized more coarsely than the rim and cup.
    spacing = [min(.003,dz),dz,dz,dz,dz,dz,dz,min(.004,dz),min(.003,dz)]
    z = np.concatenate([np.linspace(a,b,int(np.ceil((b-a)/h))+1)[:-1]
                        for a,b,h in zip(anchors[:-1],anchors[1:],spacing)] + [x[-1:]])
    return annular_mesh(z,np.interp(z,x,r),dr)


def initial_settings():
    notes=[]
    for name,midi,pos in SCALE:
        f=440*2**((midi-69)/12)
        notes.append(dict(note=name,midi=midi,position=pos,target_hz=f,
                          lip_frequency_hz=f*(1.135 if midi<=51 else 1.095),
                          blowing_pressure_pa=3600+(midi-46)*105,
                          rest_y_m=0. if midi>=58 else .000065,trim_m=0.))
    return notes


def synthesize(device, note, duration, *, hold=None, feedback=True, dz=.006, dr=.0045, snapshot=False):
    grid=mesh(note['position'],note.get('trim_m',0.),dz,dr)
    lips=LIPS | note.get('lip_parameters',{})
    with Run(device,grid,duration,lip_hz=note['lip_frequency_hz'],pressure=note['blowing_pressure_pa'],
             rest_y=note['rest_y_m'],hold=hold,feedback=feedback,attack=note.get('attack_s',.035),
             attack_boost=note.get('attack_boost',0.),attack_decay=note.get('attack_decay_s',.024),
             lip_parameters=lips) as simulation:
        if snapshot:
            simulation.until(.42)
            snapshot_field=simulation.field()
        simulation.until()
        trace=simulation.history()
        details=dict(dt_s=simulation.dt,steps=simulation.steps,cells=len(grid['cells']),
                     rings=[int(grid['rings'].min()),int(grid['rings'].max())],
                     fallback_count=int(simulation.fallbacks),tube_length_m=grid['z'][-1])
    if snapshot:
        np.savez_compressed(OUT/f'{note["note"]}_field.npz',state=snapshot_field,
                            z=grid['z'],radius=grid['radius'],rings=grid['rings'],
                            offsets=grid['offsets'],at_s=.42)
    return trace,details


def period(waveform, fs, start, end):
    """Find a periodicity without telling the estimator the desired note."""
    y=waveform[round(start*fs):round(end*fs)].astype(float)
    y-=y.mean()
    rms=float(np.std(y))
    if rms<1e-8:
        return 0.,0.,rms
    ac=signal.fftconvolve(y,y[::-1])[len(y)-1:]
    ac/=np.maximum(1,len(y)-np.arange(len(y)))
    lo,hi=int(fs/950),int(fs/70)
    peaks=signal.find_peaks(ac[lo:hi+1])[0]+lo
    if not len(peaks):
        return 0.,0.,rms
    # Compare the actual interpolated peak heights before selecting the first
    # period. This avoids octave errors when a narrow peak lies between samples.
    delta=.5*(ac[peaks-1]-ac[peaks+1])/(ac[peaks-1]-2*ac[peaks]+ac[peaks+1])
    heights=ac[peaks]-.25*(ac[peaks-1]-ac[peaks+1])*delta
    best=np.flatnonzero(heights>=max(heights)*.99)[0]
    return float(fs/(peaks[best]+delta[best])),float(heights[best]/ac[0]),rms


def analyze(note,trace,start,end):
    y=radiation(trace)
    pitch,coherence,rms=period(y,48000,start,end)
    window=trace[round(start*TRACE_FS):round(end*TRACE_FS)].astype(float)
    cup_area = (LIPS | note.get('lip_parameters',{}))['cup_area']
    lip_p=window[:,0]-window[:,8]-RHO*C/cup_area*window[:,1]
    return dict(note=note['note'],target_hz=note['target_hz'],frequency_hz=pitch,
                cents=float(1200*np.log2(pitch/note['target_hz'])) if pitch else None,
                periodicity=coherence,window_s=[start,end],radiated_rms_pa=rms,
                mouthpiece_rms_pa=float(np.std(window[:,0])),
                incoming_characteristic_rms_pa=float(np.std(window[:,8])),
                boundary_residual_pa=float(max(abs(lip_p))))


def run(probe=None,only=None):
    os.makedirs(OUT,exist_ok=True)
    device=Device()
    notes=json.loads(PRESET.read_text()) if PRESET.exists() else initial_settings()
    simulation_files = [ROOT/'src'/name for name in ('axisymmetric_euler.cu','axisymmetric_cuda.py',
                                                     'axisymmetric_mesh.py','trombone_axisymmetric.py')]
    signature=hashlib.sha256(b''.join(path.read_bytes() for path in [*simulation_files,BORE,MOUTHPIECE])).digest()
    report=[];waves=[];offset=0.
    for note in notes:
        name=note['note']
        if only and name not in only:
            continue
        duration=probe or (1.1 if name=='Bb4' else .72)
        end=duration if probe else duration-.16
        key=hashlib.sha256(signature+json.dumps([note,duration,bool(probe)],sort_keys=True).encode()).hexdigest()
        path=OUT/f'{name}_{"probe" if probe else "coupled"}.npz'
        begin=time.perf_counter()
        cached=False
        if path.exists():
            with np.load(path) as saved:
                if str(saved['key'])==key:
                    tr=saved['trace'];details=json.loads(str(saved['details']));cached=True
        if not cached:
            tr,details=synthesize(device,note,duration,hold=duration if probe else None,
                                  snapshot=probe is None and name=='Bb3')
            details['compute_s']=time.perf_counter()-begin
            np.savez(path,key=key,trace=tr,details=json.dumps(details),fs=TRACE_FS)
        row=analyze(note,tr,max(.14,end-.22),end)
        row.update(start_s=offset,duration_s=duration,position=note['position'],preset=note,**details)
        report.append(row);print(json.dumps(row),flush=True)
        (OUT/('probe.json' if probe else 'render_status.json')).write_text(json.dumps(report,indent=2))
        if not probe:
            waves.extend([radiation(tr),np.zeros(2880)])
        offset+=duration+.06
    if probe is None and not only:
        audio=np.concatenate(waves)
        if not np.isfinite(audio).all() or len(report)!=len(SCALE):
            raise RuntimeError('The scale render is incomplete or contains nonfinite samples.')
        gain=.89/np.max(abs(audio))
        pcm=np.round(32767*gain*audio).astype('i2')
        wavfile.write(OUT/'trombone_Bb2_Bb4_major_f_axisymmetric.wav',48000,pcm)
        result=dict(model='Conn 8H + 5AL / axisymmetric Euler + two-dimensional lips',
                    notes=report,rate_hz=48000,duration_s=len(audio)/48000,gain_fs_per_pa=float(gain),
                    clipped_samples=int(np.count_nonzero(abs(pcm.astype('i4'))>=32767)),gpu=device.name)
        (OUT/'measurements.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',type=float);parser.add_argument('--only',nargs='+')
    args=parser.parse_args();run(args.probe,args.only)
