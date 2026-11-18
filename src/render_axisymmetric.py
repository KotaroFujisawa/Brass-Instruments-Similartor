"""Render a concert Bb-major scale by two-way coupling the Euler grid and lips."""
import argparse
import hashlib
import os
import time

from axisymmetric_cuda import Device, Run
from axisymmetric_mesh import trumpet_mesh
from trumpet import C, RHO, ROOT, json, measure, np, signal
from render_bidirectional import measure_audio
from scipy.io import wavfile

OUT = ROOT/'output/axisymmetric'
PRESET = ROOT/'data/performance_axisymmetric.json'
TRACE_FS = 192000


def settings():
    preset = PRESET if PRESET.exists() else ROOT/'data/performance_bidirectional.json'
    notes = json.loads(preset.read_text())
    assert [s['midi'] for s in notes] == [58,60,62,63,65,67,69,70,72,74,75,77,79,81,82]
    return notes


def synthesize(device, s, duration, *, hold=None, dz=.0035, dr=.003, feedback=True, snapshot=False):
    mesh=trumpet_mesh(s['valves_semitones'],dz,dr,s.get('tuning_m',0.))
    with Run(device,mesh,duration,lip_hz=s['lip_frequency_hz'],pressure=s['blowing_pressure_pa'],
             rest_y=s['rest_y_m'],attack=s.get('attack_s',.045),hold=hold,feedback=feedback) as sim:
        if snapshot:
            sim.until(.41);field=sim.field()
        sim.until();trace=sim.history()
        diagnostics=dict(cells=len(mesh['cells']),steps=sim.steps,dt_s=sim.dt,
                         min_rings=int(mesh['rings'].min()),max_rings=int(mesh['rings'].max()),
                         reconstruction_fallbacks=int(sim.fallbacks),wall_loss_fit_error=sim.loss_error)
    if snapshot:
        diagnostics['snapshot_s']=.41
        np.savez_compressed(OUT/f'{s["note"]}_field.npz',field=field, z=mesh['z'],radius=mesh['radius'],
                            rings=mesh['rings'],offsets=mesh['offsets'])
    return trace,diagnostics


def radiation(trace):
    # Compact monopole at one metre: pressure from the total outward volume flow.
    output=RHO/(4*np.pi)*np.gradient(trace[:,9].astype('f8'))*TRACE_FS
    return signal.resample_poly(output,1,4)


def analyze(s, trace, start, end):
    yn=radiation(trace)
    audio=measure_audio(yn,start,end)
    steady=trace[round(start*TRACE_FS):round(end*TRACE_FS)].astype('f8')
    mouth=measure(steady,fs=TRACE_FS,start=0.)
    pressure_error=steady[:,0]-steady[:,8]-RHO*C/(.00023)*steady[:,1]
    pitch=audio['frequency'];target=s['target_hz']
    return dict(note=s['note'],target_hz=target,audio_frequency_hz=pitch,
        cents=float(1200*np.log2(pitch/target)) if pitch>0 else None,
        periodicity=audio['periodicity'],mouthpiece_frequency_hz=mouth['frequency'],
        mouthpiece_rms_pa=mouth['rms'],return_characteristic_rms_pa=float(np.std(steady[:,8])),
        lip_y_rms_mm=float(np.std(steady[:,3])*1000),boundary_residual_max_pa=float(np.max(abs(pressure_error))),
        radiated_rms_pa=audio['rms'],blowing_pressure_pa=s['blowing_pressure_pa'])


def run(dz,dr,only,probe):
    os.makedirs(OUT,exist_ok=True)
    device=Device()
    signature=hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in [
        'src/render_axisymmetric.py','src/axisymmetric_euler.cu','src/axisymmetric_mesh.py',
        'src/axisymmetric_cuda.py','src/trumpet.py','data/model.json','data/york_trumpet_fig7_52.csv'])).digest()
    waves=[];report=[];start=0.
    for s in settings():
        if only and s['note'] not in only:continue
        t=time.perf_counter();name=s['note'];duration=1.1 if name=='Bb5' else .72
        if probe:duration=probe
        key=hashlib.sha256(signature+json.dumps([s,duration,dz,dr],sort_keys=True).encode()).hexdigest()
        file=OUT/(f'{name}_coupled.npz' if not probe else f'{name}_probe.npz')
        old=dict(np.load(file)) if file.exists() else {}
        if str(old.get('key',''))==key:
            tr=old['trace'];details=json.loads(str(old['details']));wall_s=0.
        else:
            tr,details=synthesize(device,s,duration,hold=duration if probe else None,dz=dz,dr=dr,
                                 snapshot=not probe and name=='Bb4')
            wall_s=time.perf_counter()-t
            np.savez(file,key=key,trace=tr,fs=TRACE_FS,details=json.dumps(details))
        end=duration if probe else duration-.16
        row=analyze(s,tr,max(.12,end-.24),end)
        row.update(start_s=start,duration_s=duration,compute_wall_s=wall_s,**details)
        print(json.dumps(row),flush=True);report.append(row)
        waves.extend([radiation(tr),np.zeros(2880)]);start+=duration+.06
        file=OUT/('probe_measurements.json' if probe else 'measurements_in_progress.json')
        file.write_text(json.dumps(report,indent=2))
    if not probe and not only:
        audio=np.concatenate(waves)
        assert len(report)==15 and np.isfinite(audio).all() and np.max(abs(audio))>0
        gain=.89/np.max(abs(audio));pcm=np.round(gain*audio*32767).astype('i2')
        file=OUT/'trumpet_Bb3_Bb5_major_f_axisymmetric.wav'
        wavfile.write(file,48000,pcm)
        result=dict(model='axisymmetric 2D compressible Euler + Adachi-Sato 2D lips',gpu=device.name,
                    duration_s=len(audio)/48000,sample_rate_hz=48000,bit_depth=16,
                    common_gain_fs_per_pa=float(gain),clipped_samples=int(np.sum(abs(pcm.astype('i4'))>=32767)),
                    max_pitch_error_cents=max(abs(x['cents']) for x in report),notes=report)
        (OUT/'measurements.json').write_text(json.dumps(result,indent=2))
        print('Saved',file,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',nargs='+');parser.add_argument('--probe',type=float)
    parser.add_argument('--dz',type=float,default=.0035);parser.add_argument('--dr',type=float,default=.003)
    args=parser.parse_args();run(args.dz,args.dr,args.only,args.probe)
