"""Synthesize the three-octave forte trumpet scale with the digitized Bach 5C."""
import argparse
import hashlib
import json
import os
import tempfile
import time
from functools import cache

from axisymmetric_cuda import Device, Run
from axisymmetric_mesh import refined_trumpet_mesh
from compare_trumpet_bell import harmonic_rms, write_wav
from render_bidirectional import measure_audio
from render_trumpet_bell import build_load
from render_trumpet_high import AUDIO_FS, GAP_S, MODEL, ROOT, TRACE_FS, settings
from render_trumpet_mouthpiece import MOUTHPIECE_SCALE
from trumpet import C, RHO, measure, np, signal

OUT = ROOT/'output/trumpet_bach5c'
PRESET = ROOT/'data/performance_trumpet_bach5c.json'
WAV_NAME = 'trumpet_Bach5C_Bb3_Bb6_major_f_strong_attack.wav'
FRAMES = 1262880
RIM_SCALE = 2.01


def physical_note(note):
    mesh = refined_trumpet_mesh(
        note['valves_semitones'],tuning=note.get('tuning_m',0.),
        mouthpiece_profile='detailed',mouthpiece_scale=MOUTHPIECE_SCALE,
        mouthpiece_rim_scale=note.get('mouthpiece_rim_scale',RIM_SCALE),
    )
    area = np.pi*mesh['radius'][0]**2
    lips = note.get('lip_parameters',{}) | {'cup_area':area}
    return note | {'lip_parameters':lips},mesh


@cache
def source_signature():
    names = ('axisymmetric_cuda.py','axisymmetric_euler.cu','thermoviscous.py','bell_exterior.py',
             'render_trumpet_mouthpiece.py')
    return hashlib.sha256(b''.join((ROOT/'src'/name).read_bytes() for name in names)).digest()


def save_trace(path, trace, velocity, details):
    with tempfile.NamedTemporaryFile(dir=path.parent,suffix='.npz',delete=False) as file:
        name = file.name
        np.savez(file,trace=trace,velocity=velocity,details=json.dumps(details))
    os.replace(name,path)


def solve(device, note, seconds, initial_check=True):
    note,mesh = physical_note(note)
    load = build_load(mesh)
    key = hashlib.sha256(source_signature()+json.dumps(
        [note,seconds,MODEL],sort_keys=True).encode())
    for component in ('z','radius','rings'):
        key.update(mesh[component].tobytes())
    path = OUT/'traces'/f"{note['note']}_{key.hexdigest()[:16]}.npz"
    if path.exists():
        with np.load(path) as data:
            return (note,data['trace'],data['velocity'],json.loads(str(data['details'])),load)
    began = time.perf_counter()
    lips = MODEL['lip_parameters'] | note['lip_parameters']
    model = MODEL | {'lip_parameters':lips}
    with Run(device,mesh,seconds,lip_hz=note['lip_frequency_hz'],
             pressure=note['blowing_pressure_pa'],rest_y=note['rest_y_m'],
             attack=note['attack_s'],attack_boost=note['attack_boost'],
             attack_decay=note['attack_decay_s'],guard=note.get('wave_speed_reserve',1.58),
             hold=seconds-.14 if seconds>.8 else seconds,exterior_load=load,record_bell=True,**model) as run:
        if initial_check and seconds>.8:
            run.until(.40)
            trace = run.history()
            pilot = measure(trace,fs=TRACE_FS,start=.25,max_hz=2450)
            relative = pilot['frequency']/note['target_hz']
            if .48<relative<.52:
                relative *= 2
            if pilot['rms']<45 or not .982<relative<1.018:
                raise ValueError(f"attack settled at {pilot['frequency']:.2f} Hz, {pilot['rms']:.1f} Pa RMS")
        run.until()
        trace,velocity = run.history(),run.bell_history()
        details = dict(cells=len(mesh['cells']),time_step_s=run.dt,steps=run.steps,
                       bell_rings=int(mesh['rings'][-1]),reconstruction_fallbacks=int(run.fallbacks),
                       compute_wall_s=time.perf_counter()-began,
                       bell_flux_residual_m3_s=float(np.max(abs(
                           velocity@mesh['faces'][mesh['bell']:]['area']-trace[:,9]))))
    if not (np.isfinite(trace).all() and np.isfinite(velocity).all()):
        raise RuntimeError('Nonfinite CUDA trace')
    save_trace(path,trace,velocity,details)
    return note,trace,velocity,details,load


def levels(wave,window):
    start,end = (round(t*AUDIO_FS) for t in window)
    return float(np.std(wave[start:end]))


def inspect(note, trace, wave, duration):
    end = duration-.16 if duration>.8 else duration
    start = max(.12,end-.20)
    sustained = wave[round(start*AUDIO_FS):round(end*AUDIO_FS)]
    audio = measure_audio(wave,start,end,max_hz=2500)
    frequency = audio['frequency']
    physical_period_hz = frequency
    if frequency and .485<frequency/note['target_hz']<.515:
        harmonics = harmonic_rms(sustained,AUDIO_FS,frequency,12)
        odd = max((h['rms_pa'] for h in harmonics if h['order']%2),default=0.)
        even = max((h['rms_pa'] for h in harmonics if not h['order']%2),default=0.)
        if even and odd/even<.15:
            frequency *= 2
    early = measure_audio(wave,.21,min(.41,end),max_hz=2500)
    steady_trace = trace[round(start*TRACE_FS):round(end*TRACE_FS)]
    area = note['lip_parameters']['cup_area']
    boundary = steady_trace[:,0]-steady_trace[:,8]-(RHO*C/area)*steady_trace[:,1]
    return dict(note=note['note'],target_hz=note['target_hz'],audio_frequency_hz=frequency or None,
                audio_period_frequency_hz=physical_period_hz or None,
                cents=float(1200*np.log2(frequency/note['target_hz'])) if frequency else None,
                periodicity=audio['periodicity'],early_frequency_hz=early['frequency'] or None,
                radiated_rms_pa=audio['rms'],mouthpiece_rms_pa=float(np.std(steady_trace[:,0])),
                boundary_residual_pa=float(np.max(abs(boundary))),
                attack_rms_pa=levels(wave,(.012,.08)),
                peak_blowing_pressure_pa=float(np.max(trace[:,6])),preset=note,duration_s=duration,
                sustained_window_s=[start,end])


def render(only=None, probe=None, *, initial_check=True):
    os.makedirs(OUT/'traces',exist_ok=True)
    overrides = json.loads(PRESET.read_text()) if PRESET.exists() else {}
    notes = settings()
    device = Device()
    rows,errors,waves = [],[],[]
    offset = 0
    for original in notes:
        name = original['note']
        if only and name not in only:
            continue
        override = overrides.get(name,{})
        note = original | override | {'lip_parameters': original.get('lip_parameters',{}) |
                                                            override.get('lip_parameters',{})}
        seconds = probe or (1.65 if name=='Bb6' else 1.08)
        try:
            note,trace,velocity,details,load = solve(device,note,seconds,initial_check)
            modal = velocity@(load.basis*load.weights[:,None])
            radiated = load.radiate(modal,TRACE_FS,angles_deg=(0.,),distance_m=1.)[:,0]
            wave = signal.resample_poly(radiated,AUDIO_FS,TRACE_FS)
            row = inspect(note,trace,wave,seconds)
        except (ValueError,RuntimeError) as error:
            print(f'{name}: {error}',flush=True)
            errors.append(dict(note=name,error=str(error),preset=note))
            continue
        row.update(details,start_s=offset/AUDIO_FS)
        rows.append(row)
        offset += len(wave)+round(GAP_S*AUDIO_FS)
        waves.append(wave)
        print(f"{name:>3}  {row['audio_frequency_hz']} Hz  {row['cents']} cents  "
              f"p={row['periodicity']:.3f}  cup={row['mouthpiece_rms_pa']:.0f} Pa  "
              f"{details['compute_wall_s']:.1f} s",flush=True)
        if not only and not probe:
            (OUT/'progress.json').write_text(json.dumps(dict(notes=rows,errors=errors),indent=2))
    suffix = '_'+('_'.join(only) if only else 'all')
    if only or probe:
        (OUT/f"measurements{'_probe' if probe else ''}{suffix}.json").write_text(
            json.dumps(dict(notes=rows,errors=errors),indent=2))
        return
    failed = [row['note'] for row in rows if row['cents'] is None or abs(row['cents'])>15
              or row['periodicity']<.90 or row['radiated_rms_pa']<.12]
    if errors or failed or len(rows)!=22:
        raise RuntimeError(f"Review these notes: {', '.join([e['note'] for e in errors]+failed)}")
    gap = np.zeros(round(GAP_S*AUDIO_FS))
    audio = np.concatenate([part for wave in waves for part in (wave,gap)])
    assert len(audio)==FRAMES and np.isfinite(audio).all() and np.max(abs(audio))>0
    gain = float(.89/np.max(abs(audio)))
    write_wav(OUT,WAV_NAME,audio,gain)
    report = dict(gpu=device.name,duration_s=len(audio)/AUDIO_FS,sample_rate_hz=AUDIO_FS,bit_depth=16,
                  channels=1,common_gain_fs_per_pa=gain,clipped_samples=0,
                  mouthpiece='Bach 5C, Jackson & Scavone (2025) Fig. 2',
                  bell='coupled spherical exterior; receiver at 0 deg and 1 m',
                  max_sustained_pitch_error_cents=max(abs(row['cents']) for row in rows),notes=rows)
    (OUT/'measurements.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',nargs='+',choices=[note['note'] for note in settings()])
    parser.add_argument('--probe',type=float)
    parser.add_argument('--no-initial-check',action='store_true')
    args = parser.parse_args()
    if args.probe is not None and not .36<=args.probe<=.8:
        parser.error('--probe duration must be between .36 s and .8 s')
    render(args.only,args.probe,initial_check=not args.no_initial_check)
