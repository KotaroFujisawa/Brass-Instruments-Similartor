"""Render the Conn 8H and reconstructed Wick 5AL forte B-flat major scale."""
import argparse
import hashlib
import json
import os
import tempfile
import time
from functools import cache

from axisymmetric_cuda import Device, Run
from compare_trumpet_bell import harmonic_rms, write_wav
from render_bidirectional import measure_audio
from render_trombone_attack import ATTACK_PRESET
from render_trumpet_bell import build_load
from trombone_axisymmetric import LIPS, ROOT, TRACE_FS, refined_trombone_mesh
from trumpet import C, RHO, measure, np, signal

OUT = ROOT/'output/trombone_detailed_scale'
PRESET = ROOT/'data/performance_trombone_detailed_scale.json'
WAV_NAME = 'trombone_5AL_Bb2_Bb5_major_f_strong_attack.wav'
LOUDER_OUT = ROOT/'output/trombone_louder_forte'
LOUDER_PRESET = ROOT/'data/performance_trombone_louder_forte.json'
LOUDER_WAV_NAME = 'trombone_5AL_Bb2_Bb5_major_louder_f_strong_attack.wav'
STRONGER_OUT = ROOT/'output/trombone_stronger_forte'
STRONGER_PRESET = ROOT/'data/performance_trombone_stronger_forte.json'
STRONGER_WAV_NAME = 'trombone_5AL_Bb2_Bb5_major_stronger_f_strong_attack.wav'
AUDIO_FS = 48000
GAP = .09
UPPER_SCALE = dict(C5=72, D5=74, Eb5=75, F5=77, G5=79, A5=81, Bb5=82)
MODEL = dict(thermoviscous=True, radiation_order=2, reservoir_inflow=True)
PITCH_RANGE = dict(min_hz=55, max_hz=1850)


def settings(variant='standard'):
    notes = json.loads(ATTACK_PRESET.read_text())
    overrides = json.loads(PRESET.read_text())
    dynamics = []
    if variant in ('louder','stronger'):
        dynamics.append((LOUDER_PRESET,dict(pressure_factor=1.15)))
    if variant=='stronger':
        dynamics.append((STRONGER_PRESET,dict(
            pressure_factor=1.25,attack_s=.0065,attack_boost=.45,attack_decay_s=.028,
        )))
    dynamics = [(json.loads(path.read_text()) if path.exists() else {},defaults)
                for path,defaults in dynamics]
    for name, midi in UPPER_SCALE.items():
        notes.append(dict(note=name, midi=midi, target_hz=440*2**((midi-69)/12),
                          attack_s=.012, attack_boost=.22, attack_decay_s=.024,
                          rest_y_m=0., trim_m=0., wave_speed_reserve=1.75))
    for note in notes:
        override = overrides.get(note['note'],{})
        note = note | override | {
            'lip_parameters': note.get('lip_parameters',{}) | override.get('lip_parameters',{}),
        }
        for changes,defaults in dynamics:
            change = defaults | changes.get(note['note'],{})
            factor = change.pop('pressure_factor',1.15)
            note = note | change | {
                'lip_parameters': note['lip_parameters'] | change.get('lip_parameters',{}),
                'blowing_pressure_pa': note['blowing_pressure_pa']*factor,
                'pressure_factor': note.get('pressure_factor',1.)*factor,
            }
        yield note


@cache
def solver_signature():
    names = (
        'axisymmetric_cuda.py','axisymmetric_euler.cu','thermoviscous.py',
        'bell_exterior.py','render_trumpet_bell.py','axisymmetric_mesh.py',
        'trombone_axisymmetric.py','render_trombone_detailed_scale.py',
        'trumpet.py','render_bidirectional.py','compare_trumpet_bell.py',
    )
    return hashlib.sha256(b''.join((ROOT/'src'/name).read_bytes() for name in names)).digest()


def inspect(note, trace, wave, duration, mesh):
    end = duration-.16 if duration>.8 else duration
    start = max(.11,end-.20)
    audio = measure_audio(wave,start,end,**PITCH_RANGE)
    frequency = period_frequency = audio['frequency']
    odd_to_even = None
    if .485 < frequency/note['target_hz'] < .515:
        sustained = wave[round(start*AUDIO_FS):round(end*AUDIO_FS)]
        harmonics = harmonic_rms(sustained,AUDIO_FS,frequency,12)
        odd = max((h['rms_pa'] for h in harmonics if h['order']%2),default=0.)
        even = max((h['rms_pa'] for h in harmonics if not h['order']%2),default=0.)
        if even:
            odd_to_even = odd/even
            if odd_to_even < .15:
                frequency *= 2
    steady = trace[round(start*TRACE_FS):round(end*TRACE_FS)]
    area = np.pi*mesh['radius'][0]**2
    boundary = steady[:,0]-steady[:,8]-(RHO*C/area)*steady[:,1]
    early = measure_audio(wave,.115,min(.31,end),**PITCH_RANGE)
    attack = wave[round(.012*AUDIO_FS):round(.08*AUDIO_FS)]
    return dict(
        note=note['note'],target_hz=note['target_hz'],audio_frequency_hz=frequency or None,
        audio_period_frequency_hz=period_frequency or None,subharmonic_odd_to_even=odd_to_even,
        cents=float(1200*np.log2(frequency/note['target_hz'])) if frequency else None,
        periodicity=audio['periodicity'],early_frequency_hz=early['frequency'] or None,
        radiated_rms_pa=audio['rms'],mouthpiece_rms_pa=float(np.std(steady[:,0])),
        incoming_characteristic_rms_pa=float(np.std(steady[:,8])),
        boundary_residual_pa=float(np.max(abs(boundary))),
        attack_rms_pa=float(np.std(attack)),
        peak_blowing_pressure_pa=float(np.max(trace[:,6])),preset=note,
        duration_s=duration,sustained_window_s=[start,end],
    )


def solve(device, note, duration, initial_check, *, output=None):
    mesh = refined_trombone_mesh(note['position'],note.get('trim_m',0.))
    lips = LIPS | note['lip_parameters'] | {'cup_area':np.pi*mesh['radius'][0]**2}
    key = hashlib.sha256(solver_signature()+json.dumps(
        [note,lips,duration,MODEL],sort_keys=True).encode())
    for name in ('z','radius','rings'):
        key.update(mesh[name].tobytes())
    path = (OUT if output is None else output)/'cache'/f"{note['note']}_{key.hexdigest()[:16]}.npz"
    if path.exists():
        with np.load(path) as data:
            return data['wave'],json.loads(str(data['row']))

    began = time.perf_counter()
    load = build_load(mesh)
    with Run(device,mesh,duration,lip_hz=note['lip_frequency_hz'],
             pressure=note['blowing_pressure_pa'],rest_y=note['rest_y_m'],
             attack=note['attack_s'],attack_boost=note['attack_boost'],
             attack_decay=note['attack_decay_s'],guard=note.get('wave_speed_reserve',1.58),
             hold=duration-.14 if duration>.8 else duration,lip_parameters=lips,
             exterior_load=load,record_bell=True,**MODEL) as run:
        if initial_check and duration>.8:
            run.until(.38)
            pilot = measure(run.history(),fs=TRACE_FS,start=.21,**PITCH_RANGE)
            relative = pilot['frequency']/note['target_hz']
            if .48<relative<.52:
                relative *= 2
            if pilot['rms']<30 or not .968<relative<1.032:
                raise ValueError(f"Attack reached {pilot['frequency']:.2f} Hz, {pilot['rms']:.1f} Pa cup RMS")
        run.until()
        trace,velocity = run.history(),run.bell_history()
        details = dict(cells=len(mesh['cells']),axial_cells=len(mesh['z'])-1,
                       time_step_s=run.dt,steps=run.steps,bell_rings=int(mesh['rings'][-1]),
                       bore_length_m=float(mesh['z'][-1]),bell_radius_m=float(mesh['radius'][-1]),
                       reconstruction_fallbacks=int(run.fallbacks),
                       compute_wall_s=time.perf_counter()-began)
    if not (np.isfinite(trace).all() and np.isfinite(velocity).all()):
        raise RuntimeError('Nonfinite solver output')
    areas = mesh['faces'][mesh['bell']:]['area']
    details['bell_flux_residual_m3_s'] = float(np.max(abs(velocity@areas-trace[:,9])))
    if details['bell_flux_residual_m3_s'] >= 2e-7:
        raise RuntimeError('Bell flow does not match the annular velocities')
    modal = velocity@(load.basis*load.weights[:,None])
    radiated = load.radiate(modal,TRACE_FS,angles_deg=(0.,),distance_m=1.)[:,0]
    wave = signal.resample_poly(radiated,AUDIO_FS,TRACE_FS)
    row = inspect(note,trace,wave,duration,mesh) | details
    with tempfile.NamedTemporaryFile(dir=path.parent,suffix='.npz',delete=False) as file:
        tmp_path = file.name
        np.savez(file,wave=wave,row=json.dumps(row,allow_nan=False))
    os.replace(tmp_path,path)
    return wave,row


def render(only, probe, initial_check, *, variant='standard'):
    output,wav_name,reference_path = {
        'standard': (OUT,WAV_NAME,None),
        'louder': (LOUDER_OUT,LOUDER_WAV_NAME,OUT/'measurements.json'),
        'stronger': (STRONGER_OUT,STRONGER_WAV_NAME,LOUDER_OUT/'measurements.json'),
    }[variant]
    notes = [note for note in settings(variant) if not only or note['note'] in only]
    os.makedirs(output/'cache',exist_ok=True)
    device = Device()
    rows,errors,waves = [],[],[]
    offset = 0
    for note in notes:
        name = note['note']
        duration = probe or (1.65 if name=='Bb5' else 1.08)
        try:
            wave,row = solve(device,note,duration,initial_check,output=output)
        except (RuntimeError,ValueError) as error:
            print(f'{name}: {error}',flush=True)
            errors.append(dict(note=name,error=str(error),preset=note))
            continue
        rows.append(row | {'start_s': offset/AUDIO_FS})
        waves.append(wave)
        offset += len(wave)+round(GAP*AUDIO_FS)
        print(f"{name:>3}  {row['audio_frequency_hz']} Hz  {row['cents']} cents  "
              f"p={row['periodicity']:.3f}  bell={row['radiated_rms_pa']:.3f} Pa  "
              f"cup={row['mouthpiece_rms_pa']:.0f} Pa  {row['compute_wall_s']:.1f} s",flush=True)
    if probe or only:
        suffix = '_'.join(only) if only else 'all'
        label = f'probe_{probe:g}' if probe else 'part'
        (output/f'{label}_{suffix}.json').write_text(
            json.dumps(dict(notes=rows,errors=errors),indent=2,allow_nan=False)+'\n')
        return
    (output/'progress.json').write_text(
        json.dumps(dict(notes=rows,errors=errors),indent=2,allow_nan=False)+'\n')
    failed = [row['note'] for row in rows if row['cents'] is None or abs(row['cents'])>15
              or row['periodicity']<.90 or row['radiated_rms_pa']<.12]
    if errors or failed or len(rows)!=22:
        raise RuntimeError(f"Review these notes: {', '.join([e['note'] for e in errors]+failed)}")
    gap = np.zeros(round(GAP*AUDIO_FS))
    # Fade the final 5 ms of ring-down to avoid a discontinuity at the silent gap.
    release = .5+.5*np.cos(np.linspace(0.,np.pi,round(.005*AUDIO_FS)))
    for wave in waves:
        wave[-len(release):] *= release
    audio = np.concatenate([part for wave in waves for part in (wave,gap)])
    if len(audio)!=1262880 or not np.isfinite(audio).all() or np.max(abs(audio))==0:
        raise RuntimeError('Invalid assembled scale')
    peak = float(np.max(abs(audio)))
    gain = .89/peak
    if reference_path is not None and reference_path.exists():
        reference = json.loads(reference_path.read_text())
        gain = min(.97/peak,reference['common_gain_fs_per_pa'])
        by_note = {row['note']:row for row in reference['notes']}
        for row in rows:
            previous = by_note[row['note']]['radiated_rms_pa']
            row['change_in_sustained_level_db'] = float(20*np.log10(row['radiated_rms_pa']/previous))
    write_wav(output,wav_name,audio,gain)
    report = dict(
        gpu=device.name,duration_s=len(audio)/AUDIO_FS,sample_rate_hz=AUDIO_FS,bit_depth=16,
        channels=1,common_gain_fs_per_pa=gain,clipped_samples=0,
        mouthpiece='Denis Wick 5AL reconstructed from The Science of Brass Instruments, Fig. 4.32b',
        bore='Conn 8H reconstructed from The Science of Brass Instruments, Fig. 4.22',
        bell='coupled spherical exterior; receiver at 0 deg and 1 m',
        output_tail_fade_s=.005,
        max_sustained_pitch_error_cents=max(abs(row['cents']) for row in rows),notes=rows,
    )
    (output/'measurements.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',nargs='+',choices=[note['note'] for note in settings()])
    parser.add_argument('--probe',type=float)
    levels = parser.add_mutually_exclusive_group()
    levels.add_argument('--louder',dest='variant',action='store_const',const='louder',
                        help='Render the louder forte variant.')
    levels.add_argument('--stronger',dest='variant',action='store_const',const='stronger',
                        help='Raise the pressure further and accent the attack.')
    parser.set_defaults(variant='standard')
    parser.add_argument('--no-initial-check',action='store_true')
    args = parser.parse_args()
    if args.probe is not None and not .23<=args.probe<=.8:
        parser.error('--probe duration must be between .23 s and .8 s')
    render(args.only,args.probe,not args.no_initial_check,variant=args.variant)
