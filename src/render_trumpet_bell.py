"""Compare finite trumpet bell radiation and the refined bore on three fortes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

from axisymmetric_cuda import Device, Run
from axisymmetric_mesh import refined_trumpet_mesh, trumpet_mesh
from render_axisymmetric import TRACE_FS, analyze, radiation
from render_trumpet_high import AUDIO_FS, MODEL, ROOT, settings
from trumpet import np, signal

OUT = ROOT/'output/trumpet_bell'
NOTES = ('Bb4','Bb5','Bb6')
VARIANTS = ('previous','fine_taper','exterior')
DURATION = .80
WINDOW = (.38,.62)


def build_load(mesh):
    from bell_exterior import exterior_model

    edges = np.linspace(0.,mesh['radius'][-1],int(mesh['rings'][-1])+1)
    return exterior_model(edges)


def simulate(device, note, mesh, load, duration):
    start = time.perf_counter()
    model = MODEL | {'lip_parameters': MODEL['lip_parameters'] | note.get('lip_parameters',{})}
    with Run(device,mesh,duration,lip_hz=note['lip_frequency_hz'],
             pressure=note['blowing_pressure_pa'],rest_y=note['rest_y_m'],
             attack=note['attack_s'],attack_boost=note['attack_boost'],
             attack_decay=note['attack_decay_s'],guard=note.get('wave_speed_reserve',1.58),
             exterior_load=load,record_bell=load is not None,**model) as run:
        run.until()
        trace = run.history()
        bell_velocity = run.bell_history() if load is not None else np.zeros((0,0))
        details = dict(cells=len(mesh['cells']),axial_cells=len(mesh['z'])-1,steps=run.steps,
                       bell_rings=int(mesh['rings'][-1]),time_step_s=run.dt,
                       bell_radius_m=float(mesh['radius'][-1]),
                       bore_length_m=float(mesh['z'][-1]),
                       reconstruction_fallbacks=int(run.fallbacks),
                       compute_wall_s=time.perf_counter()-start)
    if load is not None:
        face_areas = mesh['faces'][mesh['bell']:]['area']
        flow_error = np.max(abs(bell_velocity@face_areas-trace[:,9]))
        if flow_error>=2e-7:
            raise RuntimeError(f'Bell flux and annular velocities disagree ({flow_error} m3/s)')
        details['bell_flux_residual_m3_s'] = float(flow_error)
    if not (np.isfinite(trace).all() and np.isfinite(bell_velocity).all()):
        raise RuntimeError(f'Nonfinite pressure or flow in {note["note"]}')
    return trace,bell_velocity,details


def source_signature(variant):
    names = ['src/render_trumpet_bell.py','src/axisymmetric_cuda.py','src/axisymmetric_euler.cu',
             'src/axisymmetric_mesh.py','src/trumpet.py','src/thermoviscous.py',
             'data/model.json','data/york_trumpet_fig7_52.csv']
    if variant=='exterior':
        names.append('src/bell_exterior.py')
    return hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in names)).digest()


def make_record(device, note, variant, output, duration=DURATION, axial_scale=1.):
    os.makedirs(output/'traces',exist_ok=True)
    options = dict(tuning=note.get('tuning_m',0.))
    if variant=='previous':
        mesh = trumpet_mesh(note['valves_semitones'],**options)
    else:
        mesh = refined_trumpet_mesh(note['valves_semitones'],axial_scale=axial_scale,**options)
    load = build_load(mesh) if variant=='exterior' else None
    key = hashlib.sha256(source_signature(variant)+json.dumps(
        [note,variant,duration,axial_scale,MODEL],sort_keys=True).encode()).hexdigest()
    cache = output/'traces'/f'{variant}_{note["note"]}_{key[:16]}.npz'
    if not cache.exists():
        trace,velocity,details = simulate(device,note,mesh,load,duration)
        staging = cache.with_suffix('.pending.npz')
        np.savez(staging,trace=trace,bell_velocity_m_s=velocity,diagnostics=json.dumps(details))
        staging.replace(cache)
    with np.load(cache) as data:
        trace,velocity = data['trace'],data['bell_velocity_m_s']
        details = json.loads(str(data['diagnostics']))
    if load is None:
        audio = {0:radiation(trace)}
    else:
        modal = velocity@(load.basis*load.weights[:,None])
        pressures = load.radiate(modal,TRACE_FS,angles_deg=(0,45),distance_m=1.)
        down = signal.resample_poly(pressures,AUDIO_FS,TRACE_FS,axis=0)
        audio = {0:down[:,0],45:down[:,1]}
    details['preset'] = note
    comparison_window = WINDOW if duration==DURATION else (duration-.34,duration-.16)
    old_radiation = analyze(note,trace,*comparison_window,max_hz=2500)
    print(f'{variant:11} {note["note"]}  {old_radiation["mouthpiece_frequency_hz"]:.2f} Hz  '
          f'{old_radiation["mouthpiece_rms_pa"]:.0f} Pa RMS  '
          f'{details["cells"]} cells  {details["compute_wall_s"]:.1f} s',flush=True)
    return dict(trace=trace,audio=audio,target_hz=note['target_hz'],duration_s=duration,
                window_s=comparison_window,diagnostics=details)


def render(output, notes, variants, duration, axial_scale):
    device, records = Device(), {}
    presets = {s['note']:s for s in settings()}
    for variant in variants:
        for name in notes:
            records[name,variant] = make_record(device,presets[name],variant,output,duration,axial_scale)
    if tuple(notes)==NOTES and tuple(variants)==VARIANTS:
        from compare_trumpet_bell import write_comparison

        write_comparison(records,output)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--notes',nargs='+',choices=NOTES,default=NOTES)
    parser.add_argument('--variants',nargs='+',choices=VARIANTS,default=VARIANTS)
    parser.add_argument('--duration',type=float,default=DURATION)
    parser.add_argument('--axial-scale',type=float,default=1.)
    options = parser.parse_args()
    if options.duration<=.4:
        parser.error('--duration must exceed 0.4 s for comparison after the attack')
    render(options.output,options.notes,options.variants,options.duration,options.axial_scale)
