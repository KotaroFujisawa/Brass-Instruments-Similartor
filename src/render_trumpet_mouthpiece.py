"""Compare coupled playing with two mouthpiece shapes and a finer mouthpiece grid."""
import argparse
import hashlib
import json
import os

from compare_trumpet_bell import AUDIO_FS, analyze_record, write_wav
from render_trumpet_bell import ROOT, TRACE_FS, Device, Run, build_load, settings, simulate
from axisymmetric_mesh import refined_trumpet_mesh
from trumpet import C, RHO, np, signal
from scipy.fft import rfft, rfftfreq

OUT = ROOT/'output/trumpet_mouthpiece'
NOTES = ('Bb4','Bb6')
VARIANTS = ('legacy','legacy_fine','detailed')
DURATION = .8
WINDOW = (.38,.62)
MOUTHPIECE_SCALE = .6


def mesh_for(variant, note):
    return refined_trumpet_mesh(
        note['valves_semitones'], tuning=note.get('tuning_m',0.),
        mouthpiece_profile='detailed' if variant=='detailed' else 'legacy',
        mouthpiece_scale=1. if variant=='legacy' else MOUTHPIECE_SCALE,
    )


def signature(mesh, note):
    files = ('axisymmetric_cuda.py','axisymmetric_euler.cu','thermoviscous.py',
             'bell_exterior.py','render_trumpet_bell.py','render_trumpet_high.py',
             'render_trumpet_mouthpiece.py')
    key = hashlib.sha256(b''.join((ROOT/'src'/name).read_bytes() for name in files))
    for name in ('z','radius','rings'):
        key.update(mesh[name].tobytes())
    key.update(json.dumps([note,DURATION],sort_keys=True).encode())
    return key.hexdigest()[:16]


def run_note(device, note, variant):
    mesh = mesh_for(variant,note)
    load = build_load(mesh)
    area = np.pi*mesh['radius'][0]**2
    physical_note = note | {'lip_parameters':note.get('lip_parameters',{}) | {'cup_area':area}}
    cache = OUT/'traces'/f"{variant}_{note['note']}_{signature(mesh,physical_note)}.npz"
    if not cache.exists():
        trace,velocity,details = simulate(device,physical_note,mesh,load,DURATION)
        pending = cache.with_suffix('.pending.npz')
        np.savez(pending,trace=trace,velocity=velocity,details=json.dumps(details))
        pending.replace(cache)
    with np.load(cache) as data:
        trace,velocity = data['trace'],data['velocity']
        details = json.loads(str(data['details']))
    modal_velocity = velocity@(load.basis*load.weights[:,None])
    pressure = load.radiate(modal_velocity,TRACE_FS,angles_deg=(0,45),distance_m=1.)
    scaled = signal.resample_poly(pressure,AUDIO_FS,TRACE_FS,axis=0)
    start,end = (round(x*TRACE_FS) for x in WINDOW)
    sustain = trace[start:end]
    details.update(cup_inlet_area_m2=area,
                   boundary_residual_pa=float(np.max(abs(sustain[:,0]-sustain[:,8]-(RHO*C/area)*sustain[:,1]))),
                   preset=physical_note)
    row,audio = analyze_record(note['note'],variant,dict(
        trace=trace,audio={0:scaled[:,0],45:scaled[:,1]},target_hz=note['target_hz'],
        duration_s=DURATION,window_s=WINDOW,diagnostics=details))
    print(f"{note['note']} {variant:11} pitch={row['frequency_hz']} Hz "
          f"cup={row['cup_rms_pa']:.1f} Pa  {details['cells']} cells  "
          f"{details['compute_wall_s']:.1f} s",flush=True)
    return row,audio


def pulse_impedance(device, variant, note):
    """Find acoustic input impedance of the same exterior-loaded Euler grid."""
    mesh = mesh_for(variant,note)
    load = build_load(mesh)
    amplitude = 60.
    with Run(device,mesh,.48,mode=3,pressure=amplitude,thermoviscous=True,
             reservoir_inflow=True,exterior_load=load,guard=1.58) as run:
        run.until()
        trace = run.history()
        dt = run.dt
    area = np.pi*mesh['radius'][0]**2
    time = (np.arange(len(trace))+1)/TRACE_FS-dt
    drive_p = amplitude*np.exp(-.5*((time-.00045)/.00007)**2)
    volume_flow = area/(RHO*C)*drive_p
    pressure = trace[:,8]+drive_p
    nfft = 2**19
    frequency = rfftfreq(nfft,1/TRACE_FS)
    band = (frequency>=40)&(frequency<=2600)
    z = rfft(pressure,nfft)[band]/rfft(volume_flow,nfft)[band]
    details = dict(cells=len(mesh['cells']),time_step_s=dt,
                   tail_mean_pressure_pa=float(np.mean(pressure[-TRACE_FS//20:])))
    print(f'Pulse {variant:11} {details}',flush=True)
    return frequency[band],z,details


def record_impedance(device, variants, note):
    for variant in variants:
        f,z,details = pulse_impedance(device,variant,note)
        np.savez(OUT/f'impedance_{variant}.npz',frequency_hz=f,impedance_pa_s_m3=z,
                 diagnostics=json.dumps(details))


def write_comparison(rows, waves):
    ordered = [(note,variant) for note in NOTES for variant in VARIANTS]
    peak = max(float(np.max(abs(waves[key][0]))) for key in ordered)
    gain = .89/peak
    gap = np.zeros(round(.15*AUDIO_FS))
    separator = np.zeros(round(.5*AUDIO_FS))
    blocks = []
    for variant in VARIANTS:
        blocks.append(np.concatenate([waves[NOTES[0],variant][0],gap,waves[NOTES[1],variant][0]]))
    comparison = np.concatenate([blocks[0],separator,blocks[1],separator,blocks[2]])
    write_wav(OUT,'trumpet_mouthpiece_comparison.wav',comparison,gain)
    write_wav(OUT,'trumpet_mouthpiece_detailed.wav',blocks[-1],gain)
    report = dict(notes=[rows[key] for key in ordered],variant_order=VARIANTS,
                  note_order=NOTES,common_gain_fs_per_pa=gain,duration_s=len(comparison)/AUDIO_FS,
                  sample_rate_hz=AUDIO_FS)
    (OUT/'mouthpiece_playing.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


def render(notes, variants, impedance):
    os.makedirs(OUT/'traces',exist_ok=True)
    device = Device()
    performance = {note['note']:note for note in settings()}
    if impedance:
        record_impedance(device,variants,dict(note='open',valves_semitones=[]))
        return
    rows,waves = {},{}
    for variant in variants:
        for note in notes:
            rows[note,variant],waves[note,variant] = run_note(device,performance[note],variant)
    if tuple(notes)==NOTES and tuple(variants)==VARIANTS:
        write_comparison(rows,waves)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--notes',nargs='+',choices=NOTES,default=NOTES)
    parser.add_argument('--variants',nargs='+',choices=VARIANTS,default=VARIANTS)
    parser.add_argument('--impedance',action='store_true')
    args = parser.parse_args()
    render(args.notes,args.variants,args.impedance)
