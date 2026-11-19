"""Render the longer Bb-major scale with the updated coupled trumpet model."""
import argparse
import hashlib
import os
import time

from render_axisymmetric import Device, Run, TRACE_FS, analyze, radiation, settings, trumpet_mesh
from trumpet import ROOT, json, np
from scipy.io import wavfile

OUT = ROOT/'output/trumpet_modern'
PRESET = ROOT/'data/performance_trumpet_modern.json'
WAV_NAME = 'trumpet_Bb3_Bb5_major_f_modern_1p5x.wav'
AUDIO_FS = 48000
GAP_S = .09
LIPS = dict(gap_cubic_per_m2=8e4)
MODEL = dict(thermoviscous=True, radiation_order=2, reservoir_inflow=True, lip_parameters=LIPS)


def performance():
    changes = json.loads(PRESET.read_text()) if PRESET.exists() else {}
    return [s | changes.get(s['note'], {}) for s in settings()]


def simulate(device, s, duration, *, probe=False, feedback=True, **model):
    mesh = trumpet_mesh(s['valves_semitones'], tuning=s.get('tuning_m', 0.))
    t0 = time.perf_counter()
    hold = duration if probe else duration-.14
    lip_parameters = model.pop('lip_parameters', {}) | s.get('lip_parameters', {})
    with Run(device, mesh, duration, lip_hz=s['lip_frequency_hz'],
             pressure=s['blowing_pressure_pa'], rest_y=s['rest_y_m'],
             attack=s.get('attack_s', .045), attack_boost=s.get('attack_boost', 0.),
             attack_decay=s.get('attack_decay_s', .024), lip_parameters=lip_parameters,
             hold=hold, feedback=feedback, guard=s.get('wave_speed_reserve', 1.58), **model) as run:
        run.until()
        trace = run.history()
        if not np.isfinite(trace).all():
            raise RuntimeError(f"Nonfinite boundary history for {s['note']}")
        diagnostics = dict(
            cells=int(len(mesh['cells'])), dt_s=run.dt, steps=run.steps,
            reconstruction_fallbacks=int(run.fallbacks), wall_loss_fit_error=run.loss_error,
            compute_wall_s=time.perf_counter()-t0,
        )
    return trace, diagnostics


def cached_note(device, s, duration, probe, signature):
    key = hashlib.sha256(signature+json.dumps([s,duration,probe,MODEL], sort_keys=True).encode()).hexdigest()
    cache = OUT/'traces'/f"{s['note']}_{key[:16]}.npz"
    if not cache.exists():
        trace, details = simulate(device, s, duration, probe=probe, **MODEL)
        np.savez(cache, trace=trace, details=json.dumps(details), sample_rate_hz=TRACE_FS)
    with np.load(cache) as data:
        return data['trace'], json.loads(str(data['details']))


def render(only=None, probe=None):
    os.makedirs(OUT/'traces', exist_ok=True)
    device = Device()
    signature = hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in (
        'src/axisymmetric_cuda.py', 'src/axisymmetric_euler.cu', 'src/axisymmetric_mesh.py',
        'src/thermoviscous.py', 'src/trumpet.py', 'data/model.json', 'data/york_trumpet_fig7_52.csv',
    ))).digest()
    notes = performance()
    report, waves, offset = [], [], 0
    for s in notes:
        if only and s['note'] not in only:
            continue
        duration = probe or (1.65 if s['note'] == 'Bb5' else 1.08)
        trace, diagnostics = cached_note(device, s, duration, probe is not None, signature)
        end = duration if probe else duration-.16
        row = analyze(s, trace, max(.12, end-.20), end)
        early = analyze(s, trace, .21, min(.41,end))
        row.update(
            start_s=offset/AUDIO_FS, duration_s=duration, early_cents=early['cents'],
            early_periodicity=early['periodicity'], lip_frequency_hz=s['lip_frequency_hz'],
            tuning_m=s.get('tuning_m',0.), **diagnostics,
        )
        wave = radiation(trace)
        waves.extend([wave, np.zeros(round(GAP_S*AUDIO_FS))])
        offset += len(wave)+round(GAP_S*AUDIO_FS)
        report.append(row)
        print(json.dumps(row), flush=True)
    if probe or only:
        return
    audio = np.concatenate(waves)
    assert len(report)==15 and len(audio)==869760 and np.isfinite(audio).all() and np.max(abs(audio))>0
    gain = .89/np.max(abs(audio))
    pcm = np.round(gain*audio*32767).astype('i2')
    wavfile.write(OUT/WAV_NAME, AUDIO_FS, pcm)
    result = dict(gpu=device.name, duration_s=len(audio)/AUDIO_FS, sample_rate_hz=AUDIO_FS,
                  bit_depth=16, common_gain_fs_per_pa=gain,
                  clipped_samples=int(np.sum(abs(pcm.astype('i4'))>=32767)),
                  model_parameters=MODEL, notes=report)
    (OUT/'measurements.json').write_text(json.dumps(result, indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', nargs='+')
    parser.add_argument('--probe', type=float, help='Hold a short trial for the specified time in seconds.')
    args = parser.parse_args()
    render(args.only,args.probe)
