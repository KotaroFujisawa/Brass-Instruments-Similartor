"""Render a three-octave concert B-flat major trumpet scale with accented notes."""
import argparse
import hashlib
import json
import os
import time

from render_axisymmetric import analyze, radiation
from render_trumpet_modern import AUDIO_FS, GAP_S, MODEL, ROOT, TRACE_FS, Device, performance, simulate
from trumpet import np
from scipy.io import wavfile

OUT = ROOT/'output/trumpet_high'
PRESET = ROOT/'data/performance_trumpet_high.json'
WAV_NAME = 'trumpet_Bb3_Bb6_major_f_strong_attack.wav'
ATTACK = dict(attack_s=.014, attack_boost=.22, attack_decay_s=.027)
UPPER_MIDI = [84,86,87,89,91,93,94]


def settings():
    preset = json.loads(PRESET.read_text())
    lower = [s | ATTACK | preset['overrides'].get(s['note'], {}) for s in performance()]
    upper = [ATTACK | s for s in preset['upper_octave']]
    return lower+upper


def note_trace(device, note, duration, probe, signature):
    key = hashlib.sha256(signature+json.dumps([note,duration,probe,MODEL],sort_keys=True).encode()).hexdigest()
    cache = OUT/'traces'/f'{note["note"]}_{key[:16]}.npz'
    if not cache.exists():
        trace, details = simulate(device,note,duration,probe=probe,**MODEL)
        np.savez(cache,trace=trace,details=json.dumps(details),sample_rate_hz=TRACE_FS)
    with np.load(cache) as data:
        return data['trace'], json.loads(str(data['details']))


def verify_pitch(note, trace, end):
    late = analyze(note,trace,end-.20,end,max_hz=2450)
    early = analyze(note,trace,.17,.37,max_hz=2450)
    return late,early


def render(only=None, probe=None):
    notes = settings()
    if not only and [s['midi'] for s in notes[-7:]]!=UPPER_MIDI:
        raise ValueError('The preset must contain the full upper octave.')
    os.makedirs(OUT/'traces',exist_ok=True)
    signature = hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in (
        'src/render_trumpet_modern.py','src/axisymmetric_cuda.py','src/axisymmetric_euler.cu',
        'src/axisymmetric_mesh.py','src/thermoviscous.py','src/trumpet.py',
        'data/model.json','data/york_trumpet_fig7_52.csv',
    ))).digest()
    device = Device()
    waves,report,offset = [],[],0
    for s in notes:
        if only and s['note'] not in only:
            continue
        duration = probe or (1.65 if s['midi']==94 else 1.08)
        started = time.perf_counter()
        trace,details = note_trace(device,s,duration,probe is not None,signature)
        late,early = verify_pitch(s,trace,duration if probe else duration-.16)
        wave = radiation(trace)
        pm_peak = float(np.max(trace[:,6]))
        attack_rms = float(np.std(wave[round(.012*AUDIO_FS):round(.08*AUDIO_FS)]))
        late.update(
            start_s=offset/AUDIO_FS,duration_s=duration,early_cents=early['cents'],
            early_periodicity=early['periodicity'],attack_rms_pa=attack_rms,
            peak_blowing_pressure_pa=pm_peak,preset=s,**details,
        )
        print(f"{s['note']:>3} {late['audio_frequency_hz']:7.2f} Hz  "
              f"{str(late['cents']):>19} cents  p={late['periodicity']:.3f}  "
              f"early={early['cents']}  {time.perf_counter()-started:.1f} s",flush=True)
        offset += len(wave)+round(GAP_S*AUDIO_FS)
        waves.extend([wave,np.zeros(round(GAP_S*AUDIO_FS))])
        report.append(late)
    if only or probe:
        (OUT/('probe_measurements.json' if probe else 'partial_measurements.json')).write_text(json.dumps(report,indent=2))
        return
    invalid = [row['note'] for row in report if row['cents'] is None or abs(row['cents'])>15
               or row['radiated_rms_pa']<.12 or row['periodicity']<.90]
    if invalid:
        raise RuntimeError(f'Sustained pitch or level needs review for {", ".join(invalid)}')
    audio = np.concatenate(waves)
    expected_frames = round((21*1.08+1.65+22*GAP_S)*AUDIO_FS)
    assert len(audio)==expected_frames and np.isfinite(audio).all() and np.max(abs(audio))>0
    gain = .89/np.max(abs(audio))
    pcm = np.round(gain*audio*32767).astype('i2')
    wavfile.write(OUT/WAV_NAME,AUDIO_FS,pcm)
    result = dict(gpu=device.name,duration_s=len(audio)/AUDIO_FS,sample_rate_hz=AUDIO_FS,bit_depth=16,
                  common_gain_fs_per_pa=gain,clipped_samples=int(np.count_nonzero(abs(pcm.astype('i4'))>=32767)),
                  max_sustained_pitch_error_cents=max(abs(row['cents']) for row in report),
                  notes=report)
    (OUT/'measurements.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',nargs='+')
    parser.add_argument('--probe',type=float)
    args = parser.parse_args()
    if args.probe is not None and args.probe<.37:
        parser.error('--probe must be at least 0.37 seconds for pitch measurements')
    render(args.only,args.probe)
