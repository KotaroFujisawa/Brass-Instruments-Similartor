"""Render the same trombone performance with a shorter, accented pressure rise."""
import argparse
import os
from pathlib import Path
import time

from trombone_axisymmetric import (Device, OUT as REFERENCE, PRESET, analyze,
                                   radiation, synthesize)
from trumpet import ROOT, json, np
from scipy.io import wavfile

OUT = ROOT/'output/trombone_attack'
ATTACK_PRESET = ROOT/'data/performance_trombone_attack.json'
WAV_NAME = 'trombone_Bb2_Bb4_major_f_strong_attack.wav'
AUDIO_RATE = 48000

# Rise time, temporary excess pressure, time for the excess to subside.
ARTICULATIONS = [(.012,.22,.024),(.008,.27,.022),(.017,.16,.027),(.023,.13,.020)]


def onset(waveform, late_rms):
    width=round(.016*AUDIO_RATE)
    cumulative=np.r_[0.,np.cumsum(waveform[:round(.22*AUDIO_RATE)]**2)]
    rms=np.sqrt(np.maximum(0.,(cumulative[width:]-cumulative[:-width])/width))
    continuous=np.convolve(np.asarray(rms>=.5*late_rms,'f8'),np.ones(width),'valid')
    indices=np.flatnonzero(continuous>=width)
    return (indices[0]+width/2)/AUDIO_RATE if len(indices) else .22


def early_rms(waveform):
    return float(np.std(waveform[round(.012*AUDIO_RATE):round(.080*AUDIO_RATE)]))


def render(destination=OUT, preset=None):
    os.makedirs(destination,exist_ok=True)
    reference_report=json.loads((REFERENCE/'measurements.json').read_text())
    ref_notes={r['note']:r for r in reference_report['notes']}
    notes=json.loads((PRESET if preset is None else preset).read_text())
    if [n['note'] for n in notes]!=list(ref_notes):
        raise ValueError('The preset must describe the same ordered scale as the reference.')
    device=Device(); waves=[]; report=[]; chosen=[]
    for base in notes:
        name=base['note'];reference=ref_notes[name]
        duration=reference['duration_s'];start,end=reference['window_s']
        with np.load(REFERENCE/f'{name}_coupled.npz') as saved:
            normal=radiation(saved['trace'])
        normal_attack=onset(normal,reference['radiated_rms_pa'])
        normal_early=early_rms(normal)
        if preset is None:
            alternatives=[dict(base,attack_s=attack,attack_boost=boost,attack_decay_s=decay)
                          for attack,boost,decay in ARTICULATIONS]
        else:
            alternatives=[base]
        best=None
        for note in alternatives:
            begin=time.perf_counter()
            try:
                trace,diag=synthesize(device,note,duration)
            except RuntimeError as exc:
                print(name,note['attack_s'],str(exc),flush=True)
                continue
            metrics=analyze(note,trace,start,end)
            wave=radiation(trace)
            first=onset(wave,metrics['radiated_rms_pa'])
            rms=early_rms(wave)
            gain_db=float(20*np.log10(max(rms,1e-12)/max(normal_early,1e-12)))
            earlier=1000*(normal_attack-first)
            metrics.update(attack_s=note['attack_s'],attack_boost=note['attack_boost'],
                attack_decay_s=note['attack_decay_s'],reference_onset_ms=normal_attack*1000,
                onset_ms=first*1000,onset_advance_ms=earlier,early_level_delta_db=gain_db,
                compute_s=time.perf_counter()-begin,fallback_count=diag['fallback_count'])
            print(json.dumps(metrics),flush=True)
            if (metrics['cents'] is None or abs(metrics['cents'])>14 or metrics['periodicity']<.92
                or abs(metrics['cents']-reference['cents'])>10
                or metrics['radiated_rms_pa']<.60*reference['radiated_rms_pa']):
                continue
            score=gain_db+earlier/5.
            if best is None or score>best[0]:
                best=(score,note,wave,metrics)
            if gain_db>=1.8 and earlier>=3.:
                break
        if best is None:
            raise RuntimeError(f'{name} did not settle on the intended pitch with these attacks.')
        _,note,wave,metrics=best
        chosen.append(note);report.append(metrics);waves.extend([wave,np.zeros(2880)])
        # Save progress for review without replacing the reference version.
        (destination/'measurements_in_progress.json').write_text(json.dumps(report,indent=2))
    audio=np.concatenate(waves)
    peak=float(np.max(abs(audio)))
    gain=min(reference_report['gain_fs_per_pa'],.95/peak)
    pcm=np.round(audio*gain*32767).astype('i2')
    wavfile.write(destination/WAV_NAME,AUDIO_RATE,pcm)
    result=dict(gpu=device.name,notes=report,duration_s=len(pcm)/AUDIO_RATE,
        gain_fs_per_pa=gain,reference_gain_fs_per_pa=reference_report['gain_fs_per_pa'],
        max_pitch_error_cents=max(abs(r['cents']) for r in report),
        median_onset_advance_ms=float(np.median([r['onset_advance_ms'] for r in report])),
        median_early_level_delta_db=float(np.median([r['early_level_delta_db'] for r in report])),
        clipped_samples=int(np.count_nonzero(abs(pcm.astype('i4'))>=32767)))
    (destination/'measurements.json').write_text(json.dumps(result,indent=2))
    if preset is None:
        ATTACK_PRESET.write_text(json.dumps(chosen,indent=2))
    print('Verified audio',len(pcm),'samples',result['median_onset_advance_ms'],
          'ms earlier (median)',result['median_early_level_delta_db'],'dB early level increase (median)',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preset',type=Path)
    args=parser.parse_args();render(preset=args.preset)
