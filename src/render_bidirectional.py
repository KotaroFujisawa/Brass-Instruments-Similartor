"""Self-excited two-way nonlinear Bb-major scale, with beta=0 control."""
from bidirectional import *
from scipy.io import wavfile
import argparse
import hashlib
import os


def measure_audio(y,start,end,*,min_hz=150,max_hz=1200):
    # Sharp nonlinear waves have a narrow autocorrelation maximum. Interpolate
    # before selecting its first peak to avoid integer-lag octave errors.
    up=signal.resample_poly(y,4,1)
    return measure(np.c_[up,np.zeros((len(up),3))],fs=192000,start=start,end=end,
                   min_hz=min_hz,max_hz=max_hz)


def spectrum_metrics(y,fs=48000,start=.28,end=.57):
    y=y[int(start*fs):int(end*fs)];y=y-y.mean()
    power=abs(np.fft.rfft(y*np.hanning(len(y))))**2
    freq=np.fft.rfftfreq(len(y),1/fs)
    band=(freq>150)&(freq<20000)
    return dict(rms_pa=float(np.std(y)),centroid_hz=float(np.sum(freq[band]*power[band])/power[band].sum()),
                hf_fraction_above_3khz=float(power[(freq>3000)&(freq<20000)].sum()/power[band].sum()))


def run(fs=192000):
    folder=ROOT/'output/bidirectional';os.makedirs(folder,exist_ok=True)
    settings=json.loads((ROOT/'data/performance_bidirectional.json').read_text())
    assert len(settings)==15
    cache={};waves=[];refs=[];reports=[];offset=0.
    for i,s in enumerate(settings):
        valves=tuple(s['valves_semitones'])
        if valves not in cache:cache[valves]=setup(valves,fs=fs)
        cfg=cache[valves];duration=1.1 if i==14 else .72
        signature=hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in [
            'src/bidirectional.py','src/trumpet.py','src/nonlinear_propagation.py',
            'data/model.json','data/york_trumpet_fig7_52.csv'])
            +json.dumps([s,fs,duration],sort_keys=True).encode()).hexdigest()
        tracefile=folder/f'{s["note"]}_coupled.npz'
        saved=dict(np.load(tracefile)) if tracefile.exists() else {}
        if str(saved.get('signature',''))==signature:
            tr=saved['nonlinear'];lin=saved['linear'];mach=float(saved['max_mach'])
            yn=saved['nonlinear_radiated_pa'];yl=saved['linear_radiated_pa']
        else:
            tr,mach=simulate_coupled(cfg,duration,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'])
            lin,mach_lin=simulate_coupled(cfg,duration,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'],beta=0.)
            yn=radiated(tr,fs);yl=radiated(lin,fs)
            np.savez_compressed(tracefile,nonlinear=tr,linear=lin,fs=fs,
                nonlinear_radiated_pa=yn,linear_radiated_pa=yl,audio_fs=48000,
                signature=signature,max_mach=mach)
        end=duration-.15
        mn=measure_audio(yn,start=.28,end=end)
        ml=measure_audio(yl,start=.28,end=end)
        section=slice(int(.28*fs),int(end*fs))
        boundary=np.max(abs(tr[:,0]-cfg['z0']*tr[:,1]-2*tr[:,8]))
        entry=dict(note=s['note'],target_hz=s['target_hz'],frequency_hz=mn['frequency'],
            cents=float(1200*np.log2(mn['frequency']/s['target_hz'])),periodicity=mn['periodicity'],
            linear_frequency_hz=ml['frequency'],linear_periodicity=ml['periodicity'],
            nonlinear_pitch_change_cents=float(1200*np.log2(mn['frequency']/ml['frequency'])),
            nonlinear=spectrum_metrics(yn,end=end),linear=spectrum_metrics(yl,end=end),
            start_s=offset,duration_s=duration,blowing_pressure_pa=s['blowing_pressure_pa'],
            lip_frequency_hz=s['lip_frequency_hz'],max_wave_mach_sampled=mach,
            boundary_residual_max_pa=float(boundary),
            incoming_rms_pa=float(np.std(tr[section,8])),linear_incoming_rms_pa=float(np.std(lin[section,8])),
            lip_y_rms_m=float(np.std(tr[section,3])),linear_lip_y_rms_m=float(np.std(lin[section,3])))
        reports.append(entry)
        # Full-rate boundary histories support reproducible feedback diagnostics.
        waves.extend([yn,np.zeros(2880)]);refs.extend([yl,np.zeros(2880)]);offset+=duration+.06
        print(json.dumps(entry),flush=True)
    yn=np.concatenate(waves);yl=np.concatenate(refs)
    gain=.89/max(np.max(abs(yn)),np.max(abs(yl)))
    for name,y in [('trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear',yn),
                   ('trumpet_Bb3_Bb5_major_f_bidirectional_linear_reference',yl),
                   ('comparison_bidirectional_linear_then_nonlinear',np.r_[yl,np.zeros(24000),yn])]:
        wavfile.write(folder/f'{name}.wav',48000,np.round(y*gain*32767).astype(np.int16))
    wavfile.write(folder/'radiation_bidirectional_nonlinear_Pa_float.wav',48000,yn.astype(np.float32))
    report=dict(model='two-way nonlinear waveguide + Adachi-Sato 2D lips',
        nonlinear_feedback_to_lips=True,beta=1.2,propagation_fs=fs,audio_fs=48000,
        common_gain_for_linear_and_nonlinear=True,output_gain_fs_per_pa=float(gain),
        duration_s=len(yn)/48000,notes=reports)
    (folder/'measurements.json').write_text(json.dumps(report,indent=2))
    print('Completed',len(reports),'notes;',len(yn)/48000,'s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fs',type=int,default=192000)
    run(parser.parse_args().fs)
