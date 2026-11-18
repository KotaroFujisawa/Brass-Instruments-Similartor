"""Render the Bb3-Bb5 scale with prescribed-source nonlinear bore propagation."""
from nonlinear_propagation import *
from scipy.io import wavfile
import argparse
import hashlib
import os


def spectral_metrics(y,fs,start=.28,end=.56):
    v=y[int(start*fs):int(end*fs)]
    p=abs(rfft((v-v.mean())*np.hanning(len(v))))**2
    f=rfftfreq(len(v),1/fs)
    total=p[(f>150)&(f<20000)].sum()
    return dict(hf_fraction_above_3khz=float(p[(f>3000)&(f<20000)].sum()/total),
        centroid_hz=float(np.sum(f*p)/p.sum()),rms_pa=float(np.std(v)))


def run(fs=192000,dx=.01,only=None):
    outdir=ROOT/'output/nonlinear';os.makedirs(outdir,exist_ok=True)
    (ROOT/'tmp').mkdir(exist_ok=True)
    settings=json.loads((ROOT/'data/performance.json').read_text())
    depth=json.loads((ROOT/'data/model.json').read_text())['cup_depth_m']
    cache={};waves=[];refs=[];reports=[];offset=0.
    for i,s in enumerate(settings):
        if only and s['note']!=only:continue
        valves=tuple(s['valves_semitones'])
        if valves not in cache:
            geo=propagation_geometry(valves,depth)
            key=hashlib.sha256((ROOT/'src/trumpet.py').read_bytes()
                +(ROOT/'src/nonlinear_propagation.py').read_bytes()
                +(ROOT/'data/model.json').read_bytes()
                +(ROOT/'data/york_trumpet_fig7_52.csv').read_bytes()
                +str((fs,valves)).encode()).hexdigest()[:16]
            file=ROOT/'tmp'/f'incident_{key}.npz'
            if file.exists():impulse=np.load(file)['impulse']
            else:
                impulse=incident_transfer(geo,fs=fs,nfft=int(NFFT*fs/FS))
                np.savez(file,impulse=impulse)
            # Use the existing response, or compute it on a fresh checkout.
            res=resonator(valves,depth)
            cache[valves]=(geo,impulse,res)
        geo,impulse,res=cache[valves]
        duration=1.1 if i==len(settings)-1 else .72
        tr=simulate(res['r'],res['z0'],FS,duration,s['lip_frequency_hz'],
                    s['blowing_pressure_pa'],s['rest_y_m'],hold=duration-.14)
        up=fs//FS
        assert fs%FS==0
        flow=signal.resample_poly(tr[:,1],up,1)
        drive=signal.fftconvolve(flow,impulse)[:len(flow)]
        nl,info,snapshots=propagate(drive,fs,geo,beta=BETA,dx=dx,record=s['note']=='Bb4')
        lin,_,_=propagate(drive,fs,geo,beta=0,dx=dx)
        yn=radiate_incident(nl,fs,geo);yl=radiate_incident(lin,fs,geo)
        dec=fs//48000
        yn=signal.resample_poly(yn,1,dec);yl=signal.resample_poly(yl,1,dec)
        m=measure(np.c_[yn,np.zeros((len(yn),3))],fs=48000,start=.28,end=duration-.15)
        metrics_nl=spectral_metrics(yn,48000,end=duration-.15)
        metrics_lin=spectral_metrics(yl,48000,end=duration-.15)
        entry=dict(note=s['note'],target_hz=s['target_hz'],frequency_hz=m['frequency'],
            cents=float(1200*np.log2(m['frequency']/s['target_hz'])),periodicity=m['periodicity'],
            start_s=offset,duration_s=duration,propagation=info,
            nonlinear=metrics_nl,linear=metrics_lin,
            hf_fraction_ratio_db=float(10*np.log10(metrics_nl['hf_fraction_above_3khz']/metrics_lin['hf_fraction_above_3khz'])))
        reports.append(entry)
        waves.extend([yn,np.zeros(2880)]);refs.extend([yl,np.zeros(2880)])
        offset+=duration+.06
        np.savez_compressed(outdir/f'{s["note"]}_fields.npz',drive_pa=drive,
            nonlinear_bell_incident_pa=nl,linear_bell_incident_pa=lin,fs=fs,
            nonlinear_radiated_pa=yn,linear_radiated_pa=yl,audio_fs=48000,
            snapshot_x=np.array([x for x,_ in snapshots]),
            snapshot_pressure=np.array([p for _,p in snapshots]))
        print(json.dumps(entry),flush=True)
    nl=np.concatenate(waves);lin=np.concatenate(refs)
    gain=.89/max(np.max(abs(nl)),np.max(abs(lin)))
    suffix=f'_{only}' if only else ''
    wavfile.write(outdir/f'trumpet_Bb3_Bb5_major_f_nonlinear{suffix}.wav',48000,np.round(nl*gain*32767).astype(np.int16))
    wavfile.write(outdir/f'trumpet_Bb3_Bb5_major_f_linear_reference{suffix}.wav',48000,np.round(lin*gain*32767).astype(np.int16))
    if not only:
        ab=np.r_[lin,np.zeros(24000),nl]
        wavfile.write(outdir/'comparison_linear_then_nonlinear.wav',48000,np.round(ab*gain*32767).astype(np.int16))
    report=dict(model='prescribed-incident-wave generalised Burgers, one-way',
        nonlinear_feedback_to_lips=False,propagation_fs=fs,audio_fs=48000,
        beta=BETA,sound_diffusivity_m2_s=SOUND_DIFFUSIVITY,output_gain_fs_per_pa=float(gain),
        common_gain_for_linear_and_nonlinear=True,duration_s=len(nl)/48000,notes=reports)
    (outdir/f'measurements{suffix}.json').write_text(json.dumps(report,indent=2))
    print('Completed',len(reports),'notes;',len(nl)/48000,'s',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--only');p.add_argument('--fs',type=int,default=192000)
    p.add_argument('--dx',type=float,default=.01);args=p.parse_args()
    run(args.fs,args.dx,args.only)
