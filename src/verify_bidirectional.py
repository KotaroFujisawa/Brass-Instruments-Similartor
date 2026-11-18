"""Checks for the causal, genuinely closed-loop bidirectional solver."""
from bidirectional import *
from scipy.special import jv
from scipy.io import wavfile
import argparse
import os


def uniform_cfg(length=.7,radius=.006,fs=192000,losses=False):
    cfg=setup(fs=fs,losses=losses)
    n=int(round(length*fs/C))
    rad=np.full(n,radius);cfg.update(radius=rad,area=np.pi*rad**2,reflection=np.zeros(n-1),length=n*C/fs)
    cfg['loss_k']=np.zeros(n);cfg['loss_den']=np.ones(n)
    cfg['z0']=RHO*C/(np.pi*radius**2);cfg['zbell']=cfg['z0']
    return cfg


def solver_tests():
    cfg=uniform_cfg();fs=cfg['fs'];sec=.075
    t=np.arange(int(sec*fs))/fs;freq=600.;amp=5000.
    drive=amp*np.sin(2*np.pi*freq*t)
    tr,_=simulate_coupled(cfg,sec,500,0,0,source=drive,matched_end=True,mean_hz=0.)
    y=tr[int(.025*fs):int(.065*fs),11]
    harmonics=[];ts=t[int(.025*fs):int(.065*fs)]
    for h in range(1,5):harmonics.append(2*abs(np.sum(y*np.exp(-2j*np.pi*h*freq*ts)))/len(y)/amp)
    sigma=1.2*2*np.pi*freq*amp*cfg['length']/(RHO*C**3)
    h=np.arange(1,5);exact=2*jv(h,h*sigma)/(h*sigma)
    # In the time-domain speed model this retarded-time solution is asymptotic.
    fubini=dict(sigma=sigma,measured=harmonics,asymptotic=exact.tolist(),max_error=float(max(abs(np.array(harmonics)-exact))))
    # A compact pulse must travel to the bell and back before reaching the lips.
    sec=.012;t=np.arange(int(sec*fs))/fs
    pulse=np.exp(-.5*((t-.0005)/.000035)**2)
    echo,_=simulate_coupled(cfg,sec,500,0,0,beta=0,source=pulse)
    returned=echo[:,8];expected=.0005+2*cfg['length']/C
    peak_time=np.argmax(abs(returned))/fs
    pre=returned[t<expected-.0002]
    propagation=dict(expected_return_s=expected,peak_return_s=peak_time,
        early_return_max=float(max(abs(pre))),return_peak=float(max(abs(returned))))
    # Linear limit: actual time-domain reflections vs analytical discrete TMM.
    cfg=setup();sec=.12;fs=cfg['fs'];t=np.arange(int(sec*fs))/fs
    pulse=np.exp(-.5*((t-.001)/.000025)**2)
    lin,_=simulate_coupled(cfg,sec,500,0,0,beta=0,source=pulse)
    nfft=65536;freqs=np.fft.rfftfreq(nfft,1/fs)
    r_measured=np.fft.rfft(lin[:,8],nfft)/(np.fft.rfft(pulse,nfft)+1e-30)
    sel=(freqs>=150)&(freqs<=2000)
    z=impedance(cfg,freqs[sel]);r_exact=(z-cfg['z0'])/(z+cfg['z0'])
    linear=dict(reflection_rmse=float(np.sqrt(np.mean(abs(r_measured[sel]-r_exact)**2))),
        max_reflection_150_2000=float(max(abs(r_measured[sel]))),
        wall_fractional_fit_relative_error=cfg['loss_fit_max_relative_error'])
    report=dict(fubini=fubini,round_trip=propagation,linear_limit=linear)
    report['passed']=(fubini['max_error']<.015 and propagation['early_return_max']<1e-5
        and abs(peak_time-expected)<.0003 and linear['reflection_rmse']<.005
        and linear['max_reflection_150_2000']<=1.001)
    folder=ROOT/'output/bidirectional';os.makedirs(folder,exist_ok=True)
    (folder/'solver_tests.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    if not report['passed']:raise SystemExit('Solver verification failed')


def full_tests():
    from render_bidirectional import spectrum_metrics
    folder=ROOT/'output/bidirectional'
    settings=json.loads((ROOT/'data/performance_bidirectional.json').read_text())
    measurements=json.loads((folder/'measurements.json').read_text())
    boundaries=[];convergence=[]
    for s in settings:
        data=np.load(folder/f'{s["note"]}_coupled.npz');tr=data['nonlinear'];fs=int(data['fs'])
        cfg=setup(s['valves_semitones'],fs=fs)
        residual=tr[:,0]-cfg['z0']*tr[:,1]-2*tr[:,8]
        area=2*B*tr[:,3];sel=area>1e-12
        aa=RHO*(.5/area[sel]**2-1/(SCUP*area[sel])+1/SCUP**2)
        previous=np.r_[0,tr[:-1,4]]
        flow=aa*tr[sel,4]*abs(tr[sel,4])+RHO*D/area[sel]*fs*(tr[sel,4]-previous[sel])-(tr[sel,6]-tr[sel,0])
        boundaries.append(dict(note=s['note'],boundary_residual_max_pa=float(max(abs(residual))),
            flow_equation_residual_max_pa=float(max(abs(flow))),finite=bool(np.isfinite(tr).all())))
        if s['note'] in ['Bb4','G4','Bb5']:
            duration=len(tr)/fs;end=duration-.15
            finecfg=setup(s['valves_semitones'],fs=2*fs)
            fine,mach=simulate_coupled(finecfg,duration,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'])
            yn=data['nonlinear_radiated_pa'];yf=radiated(fine,2*fs)
            coarse_m=measure(tr,fs=fs,start=.28,end=end)
            fine_m=measure(fine,fs=2*fs,start=.28,end=end)
            c=spectrum_metrics(yn,end=end);f=spectrum_metrics(yf,end=end)
            row=dict(note=s['note'],coarse_fs=fs,fine_fs=2*fs,
                coarse_frequency_hz=coarse_m['frequency'],fine_frequency_hz=fine_m['frequency'],
                pitch_difference_cents=float(1200*np.log2(fine_m['frequency']/coarse_m['frequency'])),
                radiated_rms_difference_db=float(20*np.log10(f['rms_pa']/c['rms_pa'])),
                hf_fraction_difference_db=float(10*np.log10(f['hf_fraction_above_3khz']/c['hf_fraction_above_3khz'])))
            convergence.append(row);print('Convergence',row,flush=True)
            np.savez_compressed(folder/f'{s["note"]}_fine.npz',trace=fine[::8],fs=2*fs/8,radiated_pa=yf,audio_fs=48000)
    # Feedback ablation changes only the incoming-wave term in the lip boundary.
    s=next(s for s in settings if s['note']=='Bb4');cfg=setup()
    ablation,_=simulate_coupled(cfg,.72,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'],feedback=False)
    data=np.load(folder/'Bb4_coupled.npz');tr=data['nonlinear'];lin=data['linear'];fs=int(data['fs'])
    section=slice(int(.28*fs),int(.57*fs))
    diagnostic=dict(note='Bb4',ablation=measure(ablation,fs=fs,start=.28,end=.57),
        coupled=measure(tr,fs=fs,start=.28,end=.57),
        beta_zero=measure(lin,fs=fs,start=.28,end=.57),
        lip_y_change_rms_m=float(np.std(tr[section,3]-lin[section,3])),
        incoming_change_rms_pa=float(np.std(tr[section,8]-lin[section,8])),
        note_on_difference_metrics='Unaligned trajectories include the genuine pitch/phase difference.')
    fs,pcm=wavfile.read(folder/'trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear.wav')
    audio=dict(fs=fs,duration_s=len(pcm)/fs,clipped_samples=int(np.sum(abs(pcm.astype(float))>=32767)),
        max_pitch_error_cents=float(max(abs(s['cents']) for s in measurements['notes'])),
        minimum_periodicity=float(min(s['periodicity'] for s in measurements['notes'])))
    report=dict(audio=audio,boundaries=boundaries,convergence=convergence,feedback_diagnostic=diagnostic,
                subharmonic_diagnostic=subharmonic_diagnostic(folder),
                convergence_frequency_definition='Mouthpiece-pressure autocorrelation period; G4 has weak period doubling. Audio pitches are measured separately.')
    report['passed']=(audio['clipped_samples']==0 and audio['max_pitch_error_cents']<6
        and audio['minimum_periodicity']>.98
        and all(b['finite'] and b['boundary_residual_max_pa']<1e-7 and b['flow_equation_residual_max_pa']<1e-6 for b in boundaries)
        and all(abs(c['pitch_difference_cents'])<10 and abs(c['radiated_rms_difference_db'])<1 for c in convergence)
        and diagnostic['lip_y_change_rms_m']>1e-6 and diagnostic['incoming_change_rms_pa']>100)
    (folder/'verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    if not report['passed']:raise SystemExit('Full verification needs review')


def subharmonic_diagnostic(folder):
    data=np.load(folder/'G4_coupled.npz');result={}
    for name,y,fs in [('mouthpiece',data['nonlinear'][:,0],int(data['fs'])),
                      ('radiated',data['nonlinear_radiated_pa'],48000)]:
        a=y[int(.35*fs):int(.57*fs)];t=np.arange(len(a))/fs;a=a*np.hanning(len(a))
        f=392.17
        amps=[abs(np.sum(a*np.exp(-2j*np.pi*h*f*t))) for h in [.5,1,1.5,2,2.5,3]]
        result[name+'_levels_db']=[float(20*np.log10(v/amps[1])) for v in amps]
    result.update(note='G4',reference_hz=392.17,multipliers=[.5,1,1.5,2,2.5,3],window_s=[.35,.57],
        interpretation='Weak period doubling: a radiated component near 196 Hz, about 32 dB below the 392 Hz carrier. No spectral component has been removed.')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');args=parser.parse_args()
    if args.full:full_tests()
    else:solver_tests()
