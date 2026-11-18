"""Physical/numerical checks and a standalone diagnostic figure."""
from trumpet import *
from scipy.io import wavfile


def main():
    report=json.loads((ROOT/'output/measurements.json').read_text())
    perf=json.loads((ROOT/'data/performance.json').read_text())
    fs,audio=wavfile.read(ROOT/'output/trumpet_Bb3_Bb5_major_f.wav')
    checks=dict(sample_rate=int(fs),samples=len(audio),finite=bool(np.isfinite(audio).all()),
                pcm_peak=int(np.max(abs(audio.astype(float)))),clipped_samples=int(np.sum(abs(audio.astype(float))>=32767)))
    cents=np.array([n['cents'] for n in report['notes']])
    checks['max_pitch_error_cents']=float(np.max(abs(cents)))
    checks['minimum_periodicity']=float(min(n['periodicity'] for n in report['notes']))
    # Reflection filter passivity at and between FFT bins; errors relative to TMM.
    acoustic=[]
    for valves in set(tuple(p['valves_semitones']) for p in perf):
        res=dict(np.load(ROOT/'tmp'/('res_'+''.join(map(str,valves))+'.npz')))
        r=np.fft.rfft(res['r'],4*NFFT)
        mask=(res['freq']>150)&(res['freq']<2000)
        actual=res['actual_reflection']
        z= float(res['z0'])*(1+actual)/(1-actual+1e-30)
        # Reflection error is meaningful even at impedance zeros.
        acoustic.append(dict(valves=list(valves),max_reflection=float(max(abs(r))),
            minimum_tmm_resistance=float(min(res['zin'].real)),
            reflection_rmse_150_2000=float(np.sqrt(np.mean(abs(actual[mask]-res['reflection'][mask])**2))),
            acausal_energy_ratio=float(res['acausal_energy_ratio'])))
    checks['acoustic']=acoustic
    # The nonlinear flow boundary residual is checked on a full-rate new run.
    s=perf[7]
    res=dict(np.load(ROOT/'tmp/res_.npz'))
    tr=simulate(res['r'],float(res['z0']),FS,.48,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'],hold=.48)
    area=np.maximum(2*B*tr[1:,3],0)
    sel=area>1e-10
    sa=area[sel];q=tr[1:,4][sel];qprev=tr[:-1,4][sel]
    a=RHO*(.5/sa**2-1/(SCUP*sa)+1/SCUP**2)
    flow_residual=a*q*abs(q)+RHO*D/sa*(q-qprev)*FS-(tr[1:,6][sel]-tr[1:,0][sel])
    checks['flow_equation_max_residual_pa']=float(max(abs(flow_residual)))
    checks['lip_swept_flow_rms_m3_s']=float(np.std(tr[:,5]))
    # Halve dt, regenerate the same bore response and compare the selected Bb4.
    model=json.loads((ROOT/'data/model.json').read_text())
    res2=resonator((),model['cup_depth_m'],fs=2*FS,nfft=2*NFFT)
    tr2=simulate(res2['r'],res2['z0'],2*FS,.48,s['lip_frequency_hz'],s['blowing_pressure_pa'],s['rest_y_m'],hold=.48)
    m1=measure(tr,start=.30);m2=measure(tr2,fs=2*FS,start=.30)
    checks['time_step_check_Bb4']=dict(fs_96000=m1,fs_192000=m2,
        pitch_difference_cents=float(1200*np.log2(m2['frequency']/m1['frequency'])),
        pressure_rms_relative_change=float(m2['rms']/m1['rms']-1))
    checks['time_step_checks_extremes']=[]
    for setting in [perf[0],perf[-1]]:
        dur=.6 if setting['note']=='Bb3' else 1.05
        traces=[]
        for rr,rate in [(res,FS),(res2,2*FS)]:
            test=simulate(rr['r'],float(rr['z0']),rate,dur,setting['lip_frequency_hz'],
                          setting['blowing_pressure_pa'],setting['rest_y_m'],hold=dur)
            traces.append(measure(test,fs=rate,start=dur-.2))
        checks['time_step_checks_extremes'].append(dict(note=setting['note'],
            pitch_difference_cents=float(1200*np.log2(traces[1]['frequency']/traces[0]['frequency'])),
            rms_relative_change=float(traces[1]['rms']/traces[0]['rms']-1)))
    # Smaller acoustic elements check a representative impedance peak.
    lengths,rad,_,_=geometry((),model['cup_depth_m'],dx=.001)
    freqs=np.linspace(430,500,701)
    zfine,_=tmm(freqs,lengths,rad)
    lengths,rad,_,_=geometry((),model['cup_depth_m'],dx=.002)
    zcoarse,_=tmm(freqs,lengths,rad)
    checks['space_step_check_Bb4']=dict(peak_1mm_hz=float(freqs[np.argmax(abs(zfine))]),
                                       peak_2mm_hz=float(freqs[np.argmax(abs(zcoarse))]))
    checks['passed']=(checks['finite'] and not checks['clipped_samples']
        and checks['max_pitch_error_cents']<5 and checks['minimum_periodicity']>.97
        and checks['flow_equation_max_residual_pa']<1e-5
        and all(x['max_reflection']<=1.001 for x in acoustic)
        and abs(checks['time_step_check_Bb4']['pitch_difference_cents'])<10
        and all(abs(x['pitch_difference_cents'])<10 for x in checks['time_step_checks_extremes']))
    (ROOT/'output/verification.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2),flush=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(3,2,figsize=(13,11),constrained_layout=True)
    fig.suptitle('TRUMPET  /  Bb3-Bb5  /  FORTE\nAdachi-Sato 2D lips + reconstructed mouthpiece + digitised York bore',fontsize=16)
    x=res['x'];r=res['radius_points']
    axs[0,0].plot(x*1000,r*2000,color='#c07b29')
    axs[0,0].set(xlabel='Distance from lips (mm)',ylabel='Bore diameter (mm)',title='Geometry from The Science of Brass Instruments, Fig. 7.52')
    inset=axs[0,0].inset_axes([.1,.48,.44,.4]);inset.plot(x[:43]*1000,r[:43]*2000,color='#c07b29');inset.set_title('Cup / throat / backbore',fontsize=9)
    mask=res['freq']<2200
    z0=float(res['z0']);rr=res['actual_reflection']
    zz=z0*(1+rr)/(1-rr+1e-30)
    axs[0,1].plot(res['freq'][mask],abs(res['zin'][mask])/z0,label='Geometry TMM')
    axs[0,1].plot(res['freq'][mask],abs(zz[mask])/z0,'--',label='Causal simulation',alpha=.8)
    axs[0,1].set(xlabel='Frequency (Hz)',ylabel='|Zin| / Z0',title='Mouthpiece and bore input impedance');axs[0,1].legend(fontsize=8)
    tt=np.arange(len(audio))/fs
    axs[1,0].plot(tt,audio/32768,lw=.35,color='#245c8c')
    axs[1,0].set(xlabel='Time (s)',ylabel='Amplitude (FS)',title='Radiated sound - one common gain, no pitch processing')
    for note in report['notes']:
        axs[1,0].text(note['start_s']+.1,.94,note['note'],rotation=90,fontsize=8,va='top')
    axs[1,1].bar([n['note'] for n in report['notes']],cents,color='#245c8c')
    axs[1,1].axhline(0,color='black',lw=.6);axs[1,1].axhline(5,color='gray',ls=':');axs[1,1].axhline(-5,color='gray',ls=':')
    axs[1,1].tick_params(axis='x',rotation=50)
    axs[1,1].set(ylabel='Pitch error (cents)',title='Steady-state pitch vs A4 = 440 Hz')
    n0=int(.32*FS);n1=n0+int(.01*FS)
    axs[2,0].plot(tr[n0:n1,2]*1000,tr[n0:n1,3]*1000,color='#c07b29')
    axs[2,0].axhline(0,lw=.7,color='gray');axs[2,0].set(xlabel='Lip axial position x (mm)',ylabel='Lip opening coordinate y (mm)',title='Bb4: two-dimensional lip trajectory')
    trace=np.load(ROOT/'output/trace_Bb4.npz')['radiation_pa']
    sig=trace[int(.28*48000):int(.55*48000)]
    spec=abs(np.fft.rfft(sig*np.hanning(len(sig)),32768));freq=np.fft.rfftfreq(32768,1/48000)
    axs[2,1].plot(freq,20*np.log10(np.maximum(spec/spec.max(),1e-7)),color='#245c8c')
    axs[2,1].set(xlim=(0,8000),ylim=(-100,3),xlabel='Frequency (Hz)',ylabel='Relative level (dB)',title='Bb4: radiated harmonics')
    fig.savefig(ROOT/'output/diagnostics.png',dpi=150)
    plt.close(fig)
    if not checks['passed']:raise SystemExit('Validation failed; inspect verification.json')


if __name__=='__main__':main()
