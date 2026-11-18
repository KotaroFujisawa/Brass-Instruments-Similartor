"""Analytical and numerical validation of the generalised Burgers extension."""
from nonlinear_propagation import *
from scipy.special import jv
from scipy.io import wavfile
from render_nonlinear import spectral_metrics
import argparse
import os


def analytical_checks():
    n=4096;theta=2*np.pi*np.arange(n)/n
    inp=np.sin(theta)
    sigma=.6
    out,steps=burgers_step(inp,1.,sigma,2*np.pi/n)
    actual=2*abs(rfft(out))[1:7]/n
    h=np.arange(1,7);expected=2*jv(h,h*sigma)/(h*sigma)
    errors=abs(actual-expected)
    shock,_=burgers_step(inp,1.,2.,2*np.pi/n)
    shock2,_=burgers_step(np.sin(2*np.pi*np.arange(2*n)/(2*n)),1.,2.,np.pi/n)
    spectrum1=2*abs(rfft(shock))[1:17]/n
    spectrum2=abs(rfft(shock2))[1:17]/n
    return dict(fubini=dict(sigma=sigma,harmonics=list(range(1,7)),
        measured=actual.tolist(),analytic=expected.tolist(),max_absolute_error=float(errors.max()),
        nonlinear_substeps=steps),
        shock=dict(sigma=2.,finite=bool(np.isfinite(shock).all()),
            peak=float(max(abs(shock))),input_mean_square=float(np.mean(inp**2)),
            output_mean_square=float(np.mean(shock**2)),mean=float(np.mean(shock)),
            grid_refinement_harmonic_rmse=float(np.sqrt(np.mean((spectrum1-spectrum2)**2)))))


def numerical_checks():
    fs=96000;t=np.arange(int(.20*fs))/fs
    inp=300*np.sin(2*np.pi*500*t)*np.sin(np.pi*t/t[-1])**2
    geo=dict(length=.8,body_x=np.array([0.,.8]),body_radii=np.array([.005,.015]),
             inlet_diameter=.01,outlet_diameter=.03)
    # Zero-nonlinearity, no losses: exact spherical/area spreading plus flight time.
    exact,_,_=propagate(inp,fs,geo,beta=0,wall=False,bulk=False)
    delay=.8/C;inside=(t>.02)&(t<.17)
    # Fractional delay interpolation is an independent lower-order reference.
    reference=np.interp(t-delay,t,inp,left=0,right=0)/3
    spread_error=np.linalg.norm(exact[inside]-reference[inside])/np.linalg.norm(reference[inside])
    weak,_,_=propagate(inp*1e-4,fs,geo,beta=BETA)
    linear,_,_=propagate(inp*1e-4,fs,geo,beta=0)
    weak_error=np.linalg.norm(weak-linear)/np.linalg.norm(linear)
    # Wall loss check against the analytical leading sqrt(f) attenuation ratio.
    f1=500.;f2=2000.;a=.006;length=1.
    attenuation=lambda f:np.exp(-2*np.pi*1.045/(C*632.8*a)*np.sqrt(f)*length)
    return dict(lossless_flare_relative_error=float(spread_error),
                small_signal_nonlinear_vs_linear_relative_error=float(weak_error),
                wall_amplitude_500hz=float(attenuation(f1)),wall_amplitude_2000hz=float(attenuation(f2)))


def convergence_check(note='Bb4'):
    folder=ROOT/'output/nonlinear'
    field=np.load(folder/f'{note}_fields.npz');fs=int(field['fs'])
    model=json.loads((ROOT/'data/model.json').read_text())
    setting=next(x for x in json.loads((ROOT/'data/performance.json').read_text()) if x['note']==note)
    geo=propagation_geometry(setting['valves_semitones'],model['cup_depth_m'])
    drive=field['drive_pa']
    refined=signal.resample_poly(drive,2,1)
    fine,_,_=propagate(refined,2*fs,geo,beta=BETA,dx=.005)
    y=signal.resample_poly(radiate_incident(fine,2*fs,geo),1,2*fs//48000)
    old=field['nonlinear_radiated_pa']
    seg=slice(int(.3*48000),int(.55*48000))
    # Report both time-domain and spectral differences; shock fronts may be
    # sharper at the finer mesh even when the audible harmonic envelope converges.
    aa=old[seg];bb=y[seg]
    f0=setting['target_hz']
    def harmonics(v):
        ft=rfftfreq(131072,1/48000);sp=abs(rfft(v*np.hanning(len(v)),131072))
        return np.array([sp[abs(ft-k*f0)<12].max() for k in range(1,31)])
    ha=harmonics(aa);hb=harmonics(bb)
    active=(hb>hb.max()*.003)&(np.arange(1,31)*f0<12000)
    return dict(note=note,base_fs=fs,refined_fs=2*fs,refined_dx_m=.005,
        audible_waveform_relative_l2=float(np.linalg.norm(aa-bb)/np.linalg.norm(bb)),
        harmonic_db_rmse_below_12khz=float(np.sqrt(np.mean((20*np.log10(ha[active]/hb[active]))**2))),
        base=spectral_metrics(old,48000),refined=spectral_metrics(y,48000))


def main(full=False):
    report=dict(analytical=analytical_checks(),numerical=numerical_checks())
    if full:
        report['convergence']=convergence_check()
        notes=json.loads((ROOT/'output/nonlinear/measurements.json').read_text())['notes']
        strongest=max(notes,key=lambda x:x['nonlinear']['hf_fraction_above_3khz'])['note']
        report['convergence_strongest']=convergence_check(strongest)
        rate,audio=wavfile.read(ROOT/'output/nonlinear/trumpet_Bb3_Bb5_major_f_nonlinear.wav')
        report['audio']=dict(rate=int(rate),seconds=len(audio)/rate,finite=bool(np.isfinite(audio).all()),
            clipped_samples=int(np.sum(abs(audio.astype(float))>=32767)),
            max_pitch_error_cents=max(abs(n['cents']) for n in notes),
            max_input_mach=max(n['propagation']['peak_input_mach'] for n in notes),
            min_periodicity=min(n['periodicity'] for n in notes))
    report['passed']=(report['analytical']['fubini']['max_absolute_error']<.002
        and report['analytical']['shock']['output_mean_square']<.5
        and report['analytical']['shock']['peak']<=1.001
        and report['analytical']['shock']['grid_refinement_harmonic_rmse']<.005
        and report['numerical']['lossless_flare_relative_error']<.001
        and report['numerical']['small_signal_nonlinear_vs_linear_relative_error']<.001)
    if full:report['passed']=report['passed'] and report['audio']['finite'] and not report['audio']['clipped_samples'] and report['audio']['max_pitch_error_cents']<6 and report['convergence']['harmonic_db_rmse_below_12khz']<2 and report['convergence_strongest']['harmonic_db_rmse_below_12khz']<2
    dest=ROOT/'output/nonlinear';os.makedirs(dest,exist_ok=True)
    (dest/('verification.json' if full else 'solver_tests.json')).write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    if not report['passed']:raise SystemExit('Verification failed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');args=p.parse_args()
    main(args.full)
