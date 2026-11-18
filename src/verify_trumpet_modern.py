"""Check the new bell load, inlet coupling and open-fingering input impedance."""
import os

from render_trumpet_modern import MODEL, OUT, Device, Run, performance, simulate
from render_axisymmetric import TRACE_FS, analyze, trumpet_mesh
from thermoviscous import fit
from trumpet import C, RHO, json, np, signal


def impedance(device, mesh, model):
    amplitude = 30.
    with Run(device, mesh, .65, mode=3, pressure=amplitude, **model) as run:
        run.until()
        trace = run.history()
    area = np.pi*mesh['radius'][0]**2
    t = (np.arange(len(trace))+1)/TRACE_FS
    flow = area*amplitude/(RHO*C)*np.exp(-.5*((t-.00045)/.00007)**2)
    pressure = trace[:,8]+RHO*C/area*flow
    nfft = 2**19
    f = np.fft.rfftfreq(nfft, 1/TRACE_FS)
    z = np.fft.rfft(pressure, nfft)/np.fft.rfft(flow, nfft)
    band = (f>=40)&(f<=2400)
    return f[band], z[band], float(np.mean(pressure[-TRACE_FS//20:]))


def resonances(freq, impedance):
    mag = abs(impedance)
    indices, _ = signal.find_peaks(mag, prominence=2e6)
    intervals = freq[1]-freq[0]
    peaks = []
    for i in indices[:8]:
        a, b, c = mag[i-1:i+2]
        frac = .5*(a-c)/(a-2*b+c)
        peaks.append(dict(hz=float(freq[i]+frac*intervals), mpa_s_m3=float((b-.25*(a-c)*frac)/1e6)))
    return peaks


def verify():
    os.makedirs(OUT, exist_ok=True)
    device, mesh = Device(), trumpet_mesh()
    s = 1j*np.geomspace(1e-5, 50., 2000)
    reflection = -(1+.167*s)/(1+1.393*s+.457*s*s)
    assert max(abs(reflection)) <= 1.
    for order in (0,2):
        rates, weights, error = fit(order)
        assert np.all(rates>0) and np.all(weights>=0) and error<.03
    ref_f, ref_z, ref_dc = impedance(device, mesh, {})
    new_f, new_z, new_dc = impedance(device, mesh, MODEL)
    assert np.array_equal(ref_f, new_f)
    np.savez(OUT/'input_impedance.npz', frequency_hz=new_f,
             previous_pa_s_m3=ref_z, modern_pa_s_m3=new_z)

    note = next(s for s in performance() if s['note']=='Bb4')
    on, _ = simulate(device, note, .45, probe=True, **MODEL)
    off, _ = simulate(device, note, .45, probe=True, feedback=False, **MODEL)
    coupled = analyze(note, on, .25, .45)
    decoupled = analyze(note, off, .25, .45)
    report = dict(reflection_max_magnitude=float(max(abs(reflection))),
                  end_pressure_pa=dict(previous=ref_dc, modern=new_dc),
                  peaks=dict(previous=resonances(ref_f,ref_z), modern=resonances(new_f,new_z)),
                  coupled=coupled, without_return_coupling=decoupled)
    (OUT/'verification.json').write_text(json.dumps(report, indent=2))
    assert coupled['boundary_residual_max_pa']<.02 and abs(new_dc)<.1
    assert coupled['mouthpiece_rms_pa']>2*decoupled['mouthpiece_rms_pa']
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    verify()
