"""Measure the sealed-shank 5C resonance in the axisymmetric pulse solver."""
import json
import os
from types import SimpleNamespace

from axisymmetric_mesh import annular_mesh
from render_trumpet_mouthpiece import Device, OUT, Run, TRACE_FS, mesh_for
from scipy.fft import rfft, rfftfreq
from trumpet import C, RHO, np


def closed_shank(mesh):
    """Reflection b=a in every outgoing annulus, that is zero axial velocity."""
    area = mesh['faces'][mesh['bell']:]['area'].astype(float)
    rings = len(area)
    return SimpleNamespace(
        basis=np.diag(np.sqrt(area.sum()/area)), A=np.array([[-1.]]),
        B=np.zeros((1,rings)), C=np.zeros((rings,1)), D=np.eye(rings),
    )


def verify_sealed_mouthpiece():
    mounted = mesh_for('detailed',dict(valves_semitones=[]))
    z = mounted['z']
    end = int(np.argmin(abs(z-mounted['mouthpiece_length_m'])))
    if not np.isclose(z[end],mounted['mouthpiece_length_m'],atol=1e-10):
        raise ValueError('The mesh must have a face at the measured shank end.')
    mesh = annular_mesh(z[:end+1],mounted['radius'][:end+1])
    drive_amplitude = 12.
    with Run(Device(),mesh,.36,mode=3,pressure=drive_amplitude,thermoviscous=True,
             reservoir_inflow=True,exterior_load=closed_shank(mesh),record_bell=True) as run:
        run.until()
        trace = run.history()
        velocity = run.bell_history()
        dt = run.dt
    time = (np.arange(len(trace))+1)/TRACE_FS-dt
    drive_pa = drive_amplitude*np.exp(-.5*((time-.00045)/.00007)**2)
    area = np.pi*mesh['radius'][0]**2
    flow = area/(RHO*C)*drive_pa
    nfft = 2**19
    frequency = rfftfreq(nfft,1/TRACE_FS)
    band = (frequency>700)&(frequency<6800)
    f = frequency[band]
    # The charge injected by a one-sided pulse leaves a DC pressure in a
    # sealed cavity. Evaluate both signals with equal exponential damping so
    # the finite record does not turn this DC offset into spurious peaks.
    time_weight = np.exp(-time/.03)
    z = rfft((trace[:,8]+drive_pa)*time_weight,nfft)[band]/rfft(flow*time_weight,nfft)[band]
    mag = abs(z)
    found = []
    for lower,upper in ((850,1550),(2950,3900),(5100,6100)):
        indices = np.flatnonzero((f>=lower)&(f<=upper))
        i = indices[np.argmax(mag[indices])]
        if i in (indices[0],indices[-1]):
            raise RuntimeError('No isolated resonance in the expected frequency band.')
        a,b,c = mag[i-1:i+2]
        delta = .5*(a-c)/(a-2*b+c)
        found.append(float(f[i]+delta*(f[1]-f[0])))
    result = dict(simulated_resonance_hz=found,published_measurement_hz=[1216.9,3426.5,5607.1],
                  shank='sealed',adapter_included=False,cells=len(mesh['cells']),time_step_s=dt,
                  exponential_window_time_constant_s=.03,
                  maximum_shank_velocity_m_s=float(np.max(abs(velocity))))
    if len(found)!=3 or result['maximum_shank_velocity_m_s']>1e-4:
        raise RuntimeError(f'Sealed mouthpiece verification failed: {result}')
    os.makedirs(OUT,exist_ok=True)
    (OUT/'mouthpiece_sealed_euler.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps(verify_sealed_mouthpiece(),indent=2))
