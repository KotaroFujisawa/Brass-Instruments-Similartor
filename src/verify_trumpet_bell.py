"""Check the CUDA bell boundary against the continuous exterior transfer."""
import json
import os

from render_trumpet_bell import OUT, TRACE_FS, Device, Run, build_load, refined_trumpet_mesh
from trumpet import C, RHO, np
from scipy.fft import irfft, next_fast_len, rfft, rfftfreq


def verify_bell_boundary():
    grid = refined_trumpet_mesh()
    model = build_load(grid)
    with Run(Device(),grid,.16,mode=3,pressure=60.,thermoviscous=True,
             reservoir_inflow=True,exterior_load=model,record_bell=True) as run:
        run.until()
        trace = run.history()
        velocity = run.bell_history().astype(float)
    flow = velocity@grid['faces'][grid['bell']:]['area']
    modes = velocity@(model.basis*model.weights[:,None])
    size = next_fast_len(len(modes)+round(.05*TRACE_FS),real=True)
    frequency = rfftfreq(size,1/TRACE_FS)
    reflection = model.scattering(frequency)
    identity = np.eye(len(model.D))
    impedance = np.linalg.solve(identity-reflection,identity+reflection)
    pressure_modes = RHO*C*np.einsum('fij,fj->fi',impedance,rfft(modes,n=size,axis=0))
    expected_pressure = irfft(pressure_modes,n=size,axis=0)[:len(modes),0]
    error = trace[:,10]-expected_pressure
    relative_error = np.linalg.norm(error)/np.linalg.norm(expected_pressure)
    result = dict(
        gpu_cpu_bell_pressure_relative_l2=float(relative_error),
        gpu_cpu_bell_pressure_max_pa=float(max(abs(error))),
        bell_flux_max_residual_m3_s=float(max(abs(flow-trace[:,9]))),
        peak_bell_pressure_pa=float(max(abs(trace[:,10]))),
        aperture_modes=model.basis.shape[1],
        exterior_order=len(model.polynomials)-1,
        radiation_states=len(model.A),
    )
    assert relative_error<.025 and result['bell_flux_max_residual_m3_s']<2e-7
    os.makedirs(OUT,exist_ok=True)
    (OUT/'bell_boundary_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps(verify_bell_boundary(),indent=2))
