"""Check the side ports with a weak pulse, without using a playing preset."""
from axisymmetric_cuda import Device, Run, HOLE_RADIATION_RESISTANCE, tonehole_lengths
from axisymmetric_mesh import annular_mesh
from trumpet import C, json, np


def pulse(device, mesh, duration=.018):
    with Run(device,mesh,duration,mode=3,pressure=30.,radiation_order=2,
             reservoir_inflow=True,record_toneholes=True) as simulation:
        simulation.until(.5*duration)
        first = simulation.tonehole_history()
        simulation.until()
        history = simulation.tonehole_history()
        if not np.array_equal(first,history[:len(first)]):
            raise RuntimeError('Reading a partial tonehole trace changed its samples.')
        return simulation.history(), history, simulation.trace_fs


def check(device):
    x = np.linspace(0.,1.2,301)
    radii = np.full_like(x,.016)
    geometry = ((.41,.0065,.002),(.76,.007,.003))
    opened = annular_mesh(x,radii,loss=False,toneholes=geometry)
    closed = annular_mesh(x,radii,loss=False)
    patches = opened['hole_patches']
    patch_areas = opened['faces']['area'][patches['face']]*patches['open_fraction']
    integrated = np.add.reduceat(patch_areas,opened['hole_offsets'][:-1])
    areas = np.pi*np.asarray(geometry)[:,1]**2
    np.testing.assert_allclose(integrated,areas,rtol=1e-6)
    np.testing.assert_array_equal(opened['cells'],closed['cells'])

    with Run(device,opened,.003,mode=1,record_toneholes=True) as simulation:
        simulation.until()
        rest_error = float(max(np.max(abs(simulation.field())),
                               np.max(abs(simulation.tonehole_history()))))
    if rest_error>1e-8:
        raise RuntimeError('Open toneholes disturb the equilibrium state.')

    sealed_trace,no_ports,_ = pulse(device,closed)
    vented_trace,ports,fs = pulse(device,opened)
    if no_ports.shape!=(len(sealed_trace),0) or ports.shape!=(len(vented_trace),len(geometry)):
        raise RuntimeError('Port traces do not align with the recorded mouthpiece samples.')
    arrivals = [float(np.flatnonzero(abs(flow)>.1*np.max(abs(flow)))[0]/fs)
                for flow in ports.T]
    expected_arrivals = [.00045+pos/C for pos,_,_ in geometry]
    np.testing.assert_allclose(arrivals,expected_arrivals,atol=.00035,rtol=0.)
    return_wave = vented_trace[:,8]-sealed_trace[:,8]
    reflected = float(np.flatnonzero(abs(return_wave)>.1*np.max(abs(return_wave)))[0]/fs)
    if not abs(reflected-(.00045+2*geometry[0][0]/C))<.0004:
        raise RuntimeError('A tonehole pulse did not return to the mouthpiece at the expected time.')

    k = 2*np.pi*np.geomspace(10,20000,100)/C
    chimney, relaxation = tonehole_lengths(geometry[0][1],geometry[0][2],radii[0])
    impedance = 1j*k*chimney+HOLE_RADIATION_RESISTANCE*1j*k*relaxation/(1+1j*k*relaxation)
    if np.min(impedance.real)<=0.:
        raise RuntimeError('The open tonehole radiation must dissipate acoustic energy.')
    return dict(rest_max_error=rest_error, aperture_area_m2=integrated.tolist(),
                flow_arrivals_s=arrivals, mouthpiece_reflection_arrival_s=reflected,
                returned_pressure_peak_pa=float(np.max(abs(return_wave))),
                minimum_normalized_resistance=float(np.min(impedance.real)))


if __name__=='__main__':
    print(json.dumps(check(Device()),indent=2))
