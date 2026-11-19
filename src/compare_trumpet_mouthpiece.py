"""Compare linear input impedance across trumpet mouthpiece shapes and meshes."""
from pathlib import Path
import argparse
import json
import os

from trumpet import C, RHO, ROOT, SCUP, np, signal
from scipy.special import ive
from axisymmetric_mesh import refined_trumpet_mesh
from bell_exterior import exterior_model
from thermoviscous import GAMMA_MINUS_ONE, KINEMATIC_VISCOSITY, THERMAL_DIFFUSIVITY

OUT = ROOT/'output/trumpet_mouthpiece'
REFERENCE = RHO*C/SCUP
CASES = ('legacy', 'legacy_fine', 'detailed', 'detailed_finer')
LABELS = ('Original profile, original cells', 'Original profile, finer cells',
          'Digitized 5C, finer cells', 'Digitized 5C, half-size cells')
MEASURED_CLOSED_HZ = [1216.9, 3426.5, 5607.1]
PUBLISHED_BESSEL_CLOSED_HZ = [1132.3, 3346.1, 5597.9]


def tube_impedance(frequency, edges, radii, loads):
    """Propagate pressure and volume flow back through cylindrical midpoint slices.

    Uses the e^(i omega t) convention and the full Bessel functions for local
    circular-duct viscous density and isothermal-wall thermal compliance.
    `None` denotes a sealed downstream end.
    """
    frequency = np.asarray(frequency, dtype=float)
    lengths = np.diff(edges)
    radius = .5*(radii[:-1]+radii[1:])
    if (np.any(frequency<=0) or np.any(lengths<=0) or np.any(radius<=0)
            or len(lengths)!=len(radius)):
        raise ValueError('Frequencies, cylinder lengths and radii must all be positive.')
    s = 2j*np.pi*frequency
    pressure = np.column_stack([np.ones_like(s) if load is None else load for load in loads])
    flow = np.column_stack([np.zeros_like(s) if load is None else np.ones_like(s) for load in loads])

    # Reuse viscothermal terms for repeated radii, including the throat.
    unique, index = np.unique(np.round(radius,12), return_inverse=True)
    argument = unique[:,None]*np.sqrt(s[None,:]/KINEMATIC_VISCOSITY)
    viscous = 1-2*ive(1,argument)/(argument*ive(0,argument))
    argument = unique[:,None]*np.sqrt(s[None,:]/THERMAL_DIFFUSIVITY)
    thermal = 1+2*GAMMA_MINUS_ONE*ive(1,argument)/(argument*ive(0,argument))
    propagation = (s/C)*np.sqrt(thermal/viscous)
    z_characteristic = (RHO*C/(np.pi*unique[:,None]**2))/np.sqrt(viscous*thermal)

    # Divide each two-port matrix by cosh(gamma*length); impedance is unchanged.
    for length, idx in zip(lengths[::-1],index[::-1]):
        zc = z_characteristic[idx,:,None]
        coupling = np.tanh(propagation[idx,:,None]*length)
        pressure, flow = pressure+zc*coupling*flow, flow+(coupling/zc)*pressure
    return pressure/flow


def shank_load(frequency, radius):
    """Uniform tube with a causal unflanged opening (Silva et al.)."""
    s = 2j*np.pi*frequency*radius/C
    reflection = -(1+.167*s)/(1+1.393*s+.457*s*s)
    return (RHO*C/(np.pi*radius**2))*((1+reflection)/(1-reflection))


def bell_load(frequency, mesh):
    """Mean pressure per volume velocity for uniform velocity on the bell."""
    radius = mesh['radius'][-1]
    bell = exterior_model(np.linspace(0,radius,int(mesh['rings'][-1])+1))
    reflection = bell.scattering(frequency)
    identity = np.eye(bell.basis.shape[1])
    impedance = np.linalg.solve(identity-reflection,identity+reflection)
    return (RHO*C/(np.pi*radius**2))*impedance[:,0,0]


def locate_edge(mesh, key):
    z = mesh['z']
    index = int(np.argmin(abs(z-mesh[key])))
    if abs(z[index]-mesh[key])>1e-10:
        raise ValueError(f'Mesh does not contain a boundary at {key}.')
    return index


def profile_metrics(mesh):
    close = locate_edge(mesh,'mouthpiece_length_m')
    x,z = mesh['z'][:close+1],mesh['radius'][:close+1]
    volume = np.sum(np.pi*np.diff(x)*(z[:-1]**2+z[:-1]*z[1:]+z[1:]**2)/3)
    return dict(length_m=float(x[-1]), mouthpiece_axial_cells=close,
                internal_volume_ml=float(volume*1e6),
                inlet_radius_mm=float(z[0]*1e3), minimum_radius_mm=float(min(z)*1e3),
                exit_radius_mm=float(z[-1]*1e3), leadpipe_start_m=mesh['leadpipe_start_m'])


def peaks(frequency, impedance):
    magnitude = 20*np.log10(np.maximum(abs(impedance)/REFERENCE,np.finfo(float).tiny))
    indices = signal.find_peaks(magnitude,prominence=1.)[0]
    result = []
    for i in indices:
        y0,y1,y2 = magnitude[i-1:i+2]
        offset = .5*(y0-y2)/(y0-2*y1+y2)
        result.append(dict(frequency_hz=float(frequency[i]+offset*(frequency[1]-frequency[0])),
                           level_db=float(y1-.25*(y0-y2)*offset)))
    return result


def match_peaks(reference, comparison):
    """Match neighbouring resonances in increasing frequency."""
    limit = (min(np.diff([p['frequency_hz'] for p in reference]))/2
             if len(reference)>1 else 150.)
    unused = comparison.copy()
    result = []
    for anchor in reference:
        if not unused:
            break
        f = anchor['frequency_hz']
        closest = min(unused,key=lambda p:abs(p['frequency_hz']-f))
        freq = closest['frequency_hz']
        if abs(freq-f)>limit:
            continue
        unused.remove(closest)
        result.append(dict(reference_hz=f, comparison_hz=freq, shift_hz=freq-f,
                           shift_cents=float(1200*np.log2(freq/f)),
                           level_change_db=closest['level_db']-anchor['level_db']))
    return result


def plot_comparison(f_body, f_piece, curves, measured, output):
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2,2,figsize=(11.0,6.5),layout='constrained')
    styles = [('#636b73','-'),('#c87331','--'),('#007c83','-'),('#364f79',':')]
    for name,label,(color,style) in zip(CASES,LABELS,styles):
        values = curves[name]
        for ax, x, y in ((axes[0,0],f_body,values['full_bore']),
                         (axes[1,0],f_piece,values['open_shank']),
                         (axes[1,1],f_piece,values['closed_shank'])):
            ax.plot(x,20*np.log10(abs(y)/REFERENCE),label=label,color=color,ls=style,lw=1.2)
        axes[0,1].plot(f_body,np.angle(values['full_bore'],deg=True),color=color,ls=style,lw=1.05)
    for f in measured:
        axes[1,1].axvline(f,color='#a7a9ad',lw=.8,ls='--')
    axes[0,0].legend(frameon=False,fontsize=8)
    axes[0,0].set(title='Mouthpiece and York trumpet',ylabel='Input magnitude (dB re ρc/S₀)')
    axes[0,1].set(title='Mouthpiece and York trumpet',ylabel='Input phase (degrees)')
    axes[1,0].set(title='Detached, open shank',ylabel='Input magnitude (dB re ρc/S₀)')
    axes[1,1].set(title='Detached, closed shank',ylabel='Input magnitude (dB re ρc/S₀)')
    axes[1,1].plot([],[],color='#a7a9ad',lw=.8,ls='--',label='Measured 5C peaks')
    axes[1,1].legend(frameon=False,fontsize=8,loc='upper right')
    for ax in axes.flat:
        ax.set_xlabel('Frequency (Hz)')
        ax.grid(alpha=.16)
    fig.savefig(output/'mouthpiece_impedance.png',dpi=165)
    plt.close(fig)


def compare(output=OUT):
    output = Path(output)
    os.makedirs(output,exist_ok=True)
    grids = {
        'legacy': refined_trumpet_mesh(),
        'legacy_fine': refined_trumpet_mesh(mouthpiece_scale=.60),
        'detailed': refined_trumpet_mesh(mouthpiece_profile='detailed',mouthpiece_scale=.60),
        'detailed_finer': refined_trumpet_mesh(mouthpiece_profile='detailed',mouthpiece_scale=.30),
    }
    f_body = np.arange(25.,2601.,1.)
    f_piece = np.arange(30.,6501.,2.)
    result = dict(
        method='Linear plane-wave cylindrical midpoint transfer matrices with full Bessel thermoviscous terms',
        reference_impedance_pa_s_per_m3=REFERENCE,
        full_bore_termination='Uniform velocity on the equivalent spherical exterior load',
        open_shank_termination='Causal unflanged circular-pipe radiation impedance',
        closed_shank_termination='Rigid end at the shank outlet',
        measured_5c_closed_shank_source='Jackson and Scavone (2025), Table 1, doi:10.1121/2.0002165',
        measured_5c_closed_shank_peak_hz=MEASURED_CLOSED_HZ,
        published_bessel_tmm_5c_closed_shank_peak_hz=PUBLISHED_BESSEL_CLOSED_HZ,
        full_bore_geometry_comparison=(
            'The digitized 5C is joined to the existing York receiver by a 1.8 mm smooth adapter. '
            'The attached-profile change includes this connection. The detached cases end at '
            'the measured shank outlet.'
        ),
        variants={},
    )
    curves = {}
    downstream = {}
    for name, mesh in grids.items():
        x,r = mesh['z'],mesh['radius']
        end = locate_edge(mesh,'mouthpiece_length_m')
        join = locate_edge(mesh,'leadpipe_start_m')
        # Mouthpiece axial refinement leaves the downstream geometry unchanged.
        profile = ('legacy' if name.startswith('legacy') else 'detailed')
        if profile not in downstream:
            downstream[profile] = tube_impedance(f_body,x[join:],r[join:],[bell_load(f_body,mesh)])[:,0]
        full = tube_impedance(f_body,x[:join+1],r[:join+1],[downstream[profile]])[:,0]
        detached = tube_impedance(f_piece,x[:end+1],r[:end+1],[shank_load(f_piece,r[end]),None])
        curves[name] = dict(full_bore=full,open_shank=detached[:,0],closed_shank=detached[:,1])
        if any(not np.isfinite(values).all() or np.min(values.real)<-1e-6
               for values in curves[name].values()):
            raise RuntimeError(f'Nonfinite or active input impedance for {name}.')
        result['variants'][name] = dict(
            geometry=profile_metrics(mesh),
            peaks={target:peaks(f_body if target=='full_bore' else f_piece,values)
                   for target,values in curves[name].items()},
        )

    substitutions = (
        ('mouthpiece_mesh_only','legacy','legacy_fine'),
        ('profile_change','legacy_fine','detailed'),
        ('further_refinement','detailed','detailed_finer'),
    )
    result['differences'] = {
        label: {fixture: match_peaks(result['variants'][first]['peaks'][fixture],
                                     result['variants'][second]['peaks'][fixture])
                for fixture in curves[first]}
        for label,first,second in substitutions
    }
    detail = result['variants']['detailed']['peaks']['closed_shank']
    result['closed_shank_comparison_with_measured_5c'] = [
        dict(measured_hz=f,published_bessel_tmm_hz=published,calculated_hz=peak['frequency_hz'],
             error_vs_measurement_hz=peak['frequency_hz']-f,
             error_vs_published_bessel_tmm_hz=peak['frequency_hz']-published)
        for f,published,peak in zip(MEASURED_CLOSED_HZ,PUBLISHED_BESSEL_CLOSED_HZ,detail)
    ]
    np.savez(output/'mouthpiece_impedance_curves.npz',full_bore_hz=f_body,detached_hz=f_piece,
             **{name+'_'+fixture: data for name,sweep in curves.items()
                for fixture,data in sweep.items()})
    plot_comparison(f_body,f_piece,curves,MEASURED_CLOSED_HZ,output)
    (output/'mouthpiece_impedance.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    args = parser.parse_args()
    results = compare(args.output)
    for name, row in results['variants'].items():
        print(name,[(fixture,[round(p['frequency_hz'],2) for p in values[:5]])
                    for fixture,values in row['peaks'].items()])
