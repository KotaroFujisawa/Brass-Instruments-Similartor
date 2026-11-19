"""Reconstructed BB-flat contrabass tuba and BB-67C4 mouthpiece in SI units."""
from functools import cache

from axisymmetric_mesh import annular_mesh
from trumpet import C, ROOT, np
from scipy.interpolate import PchipInterpolator

BORE = ROOT/'data/alexander_bb_tuba_reconstructed.csv'
MOUTHPIECE = ROOT/'data/yamaha_bb67c4_reconstructed.csv'
BORE_SOURCE = ('Alexander CC (Campbell, Gilbert & Myers Figs. 7.46 and 6.9), '
               'adapted to BBb with model 163 bore and bell dimensions and an axial calibration')
MOUTHPIECE_SOURCE = 'Yamaha BB-67C4, reconstructed from Vennture measured dimensions and volumes'
VALVE_SEMITONES = {1: 2, 2: 1, 3: 3, 4: 5}
# Register spacing of the direct chart reconstruction is about 27.3 Hz.
# Shorten the reconstructed instrument to fit the BB-flat mode spacing. Keep
# the measured mouthpiece and nominal valve length increments independent.
AXIAL_CALIBRATION = .936
VALVE_WINDOW_M = (.82*AXIAL_CALIBRATION,1.26*AXIAL_CALIBRATION)
BB0_HZ = 440*2**((22-69)/12)
NOMINAL_LENGTH_M = C/(2*BB0_HZ)
LIPS = dict(width=.014, thickness=.0035, joint=.008, rest_x=.0019,
            mass_product=5.5, quality=3.5, cup_area=np.pi*(.0325374/2)**2)


@cache
def _profiles():
    x_bore,r_bore = np.loadtxt(BORE,delimiter=',',skiprows=1).T
    mouthpiece = np.loadtxt(MOUTHPIECE,delimiter=',',skiprows=1).T
    return (AXIAL_CALIBRATION*x_bore,r_bore),mouthpiece


def valve_extension(valves):
    """Estimate added loop length for each operated valve of the 4-valve tuba."""
    return NOMINAL_LENGTH_M*sum(2**(VALVE_SEMITONES[i]/12)-1 for i in valves)


def profile(valves=(), tuning=0.):
    """Return axial distance and inner radius, starting at the mouthpiece rim."""
    valves = tuple(valves)
    if len(set(valves))!=len(valves) or any(i not in VALVE_SEMITONES for i in valves):
        raise ValueError('Use distinct valve numbers 1–4.')
    (xb,rb),(xm,rm) = _profiles()
    added = valve_extension(valves)+tuning
    first,last = VALVE_WINDOW_M
    if not np.isfinite(added) or added<=first-last:
        raise ValueError('Tuning must leave a positive length in the valve section.')
    # Keep the conical branches and bell fixed. All length adjustments occur
    # in the 20.5 mm cylindrical valve section.
    shifted = xb+added*np.clip((xb-first)/(last-first),0.,1.)
    return np.r_[xm,xm[-1]+shifted[1:]],np.r_[rm,rb[1:]]


def refined_tuba_mesh(valves=(), tuning=0., *, axial_scale=1., dr=.0075):
    """Produce annular finite volumes for the two-way nonlinear acoustic solver."""
    return _mesh(tuple(sorted(valves)),tuning,axial_scale,dr)


@cache
def _mesh(valves,tuning,axial_scale,dr):
    if not all(np.isfinite(value) and value>0 for value in (axial_scale,dr)):
        raise ValueError('Mesh spacing must be positive and finite.')
    x,r = profile(valves,tuning)
    xm = _profiles()[1][0]
    mp_end = float(xm[-1])
    throat = float(x[np.argmin(r)])
    end = float(x[-1])
    added = valve_extension(valves)+tuning
    first,last = VALVE_WINDOW_M
    anchors = np.array([0.,.014,.030,throat,throat+.004,mp_end,
                        mp_end+first,mp_end+last+added,end-.75,end-.23,end-.060,end])
    spacing = axial_scale*np.array([.0036,.0034,.0032,.0041,.0034,.014,.014,
                                    .014,.012,.0085,.0065])
    lengths = np.diff(anchors)
    segments = np.ceil(lengths/spacing).astype(int)
    z = np.concatenate([np.linspace(a,b,n+1)[:-1]
                        for a,b,n in zip(anchors[:-1],anchors[1:],segments)] + [[end]])
    radius = PchipInterpolator(x,r)(z)
    radius[-1] = r[-1]
    mesh = annular_mesh(z,radius,dr)
    mesh.update(bore_profile=BORE.name,mouthpiece_profile=MOUTHPIECE.name,
                mouthpiece_length_m=mp_end,mouthpiece_throat_radius_m=float(r.min()),
                valve_length_m=valve_extension(valves),tuning_m=tuning,
                bore_axial_calibration=AXIAL_CALIBRATION)
    return mesh
