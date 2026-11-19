"""Digitized B-flat euphonium bore and reconstructed Boosey mouthpiece."""
from functools import cache

from axisymmetric_mesh import annular_mesh
from trumpet import C, ROOT, np
from scipy.interpolate import PchipInterpolator


BORE = ROOT/'data/euphonium_bore_fig221.csv'
MOUTHPIECE = ROOT/'data/boosey_euphonium_mouthpiece_reconstructed.csv'
BORE_SOURCE = 'Campbell, Gilbert & Myers, The Science of Brass Instruments, Fig. 2.21'
MOUTHPIECE_SOURCE = 'Boosey & Co, Royal Birmingham Conservatoire catalogue 13.15'
VALVE_SEMITONES = {1: 2, 2: 1, 3: 3, 4: 5}
VALVE_WINDOW_M = (.220,.414)
BB1_HZ = 440*2**((34-69)/12)
NOMINAL_LENGTH_M = C/(2*BB1_HZ)
RECEIVER_LENGTH_M = .004
LIPS = dict(width=.0112, thickness=.0030, joint=.0061, rest_x=.00155,
            mass_product=3.5, quality=3.8, cup_area=np.pi*(.0269/2)**2)


@cache
def _profiles():
    bore = np.loadtxt(BORE,delimiter=',',skiprows=1).T
    mouthpiece = np.loadtxt(MOUTHPIECE,delimiter=',',skiprows=1).T
    return bore,mouthpiece


def valve_extension(valves):
    """Approximate the fourth valve's compensating loops by the ideal net length."""
    semitones = [VALVE_SEMITONES[valve] for valve in valves]
    if 4 in valves:
        return NOMINAL_LENGTH_M*(2**(sum(semitones)/12)-1)
    return NOMINAL_LENGTH_M*sum(2**(offset/12)-1 for offset in semitones)


def profile(valves=(), tuning=0.):
    valves = tuple(valves)
    if len(set(valves))!=len(valves) or any(valve not in VALVE_SEMITONES for valve in valves):
        raise ValueError('Use distinct valve numbers 1–4.')
    (xb,rb),(xm,rm) = _profiles()
    added = valve_extension(valves)+tuning
    first,last = VALVE_WINDOW_M
    if not np.isfinite(added) or added<=first-last:
        raise ValueError('Tuning must leave a positive length in the valve section.')
    adjusted_xb = xb+added*np.clip((xb-first)/(last-first),0.,1.)
    return np.r_[xm,xm[-1]+RECEIVER_LENGTH_M+adjusted_xb],np.r_[rm,rb]


def refined_euphonium_mesh(valves=(), tuning=0., *, axial_scale=1., dr=.006):
    """Produce annular finite volumes without changing the mouthpiece or bell."""
    return _mesh(tuple(sorted(valves)),tuning,axial_scale,dr)


@cache
def _mesh(valves,tuning,axial_scale,dr):
    if not all(np.isfinite(value) and value>0. for value in (axial_scale,dr)):
        raise ValueError('Mesh spacing must be positive and finite.')
    x,r = profile(valves,tuning)
    mouthpiece_end = float(_profiles()[1][0][-1])
    throat = float(x[np.argmin(r)])
    receiver_end = mouthpiece_end+RECEIVER_LENGTH_M
    end = float(x[-1])
    first,last = VALVE_WINDOW_M
    added = valve_extension(valves)+tuning
    anchors = np.array([
        0.,.012,.028,throat,throat+.004,mouthpiece_end,receiver_end,
        receiver_end+first,receiver_end+last+added,end-.36,end-.077,end-.030,end,
    ])
    spacing = axial_scale*np.array([
        .0034,.0031,.0028,.0034,.0029,.0034,.008,.008,.008,.006,.005,.004,
    ])
    segments = np.ceil(np.diff(anchors)/spacing).astype(int)
    edges = np.concatenate([np.linspace(a,b,n+1)[:-1]
                            for a,b,n in zip(anchors[:-1],anchors[1:],segments)]
                           + [[end]])
    radius = PchipInterpolator(x,r)(edges)
    radius[-1] = r[-1]
    mesh = annular_mesh(edges,radius,dr)
    mesh.update(
        bore_profile=BORE.name,mouthpiece_profile=MOUTHPIECE.name,
        mouthpiece_length_m=mouthpiece_end,mouthpiece_throat_radius_m=float(r.min()),
        receiver_length_m=RECEIVER_LENGTH_M,
        valve_length_m=valve_extension(valves),tuning_m=tuning,
        bore_source=BORE_SOURCE,mouthpiece_source=MOUTHPIECE_SOURCE,
    )
    return mesh
