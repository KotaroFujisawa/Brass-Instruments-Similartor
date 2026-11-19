"""Full-double F/B-flat horn with a common leadpipe, mouthpiece and bell."""
from functools import cache

from scipy.interpolate import PchipInterpolator

from axisymmetric_mesh import annular_mesh
from trumpet import C, ROOT, np


BORE = ROOT/'data/horn_bb_bore_fig734.csv'
MOUTHPIECE = ROOT/'data/alexander_8m_mouthpiece_reconstructed.csv'
BORE_SOURCE = 'Campbell, Gilbert & Myers, The Science of Brass Instruments, Fig. 7.34'
MOUTHPIECE_SOURCE = 'Alexander 8M dimensions; profile fit to Campbell et al., Fig. 7.55c'
BELL_RADIUS_M = .155
BELL_FLARE_EXPONENT = 1.5
RECEIVER_LENGTH_M = .004
VALVE_WINDOW_M = (.7359, 1.2836)
VALVE_SEMITONES = {1: 2, 2: 1, 3: 3}
NOMINAL_LENGTH_M = {side: C/(2*440*2**((midi-69)/12))
                    for side,midi in (('Bb',34),('F',29))}
F_LOOP_M = NOMINAL_LENGTH_M['F']-NOMINAL_LENGTH_M['Bb']
LIPS = dict(width=.0078, thickness=.0021, joint=.0042, rest_x=.00125,
            mass_product=2.7, quality=3.8, cup_area=np.pi*(.0175/2)**2)


@cache
def _profiles():
    x,r = np.loadtxt(BORE,delimiter=',',skiprows=1).T
    mouthpiece = np.loadtxt(MOUTHPIECE,delimiter=',',skiprows=1).T

    # Fig. 7.34 cuts off the flare. Continue its terminal slope with a
    # generalized horn flare, anchored at the modern Model 103 bell radius.
    exponent = BELL_FLARE_EXPONENT
    slope = PchipInterpolator(x,r).derivative()(x[-1])
    flare_scale = exponent*r[-1]/slope
    bell_radius = np.geomspace(r[-1],BELL_RADIUS_M,25)[1:]
    bell_x = x[-1] + flare_scale*(1-(r[-1]/bell_radius)**(1/exponent))
    return (np.r_[x,bell_x],np.r_[r,bell_radius]),mouthpiece


def valve_extension(side, valves):
    """Sum the independent valve loops of the selected side."""
    if side not in NOMINAL_LENGTH_M:
        raise ValueError('Horn side must be F or Bb.')
    valves = tuple(valves)
    if len(set(valves))!=len(valves) or any(valve not in VALVE_SEMITONES for valve in valves):
        raise ValueError('Use distinct valve numbers 1–3.')
    return NOMINAL_LENGTH_M[side]*sum(2**(VALVE_SEMITONES[valve]/12)-1 for valve in valves)


def profile(side, valves=(), tuning=0.):
    """Return axial stations and radii in metres for one change-valve setting."""
    (xb,rb),(xm,rm) = _profiles()
    extension = valve_extension(side,valves)+tuning+(F_LOOP_M if side=='F' else 0.)
    first,last = VALVE_WINDOW_M
    if not np.isfinite(extension) or extension<=first-last:
        raise ValueError('Tuning must leave a positive length in the valve section.')
    adjusted_xb = xb+extension*np.clip((xb-first)/(last-first),0.,1.)
    return np.r_[xm,xm[-1]+RECEIVER_LENGTH_M+adjusted_xb],np.r_[rm,rb]


def refined_horn_mesh(side, valves=(), tuning=0., *, axial_scale=1., dr=.0055):
    """Produce annular cells, with finer spacing in the throat and bell."""
    return _mesh(side,tuple(sorted(valves)),tuning,axial_scale,dr)


@cache
def _mesh(side,valves,tuning,axial_scale,dr):
    if not all(np.isfinite(value) and value>0. for value in (axial_scale,dr)):
        raise ValueError('Mesh spacing must be positive and finite.')
    x,r = profile(side,valves,tuning)
    mouthpiece_end = float(_profiles()[1][0][-1])
    receiver_end = mouthpiece_end+RECEIVER_LENGTH_M
    end = float(x[-1])
    anchors = np.array([0.,.009,.028,.033,mouthpiece_end,receiver_end,
                        receiver_end+.52,end-.48,end-.135,end-.045,end])
    spacing = axial_scale*np.array([.003,.0034,.003,.004,.004,.008,.010,.006,.0045,.0035])
    segments = np.ceil(np.diff(anchors)/spacing).astype(int)
    edges = np.concatenate([np.linspace(a,b,n+1)[:-1]
                            for a,b,n in zip(anchors[:-1],anchors[1:],segments)]
                           + [[end]])
    radius = PchipInterpolator(x,r)(edges)
    radius[-1] = BELL_RADIUS_M
    mesh = annular_mesh(edges,radius,dr)
    mesh.update(
        horn_side=side,bore_profile=BORE.name,mouthpiece_profile=MOUTHPIECE.name,
        mouthpiece_length_m=mouthpiece_end,mouthpiece_throat_radius_m=.002,
        receiver_length_m=RECEIVER_LENGTH_M,valve_length_m=valve_extension(side,valves),
        side_loop_length_m=F_LOOP_M if side=='F' else 0.,tuning_m=tuning,
        bore_source=BORE_SOURCE,mouthpiece_source=MOUTHPIECE_SOURCE,
    )
    return mesh
