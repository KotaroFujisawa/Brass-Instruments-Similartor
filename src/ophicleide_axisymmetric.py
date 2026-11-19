"""Reconstructed nine-key B-flat Gautrot ophicleide, in SI units."""
from functools import cache

from axisymmetric_mesh import annular_mesh
from trumpet import ROOT, np
from scipy.interpolate import PchipInterpolator


BORE = ROOT/'data/ophicleide_gautrot_eu3590.csv'
PRESET = ROOT/'data/performance_ophicleide.json'
BORE_SOURCE = (
    'Gautrot EU3590: Campbell, Gilbert & Myers Fig. 7.34 and Campbell, Chick & Myers Fig. 2; '
    'terminal flare inferred from the bell of a comparable B-flat Gautrot instrument'
)
MOUTHPIECE_SOURCE = (
    'Reconstructed from Caussinus historical bass-ophicleide dimensions, '
    'as transcribed by Segress (1965), p. 24'
)
HOLE_SOURCE = (
    'Reconstructed keyed chimneys; hole order from the B-flat ophicleide fingering chart '
    'and spacing from Campbell, Gilbert & Myers §4.4'
)

CUP_DEPTH_M = .034
THROAT_END_M = .038
MOUTHPIECE_LENGTH_M = .072
ADAPTER_LENGTH_M = .003
BORE_START_M = MOUTHPIECE_LENGTH_M+ADAPTER_LENGTH_M
AXIAL_CALIBRATION = .985

# The first key shuts its normally open hole. Other keys open their holes.
# Entries are centerline distance from the mouthpipe entrance, radius and
# physical chimney height. The terminal flare and port dimensions are estimates.
HOLES = (
    (2.561,.0335,.0070),
    (2.421,.0235,.0065),
    (2.281,.0215,.0060),
    (2.141,.0205,.0055),
    (1.909,.0200,.0055),
    (1.746,.0180,.0050),
    (1.432,.0168,.0050),
    (1.212,.0155,.0045),
    (1.136,.0145,.0045),
)
DEFAULT_OPEN_HOLES = (1,)
LIPS = dict(width=.011,thickness=.003,joint=.0060,rest_x=.0014,
            mass_product=4.6,quality=4.5,cup_area=np.pi*.0135**2)


@cache
def _bore():
    return np.loadtxt(BORE,delimiter=',',skiprows=1).T


def open_toneholes(open_holes=DEFAULT_OPEN_HOLES, tuning=0.):
    """Return the open `(z, radius, physical chimney height)` ports."""
    if len(set(open_holes))!=len(open_holes) or any(
            not isinstance(number,int) or number<1 or number>len(HOLES) for number in open_holes):
        raise ValueError('Open holes must be distinct numbers from 1 through 9.')
    scale = AXIAL_CALIBRATION+tuning/_bore()[0][-1]
    if not np.isfinite(scale) or scale<=0:
        raise ValueError('Tuning must leave a positive bore length.')
    return tuple((BORE_START_M+scale*x,r,h) for x,r,h in (HOLES[i-1] for i in sorted(open_holes)))


def profile(tuning=0.):
    """Return reconstructed cup, throat, backbore and bore from the lip plane."""
    x_bore,r_bore = _bore()
    scale = AXIAL_CALIBRATION+tuning/x_bore[-1]
    if not np.isfinite(scale) or scale<=0:
        raise ValueError('Tuning must leave a positive bore length.')
    cup_z = np.linspace(0.,CUP_DEPTH_M,25)
    # Smooth funnel into the historical 7 mm choke, then a straight backbore.
    cup_r = .0035+.010*(1-(cup_z/CUP_DEPTH_M)**1.55)**1.25
    return (np.r_[cup_z,THROAT_END_M,MOUTHPIECE_LENGTH_M,BORE_START_M+scale*x_bore],
            np.r_[cup_r,.0035,.005,r_bore])


def refined_ophicleide_mesh(open_holes=DEFAULT_OPEN_HOLES, tuning=0., *, axial_scale=1., dr=.006):
    """Mesh the complete bore and distribute discrete side ports on its wall."""
    return _mesh(tuple(sorted(open_holes)),tuning,axial_scale,dr)


@cache
def _mesh(open_holes,tuning,axial_scale,dr):
    if not all(np.isfinite(value) and value>0 for value in (axial_scale,dr)):
        raise ValueError('Mesh spacing must be positive and finite.')
    holes = open_toneholes(open_holes,tuning)
    x,r = profile(tuning)
    bell_end = float(x[-1])
    anchors = [0.,CUP_DEPTH_M,THROAT_END_M,MOUTHPIECE_LENGTH_M,BORE_START_M,
               bell_end-.55,bell_end-.16,bell_end-.04,bell_end]
    if any(b<=a for a,b in zip(anchors[:-1],anchors[1:])):
        raise ValueError('Tuning leaves insufficient length for the bell and mouthpipe.')
    spacing = np.array([.0026,.0020,.0028,.0024,.010,.008,.006,.0045])*axial_scale
    z = np.concatenate([np.linspace(a,b,1+int(np.ceil((b-a)/dz)))[:-1]
                        for a,b,dz in zip(anchors[:-1],anchors[1:],spacing)] + [[bell_end]])
    r_z = PchipInterpolator(x,r)(z)
    mesh = annular_mesh(z,r_z,dr,toneholes=holes)
    mesh.update(bore_profile=BORE.name,mouthpiece_profile=MOUTHPIECE_SOURCE,
                mouthpiece_length_m=MOUTHPIECE_LENGTH_M,mouthpiece_throat_radius_m=.0035,
                open_hole_numbers=open_holes,tuning_m=tuning,bore_axial_calibration=AXIAL_CALIBRATION)
    return mesh
