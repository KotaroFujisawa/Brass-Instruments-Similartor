"""Axial radius profiles for the trumpet mouthpiece, in SI units."""
from trumpet import ROOT, json, mouthpiece, np

MEASURED_PROFILE = ROOT/'data/bach_5c_jackson_2025_fig2.csv'
ADAPTER_LENGTH_M = .0018


def profile(name='legacy',depth=None):
    """Return stations and radii for the isolated mouthpiece."""
    if name=='legacy':
        if depth is None:
            depth = json.loads((ROOT/'data/model.json').read_text())['cup_depth_m']
        return mouthpiece(depth)
    if name!='detailed':
        raise ValueError(f'Unknown trumpet mouthpiece profile {name!r}')

    z,radius = np.loadtxt(MEASURED_PROFILE,delimiter=',',skiprows=1,unpack=True)
    if not (np.all(np.isfinite(z)) and np.all(np.isfinite(radius)) and
            np.all(np.diff(z)>0) and np.all(radius>0)):
        raise ValueError('The mouthpiece profile requires increasing stations and positive radii.')
    return z-z[0],radius


def throat_limits(z,radius):
    """Mesh anchors around the near-constant-radius part of the measured throat."""
    minimum = np.argmin(radius)
    within_tolerance = radius<=radius[minimum]+.00004
    close = np.flatnonzero(within_tolerance)
    first, last = close[0], close[-1]
    if first==last:
        first, last = first-1,last+1
    if not 0<first<last<len(z)-1:
        raise ValueError('The profile must resolve a throat inside the mouthpiece.')
    return z[first],z[last]
