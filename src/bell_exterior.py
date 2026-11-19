"""Passive modal load and exterior sound for an equivalent spherical bell cap.

Plane annuli drive the front hemisphere of a sphere of the bell radius. The
opposite hemisphere is rigid. The mapping preserves volume flux and acoustic
power. Each retained spherical harmonic satisfies the outgoing Helmholtz
equation exactly. The spherical baffle is an approximation to the bell rim
and exterior wall geometry.
"""
from dataclasses import dataclass
from math import factorial

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.fft import irfft, next_fast_len, rfft, rfftfreq
from scipy.linalg import block_diag
from scipy.signal import residue
from scipy.special import eval_legendre, j1, jn_zeros


def spherical_polynomials(order):
    """Numerator of the modified spherical Hankel function and its derivative.

    In the dimensionless Laplace variable x, k_l(x)=exp(-x)*N(x)/x**(l+1).
    The derivative satisfies -k_l'(x)=exp(-x)*D(x)/x**(l+2).
    An overall common normalization of k_l does not affect either transfer.
    """
    numerator = np.array([factorial(order+i)/(2**i*factorial(i)*factorial(order-i))
                          for i in range(order+1)])
    denominator = np.polysub(np.polymul([1, order+1], numerator),
                             np.polymul([1, 0], np.polyder(numerator)))
    return numerator, denominator


def impedance_states(numerator, denominator):
    """Real pole-residue realization of x*N(x)/D(x) with feedthrough one."""
    residues, poles, direct = residue(np.r_[numerator, 0.], denominator)
    if not np.allclose(direct, [1.]) or np.any(poles.real>=0):
        raise RuntimeError('Unexpected poles or feedthrough in spherical impedance')
    blocks, inputs, outputs = [], [], []
    for value, pole in zip(residues, poles):
        if pole.imag>0:
            blocks.append([[pole.real, -pole.imag], [pole.imag, pole.real]])
            inputs.extend([1., 0.])
            outputs.extend([2*value.real, -2*value.imag])
        elif pole.imag==0:
            blocks.append([[pole.real]])
            inputs.append(1.)
            outputs.append(value.real)
    return block_diag(*blocks), np.array(inputs), np.array(outputs)


def aperture_modes(edges, ports):
    """Area-average rigid circular-duct modes in each of the annuli."""
    weights = np.diff(edges**2)
    shapes = np.ones((len(weights), ports))
    if ports>1:
        roots = jn_zeros(1, ports-1)
        integral = edges[:, None]*j1(edges[:, None]*roots[None, :])/roots[None, :]
        shapes[:, 1:] = 2*np.diff(integral, axis=0)/weights[:, None]
    scale = np.sqrt(weights)
    q, r = np.linalg.qr(scale[:, None]*shapes)
    return q*np.sign(np.diag(r))[None, :]/scale[:, None], weights


def sphere_projection(edges, basis, order):
    """Return modal fluxes on the sphere per unit modal aperture velocity.

    v_radial = cos(theta)*u_axial, so v_radial*dS = u_axial*dS_projected.
    The sphere modes sqrt(2l+1)*P_l(cos(theta)) are orthonormal in area
    fraction over the whole sphere.
    """
    nodes, quadrature = leggauss(order+2)
    mu = np.sqrt(np.maximum(0., 1-edges**2))
    half_widths = (mu[:-1]-mu[1:])/2
    midpoints = (mu[:-1]+mu[1:])/2
    positions = midpoints[:, None]+half_widths[:, None]*nodes[None, :]
    integral_weights = .5*half_widths[:, None]*quadrature[None, :]*positions
    integrals = np.array([
        np.sum(eval_legendre(l, positions)*np.sqrt(2*l+1)*integral_weights, axis=1)
        for l in range(order+1)
    ])
    return integrals@basis


@dataclass(frozen=True)
class BellExterior:
    """Continuous scattering matrices use dimensionless time c*t/radius.

    Incident and reflected modal pressures are a=(p+rho*c*u)/2 and
    b=(p-rho*c*u)/2. Port velocity is positive out of the bore. The matrices
    act in either SI pressure units or a common nondimensional pressure unit.
    """
    radius_m: float
    basis: np.ndarray
    weights: np.ndarray
    projection: np.ndarray
    polynomials: tuple
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    sound_speed_m_s: float
    density_kg_m3: float

    def scattering(self, frequency_hz):
        """Frequency-domain reflection including coupling between ports."""
        x = 2j*np.pi*np.atleast_1d(frequency_hz)*self.radius_m/self.sound_speed_m_s
        diagonal = np.column_stack([x*np.polyval(n, x)/np.polyval(d, x)
                                    for n, d in self.polynomials])
        normalized = 2*self.projection
        impedance = np.einsum('li,fl,lj->fij', normalized, diagonal, normalized)
        identity = np.eye(self.basis.shape[1])
        return np.linalg.solve(impedance+identity, impedance-identity)

    def receiver_transfer(self, frequency_hz, angles_deg=(0., 45.), distance_m=1.,
                          *, remove_travel_time=True):
        """Pressure transfer from aperture modal velocity, in Pa per m/s.

        Distance is measured from the equatorial center of the sphere. Removing
        travel time references arrival to propagation from the virtual surface
        at radius a, with the common (distance-a)/c delay omitted.
        """
        if distance_m<=self.radius_m:
            raise ValueError('Distance must be outside the matching sphere')
        x = 2j*np.pi*np.atleast_1d(frequency_hz)*self.radius_m/self.sound_speed_m_s
        ratio = self.radius_m/distance_m
        y = x/ratio
        mu = np.cos(np.deg2rad(np.atleast_1d(angles_deg)))
        radial = np.column_stack([
            x*ratio**(l+1)*np.polyval(n, y)/np.polyval(d, x)
            for l, (n, d) in enumerate(self.polynomials)
        ])
        if not remove_travel_time:
            radial *= np.exp(x-y)[:, None]
        angular = np.column_stack([np.sqrt(2*l+1)*eval_legendre(l, mu)
                                   for l in range(len(self.polynomials))])
        scale = self.density_kg_m3*self.sound_speed_m_s
        return scale*np.einsum('fl,al,lp->fap', radial, angular, self.projection)

    def radiate(self, modal_velocity_m_s, fs, angles_deg=(0., 45.), distance_m=1.,
                *, remove_travel_time=True):
        """Return receiver pressures at the input sample rate, columns by angle."""
        velocity = np.asarray(modal_velocity_m_s)
        if velocity.ndim!=2 or velocity.shape[1]!=self.basis.shape[1]:
            raise ValueError('Velocity must have shape (samples, aperture modes)')
        if len(velocity)==0 or fs<=0:
            raise ValueError('Supply a nonempty history and positive sample rate')
        delay = 0. if remove_travel_time else distance_m/self.sound_speed_m_s
        size = next_fast_len(len(velocity)+int(np.ceil((.05+delay)*fs)), real=True)
        transfer = self.receiver_transfer(rfftfreq(size, 1/fs), angles_deg, distance_m,
                                          remove_travel_time=remove_travel_time)
        spectrum = np.einsum('fap,fp->fa', transfer, rfft(velocity, n=size, axis=0))
        return irfft(spectrum, n=size, axis=0)[:len(velocity)]


def exterior_model(ring_edges_m, ports=4, exterior_order=10, *,
                   sound_speed_m_s=343., density_kg_m3=1.2):
    """Construct a passive, reciprocal radiation load from physical ring edges."""
    edges = np.asarray(ring_edges_m, dtype=float)
    if (edges.ndim!=1 or len(edges)<2 or edges[0]!=0
            or not np.isfinite(edges).all() or np.any(np.diff(edges)<=0)):
        raise ValueError('Ring edges must increase strictly from the symmetry axis')
    if not 1<=ports<=len(edges)-1 or not 0<=exterior_order<=12:
        raise ValueError('Invalid aperture basis size or spherical harmonic order')
    if sound_speed_m_s<=0 or density_kg_m3<=0:
        raise ValueError('Sound speed and density must be positive')
    radius = edges[-1]
    basis, weights = aperture_modes(edges/radius, ports)
    projection = sphere_projection(edges/radius, basis, exterior_order)
    # The aperture is one quarter of the sphere's area. Divide by its area
    # fraction in the power identity, or multiply each coupling by 2.
    normalized = 2*projection
    polynomials = tuple(spherical_polynomials(l) for l in range(exterior_order+1))
    blocks, input_blocks, output_blocks = [], [], []
    for coefficients, coupling in zip(polynomials, normalized):
        state, drive, response = impedance_states(*coefficients)
        blocks.append(state)
        input_blocks.append(np.outer(drive, coupling))
        output_blocks.append(np.outer(coupling, response))
    impedance_a = block_diag(*blocks)
    impedance_b = np.concatenate(input_blocks, axis=0)
    impedance_c = np.concatenate(output_blocks, axis=1)
    impedance_d = normalized.T@normalized
    identity = np.eye(ports)
    c = np.linalg.solve(identity+impedance_d, impedance_c)
    d = np.linalg.solve(identity+impedance_d, impedance_d-identity)
    a = impedance_a-impedance_b@c
    b = impedance_b@(identity-d)
    return BellExterior(radius, basis, weights, projection, polynomials, a, b, c, d,
                        sound_speed_m_s, density_kg_m3)
