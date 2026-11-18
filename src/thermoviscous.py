"""Local circular-duct thermoviscous closure for the Euler time step."""
from functools import lru_cache

from trumpet import np
from scipy.optimize import nnls
from scipy.special import ive

KINEMATIC_VISCOSITY = 1.54e-5
THERMAL_DIFFUSIVITY = 2.17e-5
GAMMA_MINUS_ONE = .4
POLES = 10


def diffusion_response(w, order):
    """Dynamic part of the dimensionless loss, with the viscosity DC term separate."""
    x = np.sqrt(w)
    return 2*x*ive(order+1, x)/ive(order, x)


@lru_cache
def fit(order):
    """Fit positive diffusive modes to the circular cross-section Bessel response."""
    rates = np.geomspace(31 if order == 2 else 6, 2e8, POLES)
    w = 1j*np.geomspace(12, 3e7, 500)
    target = diffusion_response(w, order)
    basis = w[:, None]/(w[:, None]+rates)
    scale = abs(target)
    matrix = basis/scale[:, None]
    rhs = target/scale
    weights, _ = nnls(np.r_[matrix.real, matrix.imag], np.r_[rhs.real, rhs.imag])
    error = float(np.max(abs(basis@weights-target)/scale))
    return rates, weights, error


def memory_coefficients(radius_squared, dt):
    viscous_rates, viscous_weights, viscous_error = fit(2)
    thermal_rates, thermal_weights, thermal_error = fit(0)
    viscous_step = (KINEMATIC_VISCOSITY*dt/radius_squared)[None, :]
    thermal_step = (THERMAL_DIFFUSIVITY*dt/radius_squared)[None, :]
    viscous_poles = viscous_rates[:, None]*viscous_step
    thermal_poles = thermal_rates[:, None]*thermal_step
    viscous_update = -np.expm1(-viscous_poles)
    thermal_update = -np.expm1(-thermal_poles)
    # Average the source over the step under frozen pressure and momentum.
    viscous_loss = viscous_weights[:, None]*viscous_update/viscous_rates[:, None]
    thermal_loss = GAMMA_MINUS_ONE*thermal_weights[:, None]*thermal_update/thermal_rates[:, None]
    return np.stack([viscous_update, thermal_update, viscous_loss, thermal_loss], axis=-1), {
        'viscous_relative_max': viscous_error, 'thermal_relative_max': thermal_error,
    }
