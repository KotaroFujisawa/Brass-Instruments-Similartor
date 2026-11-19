"""Passive equivalent-aperture approximation to a horn player's right hand."""
from dataclasses import dataclass

import numpy as np

from bell_exterior import BellExterior, exterior_model


# Local area ratio in the normal (partially covered) horn of Yoshikawa & Nobara,
# Acoustical Modeling of Mutes for Brass Instruments, Fig. 4(b).
# Use it as an estimated terminal area fraction, not a measurement at the rim.
HAND_OPEN_AREA_FRACTION = (50 / 84)**2


@dataclass(frozen=True)
class HornBellLoad:
    """Connect the full bell to a smaller equivalent spherical aperture.

    Pressure is continuous and volume flow is conserved across the equivalent
    area contraction. State matrices use time scaled by the full bell radius
    so the existing solver can advance them at its usual boundary time step.
    """
    aperture: BellExterior
    radius_m: float
    open_area_fraction: float
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray

    @property
    def basis(self):
        return self.aperture.basis

    @property
    def weights(self):
        return self.aperture.weights

    @property
    def sound_speed_m_s(self):
        return self.aperture.sound_speed_m_s

    @property
    def density_kg_m3(self):
        return self.aperture.density_kg_m3

    def scattering(self, frequency_hz):
        reflection = self.aperture.scattering(frequency_hz)
        identity = np.eye(self.basis.shape[1])
        mismatch = (1 - self.open_area_fraction) / (1 + self.open_area_fraction)
        return np.linalg.solve(
            identity + mismatch * reflection, reflection + mismatch * identity,
        )

    def receiver_transfer(self, frequency_hz, angles_deg=(0., 45.), distance_m=1.,
                          *, remove_travel_time=True):
        return self.aperture.receiver_transfer(
            frequency_hz, angles_deg, distance_m, remove_travel_time=remove_travel_time,
        ) / self.open_area_fraction

    def radiate(self, modal_velocity_m_s, fs, angles_deg=(0., 45.), distance_m=1.,
                *, remove_travel_time=True):
        return self.aperture.radiate(
            modal_velocity_m_s / self.open_area_fraction, fs, angles_deg, distance_m,
            remove_travel_time=remove_travel_time,
        )


def build_horn_load(grid, open_area_fraction=HAND_OPEN_AREA_FRACTION):
    """Build the ordinary hand-position load for an axisymmetric horn grid."""
    alpha = float(open_area_fraction)
    if not 0. < alpha <= 1.:
        raise ValueError('The effective horn bell opening must be positive and at most one.')
    full_radius = float(grid['radius'][-1])
    aperture = exterior_model(np.linspace(0., full_radius*np.sqrt(alpha),
                                          int(grid['rings'][-1])+1))
    identity = np.eye(aperture.basis.shape[1])
    # Eliminate the aperture's incoming and outgoing characteristics using
    # p_aperture = p_bell and u_aperture = u_bell / alpha.
    denominator = (alpha+1)*identity - (alpha-1)*aperture.D
    c = np.linalg.solve(denominator, 2*alpha*aperture.C)
    d = np.linalg.solve(denominator, (alpha+1)*aperture.D - (alpha-1)*identity)
    a = (aperture.A + (alpha-1)/(2*alpha)*aperture.B@c)/np.sqrt(alpha)
    b = aperture.B@((alpha+1)*identity + (alpha-1)*d)/(2*alpha*np.sqrt(alpha))
    return HornBellLoad(aperture, full_radius, alpha, a, b, c, d)
