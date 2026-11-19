"""Numerical checks of the spherical acoustic load and receiver transfer."""
import argparse
import json
from pathlib import Path

from axisymmetric_mesh import trumpet_mesh
from bell_exterior import exterior_model
import numpy as np
from scipy.special import eval_legendre, spherical_jn, spherical_yn


def state_response(model, frequency_hz, dtype):
    a, b, c, d = (value.astype(dtype) for value in (model.A, model.B, model.C, model.D))
    x = 2j*np.pi*frequency_hz*model.radius_m/model.sound_speed_m_s
    identity = np.eye(len(a), dtype=np.complex64 if dtype==np.float32 else np.complex128)
    return np.array([d+c@np.linalg.solve(z*identity-a, b) for z in x.astype(identity.dtype)])


def hankel_receiver(model, frequency_hz, angles_deg, distance_m=1.):
    """Independent transfer obtained directly from the outgoing spherical Hankel functions."""
    k = 2*np.pi*frequency_hz/model.sound_speed_m_s
    mu = np.cos(np.deg2rad(angles_deg))
    result = 0.
    for l, coupling in enumerate(model.projection):
        h = spherical_jn(l,k*distance_m)-1j*spherical_yn(l,k*distance_m)
        dh = (spherical_jn(l,k*model.radius_m,derivative=True)
              -1j*spherical_yn(l,k*model.radius_m,derivative=True))
        radial = -1j*h/dh*np.exp(1j*k*(distance_m-model.radius_m))
        angular = np.sqrt(2*l+1)*eval_legendre(l,mu)
        result = result+np.einsum('f,a,p->fap',radial,angular,coupling)
    return result*model.sound_speed_m_s*model.density_kg_m3


def verification(model):
    frequency = np.geomspace(.01, 48000., 240)
    exact = model.scattering(frequency)
    margin = float(np.max(np.linalg.eigvals(model.A).real))
    order = len(model.polynomials)-1
    states_error = float(np.max(abs(state_response(model, frequency, np.float64)-exact)))
    float32_error = float(np.max(abs(state_response(model, frequency, np.float32)-exact)))
    maximum_reflection = float(np.max(np.linalg.svd(exact, compute_uv=False)))
    reciprocal_error = float(np.max(abs(exact-np.swapaxes(exact, 1, 2))))
    n_ports = model.basis.shape[1]
    basis_error = float(np.max(abs(model.basis.T@(model.weights[:, None]*model.basis)-np.eye(n_ports))))
    # Apply only the constant aperture velocity and retain all pressure modes.
    impedance = np.linalg.solve(np.eye(n_ports)-exact[0], np.eye(n_ports)+exact[0])
    x = 2j*np.pi*frequency[0]*model.radius_m/model.sound_speed_m_s
    end_correction = float((impedance[0,0]/x).real)
    fundamentals = np.array([466.1637615, 932.327523, 1864.655046])
    harmonics = fundamentals[-1]*np.arange(1., 7.)
    angles = (0., 45.)
    receiver = model.receiver_transfer(harmonics, angles)
    hankel = hankel_receiver(model, harmonics, angles)
    hankel_error = float(np.max(np.linalg.norm(receiver-hankel, axis=2)/np.linalg.norm(hankel, axis=2)))
    rings = np.sqrt(np.r_[0., np.cumsum(model.weights)])*model.radius_m
    reference = exterior_model(rings, n_ports, 12, sound_speed_m_s=model.sound_speed_m_s,
                               density_kg_m3=model.density_kg_m3)
    reference_receiver = reference.receiver_transfer(harmonics, angles)
    order_error = np.linalg.norm(receiver-reference_receiver, axis=2)/np.linalg.norm(reference_receiver, axis=2)
    directivity = 20*np.log10(abs(receiver[:,1,0]/receiver[:,0,0]))
    return {
        'spherical_order': order,
        'aperture_ports': n_ports,
        'state_count': len(model.A),
        'max_real_pole_in_dimensionless_time': margin,
        'state_space_formula_max_error': states_error,
        'float32_state_space_formula_max_error': float32_error,
        'reflection_max_singular_value': maximum_reflection,
        'reflection_reciprocity_max_error': reciprocal_error,
        'aperture_basis_orthonormal_max_error': basis_error,
        'low_frequency_mean_port_end_correction_over_radius': end_correction,
        'receiver_hankel_max_relative_error': hankel_error,
        'angles_deg': list(angles),
        'bb6_harmonic_frequency_hz': harmonics.tolist(),
        'plane_mode_45deg_relative_to_0deg_db': directivity.tolist(),
        'receiver_operator_relative_error_vs_order12': order_error.tolist(),
        'passed': bool(margin<0 and states_error<1e-8 and float32_error<3e-6
                       and maximum_reflection<1+1e-10 and reciprocal_error<1e-10
                       and basis_error<1e-10 and hankel_error<1e-9),
    }


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, default=10)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    mesh = trumpet_mesh()
    edges = np.linspace(0., mesh['radius'][-1], mesh['rings'][-1]+1)
    result = verification(exterior_model(edges, exterior_order=args.order))
    output = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(output+'\n', encoding='utf-8')
    print(output)
    if not result['passed']:
        raise SystemExit(1)
