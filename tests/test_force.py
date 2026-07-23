"""Independent validation for the T04 nodal interaction-force observable."""

import numpy as np
import pytest

from acoustic_ms import (
    nodal_pair_force_magnitude,
    corrected_nodal_pair_force_magnitude,
    solve_rayleigh_nodal_interaction_forces,
    mode_index,
)
from acoustic_ms.scattering import rayleigh_scattering_coefficients
from acoustic_ms.special import spherical_hankel1


def _cartesian_oracle(positions, k, radius, energy, f1, s10):
    """Sec. 7 Cartesian oracle; it does not use translations or production force."""
    positions = np.asarray(positions, dtype=float)
    output = np.zeros((len(positions), 2))
    factor = 4 * np.pi * k * radius**3 * energy * np.sqrt(3 / (4 * np.pi))
    for target in range(len(positions)):
        for source in range(len(positions)):
            if target != source:
                displacement = positions[target] - positions[source]
                distance = np.linalg.norm(displacement)
                x = k * distance
                q = spherical_hankel1(1, x, derivative=True) / x - spherical_hankel1(1, x) / x**2
                output[target] += factor * np.real(np.conj(f1) * s10[source] * q * displacement[:2] / distance)
    return output


def _scalar_pair(ka, separation, f1, energy=1.0):
    s1 = 1j * f1 * ka**3 / 6
    a10 = np.sqrt(12 * np.pi)
    kd = ka * separation
    coupling = spherical_hankel1(0, kd) + spherical_hankel1(2, kd)
    s10 = s1 * a10 / (1 - s1 * coupling)
    q = spherical_hankel1(1, kd, derivative=True) / kd - spherical_hankel1(1, kd) / kd**2
    radial = 4 * np.pi * ka * energy * np.sqrt(3 / (4 * np.pi)) * np.real(np.conj(f1) * s10 * q)
    return s10, radial


def _pair_positions(separation, angle=0.0):
    axis = np.array([np.cos(angle), np.sin(angle), 0.0])
    return np.array([-separation / 2 * axis, separation / 2 * axis])


def test_single_particle_has_zero_interaction_field_and_force():
    result = solve_rayleigh_nodal_interaction_forces([[0, 0, 0]], .08, 1, 2, -.3, .7)
    assert np.allclose(result.local_scattered_coefficients, 0)
    assert np.allclose(result.forces_xy, 0)
    assert np.count_nonzero(result.solution.coefficients) == 1


@pytest.mark.parametrize("ka,separation,f1,angle", [(.1, 2., .8, 0.), (.07, 3.4, .4, .61), (.03, 8., 1., 1.1)])
def test_force_matches_independent_cartesian_oracle(ka, separation, f1, angle):
    positions = _pair_positions(separation, angle)
    result = solve_rayleigh_nodal_interaction_forces(positions, ka, 1, 1.7, -.5, f1)
    expected = _cartesian_oracle(positions, ka, 1, 1.7, f1, result.solution.coefficients[:, 2])
    assert np.allclose(result.forces_xy, expected, rtol=2e-12, atol=2e-13)


@pytest.mark.parametrize("f1,expected", [(.1, -.011936371917121), (.4, -.194729303800953), (.8, -.799842697325624), (1., -1.26676999261163)])
def test_scalar_pair_solution_and_contact_regression(f1, expected):
    positions = _pair_positions(2.)
    result = solve_rayleigh_nodal_interaction_forces(positions, .1, 1, 1, .2, f1)
    s10, radial = _scalar_pair(.1, 2., f1)
    direction = (positions[0] - positions[1])[:2] / 2
    assert np.allclose(result.solution.coefficients[:, 2], s10, rtol=2e-12, atol=2e-14)
    assert np.isclose(np.dot(result.forces_xy[0], direction), radial, rtol=2e-12)
    assert np.isclose(radial, expected, rtol=2e-12)
    assert radial < 0


def test_eq_28_planar_pair_and_local_symmetry():
    ka, separation, f1 = .1, 2., .8
    result = solve_rayleigh_nodal_interaction_forces(_pair_positions(separation), ka, 1, 1, 0., f1)
    b21 = result.local_scattered_coefficients[0, 7]
    expected_x = -2 * np.sqrt(30 * np.pi) / (15 * ka**2) * ka**3 * np.real(np.conj(f1) * b21)
    assert np.isclose(result.forces_xy[0, 0], expected_x, rtol=2e-12)
    assert np.allclose(result.forces_xy[1], -result.forces_xy[0], rtol=2e-12, atol=2e-14)
    assert np.allclose(result.local_scattered_coefficients[0, 5], -b21, rtol=2e-12)


def test_pair_symmetries_and_physical_scales():
    positions = _pair_positions(3., .4)
    result = solve_rayleigh_nodal_interaction_forces(positions, .08, 1, 1.2, -.5, .8)
    assert np.allclose(result.forces_xy[1], -result.forces_xy[0], rtol=2e-12, atol=2e-14)
    direction = positions[0, :2] - positions[1, :2]
    assert abs(direction[0] * result.forces_xy[0, 1] - direction[1] * result.forces_xy[0, 0]) < 2e-13
    shift = np.array([3., -2., 0.])
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions + shift, .08, 1, 1.2, -.5, .8).forces_xy, result.forces_xy)
    angle = .8; rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions @ rotation.T, .08, 1, 1.2, -.5, .8).forces_xy, result.forces_xy @ rotation[:2, :2].T)
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions, .08, 1, 2.4, -.5, .8).forces_xy, 2 * result.forces_xy)
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions, .08, 1, 1.2, 4., .8).forces_xy, result.forces_xy)


def test_zeros_radius_scale_and_invalid_energy():
    positions = _pair_positions(2.)
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions, .1, 1, 0, 0, .8).forces_xy, 0)
    assert np.allclose(solve_rayleigh_nodal_interaction_forces(positions, .1, 1, 1, 0, 0).forces_xy, 0)
    scaled = solve_rayleigh_nodal_interaction_forces(2 * positions, .05, 2, 1, 0, .8).forces_xy
    base = solve_rayleigh_nodal_interaction_forces(positions, .1, 1, 1, 0, .8).forces_xy
    assert np.allclose(scaled, 4 * base, rtol=2e-12)
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(positions, .1, 1, -1, 0, .8)


def test_silva_bruus_recovery_limits_and_attraction():
    weak = solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), .1, 1, 1, 0, 1e-7).forces_xy[0, 0]
    weak_sb = nodal_pair_force_magnitude(.1, 1, 2, 1, 1e-7)
    assert np.isclose(weak, -weak_sb, rtol=1e-8)
    far = solve_rayleigh_nodal_interaction_forces(_pair_positions(1000.), .1, 1, 1, 0, .8).forces_xy[0, 0]
    far_sb = nodal_pair_force_magnitude(.1, 1, 1000, 1, .8)
    assert np.isclose(far, -far_sb, rtol=2e-6)
    near = solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), 1e-4, 1, 1, 0, .8).forces_xy[0, 0]
    assert near > 0


@pytest.mark.parametrize("positions,k,radius,energy,f1", [
    ([[0., 0., 0.]], 0., 1., 1., .8),
    ([[0., 0., 0.]], .1, 0., 1., .8),
    ([[0., 0., .1]], .1, 1., 1., .8),
    ([[0., 0., 0.], [1., 0., 0.]], .1, 1., 1., .8),
    ([[0., 0., 0.]], .1, 1., 1., 1.1),
])
def test_force_api_reuses_solver_domain_validation(positions, k, radius, energy, f1):
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(positions, k, radius, energy, 0., f1)


def test_oblique_pair_particle_permutation_preserves_all_particle_observables():
    positions = _pair_positions(3.7, .53)
    original = solve_rayleigh_nodal_interaction_forces(positions, .08, 1., 1.2, -.4, .7)
    order = np.array([1, 0])
    permuted = solve_rayleigh_nodal_interaction_forces(positions[order], .08, 1., 1.2, -.4, .7)
    assert np.allclose(permuted.forces_xy, original.forces_xy[order], rtol=2e-12, atol=2e-14)
    assert np.allclose(permuted.local_scattered_coefficients, original.local_scattered_coefficients[order], rtol=2e-12, atol=2e-14)
    assert np.allclose(permuted.solution.coefficients, original.solution.coefficients[order], rtol=2e-12, atol=2e-14)


def test_planar_aligned_pair_has_only_allowed_local_modes():
    result = solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), .1, 1., 1., 0., .8)
    forbidden = tuple(mode_index(ell, m) for ell, m in ((0, 0), (1, -1), (1, 1), (2, -2), (2, 0), (2, 2)))
    assert np.max(abs(result.local_scattered_coefficients[:, forbidden])) < 2e-13
    assert np.allclose(result.local_scattered_coefficients[:, mode_index(2, -1)], -result.local_scattered_coefficients[:, mode_index(2, 1)], rtol=2e-12, atol=2e-14)


@pytest.mark.parametrize("energy", (np.nan, np.inf, -np.inf))
def test_nonfinite_energy_density_is_rejected_early(energy):
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), .1, 1., energy, 0., .8)


def test_lower_f1_domain_limit_is_rejected():
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), .1, 1., 1., 0., -2.01)


@pytest.mark.parametrize("lmax", (0, 2))
def test_force_api_rejects_unsupported_scattering_lmax(lmax):
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(_pair_positions(2.), .1, 1., 1., 0., .8, lmax=lmax)


@pytest.mark.parametrize("positions", ([0., 0., 0.], [[0., 0.]]))
def test_force_api_rejects_invalid_position_shape(positions):
    with pytest.raises(ValueError):
        solve_rayleigh_nodal_interaction_forces(positions, .1, 1., 1., 0., .8)
