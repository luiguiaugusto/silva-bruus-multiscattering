"""Scientific invariance and convergence tests for exact-Mie Model E."""

import numpy as np
import pytest

from acoustic_ms import (
    equilateral_trimer,
    irregular_quartet,
    solve_model_e_nodal,
)


def _rotate(vectors, angle):
    matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    result = np.array(vectors, copy=True)
    result[:, :2] = result[:, :2] @ matrix.T
    return result


def test_model_e_force_channels_reconstruct_total_force():
    result = solve_model_e_nodal(irregular_quartet(2.2), 0.1, 1.0, 1.0, 0.1, 0.8, 4)
    np.testing.assert_allclose(
        result.interaction_forces_xyz,
        result.external_scattered_forces_xyz + result.scattered_scattered_forces_xyz,
        rtol=2e-15,
        atol=2e-15,
    )
    np.testing.assert_allclose(
        result.total_forces_xyz,
        result.external_forces_xyz + result.external_scattered_forces_xyz + result.scattered_scattered_forces_xyz,
        rtol=2e-15,
        atol=2e-15,
    )
    assert result.decomposition_residual < 1e-15


def test_model_e_translation_rotation_and_permutation_covariance():
    positions = irregular_quartet(2.3)
    base = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.8, 3)
    shift = np.array([3.2, -1.7, 0.0])
    translated = solve_model_e_nodal(positions + shift, 0.1, 1.0, 1.0, 0.0, 0.8, 3)
    np.testing.assert_allclose(translated.total_forces_xyz, base.total_forces_xyz, rtol=3e-12, atol=3e-13)
    angle = 0.37
    rotated_positions = _rotate(positions, angle)
    rotated = solve_model_e_nodal(rotated_positions, 0.1, 1.0, 1.0, 0.0, 0.8, 3)
    np.testing.assert_allclose(rotated.total_forces_xyz, _rotate(base.total_forces_xyz, angle), rtol=3e-12, atol=3e-13)
    order = np.array([3, 1, 0, 2])
    permuted = solve_model_e_nodal(positions[order], 0.1, 1.0, 1.0, 0.0, 0.8, 3)
    np.testing.assert_allclose(permuted.total_forces_xyz, base.total_forces_xyz[order], rtol=3e-12, atol=3e-13)


def test_model_e_energy_scaling_and_material_matching():
    positions = equilateral_trimer(2.2)
    base = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.8, 3)
    scaled = solve_model_e_nodal(positions, 0.1, 1.0, 2.5, 0.0, 0.8, 3)
    np.testing.assert_allclose(scaled.total_forces_xyz, 2.5 * base.total_forces_xyz, rtol=2e-14, atol=2e-14)
    zero_energy = solve_model_e_nodal(positions, 0.1, 1.0, 0.0, 0.0, 0.8, 3)
    assert np.all(zero_energy.total_forces_xyz == 0.0)
    matched = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.0, 3)
    assert np.all(matched.total_forces_xyz == 0.0)


def test_model_e_pair_action_reaction_and_attraction():
    positions = np.array([[-1.05, 0.0, 0.0], [1.05, 0.0, 0.0]])
    result = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.8, 5)
    np.testing.assert_allclose(result.total_forces_xyz[0], -result.total_forces_xyz[1], rtol=3e-12, atol=3e-13)
    assert result.interaction_forces_xyz[0, 0] > 0.0
    assert result.interaction_forces_xyz[1, 0] < 0.0


def test_model_e_converges_across_successive_orders():
    positions = irregular_quartet(4.0)
    forces = [
        solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.8, order).total_forces_xyz
        for order in range(2, 8)
    ]
    errors = [
        np.linalg.norm(forces[index] - forces[index - 1]) / np.linalg.norm(forces[index])
        for index in range(1, len(forces))
    ]
    assert errors[-1] < errors[0]
    assert errors[-1] < 1e-5


def test_model_e_requires_force_resolving_order():
    with pytest.raises(ValueError, match="lmax >= 2"):
        solve_model_e_nodal(np.zeros((1, 3)), 0.1, 1.0, 1.0, 0.0, 0.8, 1)
