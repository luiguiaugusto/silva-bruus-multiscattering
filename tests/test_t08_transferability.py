"""Scientific API tests for the T08 matched baseline and predictors."""

import numpy as np

from acoustic_ms.cluster_families import compact_cluster, irregular_cluster
from acoustic_ms.comparison import compare_nodal_force_models
from acoustic_ms.model_d import solve_multipolar_nodal_interaction_forces
from acoustic_ms.scaling import maximum_geometric_coupling
from acoustic_ms.transferability import (
    matched_multipolar_pairwise_baseline,
    normalized_rms_difference,
    spectral_radius_l1,
    two_step_converged,
)


def test_b1_reproduces_historical_model_b():
    positions = irregular_cluster(4, 2.5)
    historical = compare_nodal_force_models(positions, 0.1, 1, 1, 0, 0.8)
    matched = matched_multipolar_pairwise_baseline(positions, 0.1, 1, 1, 0, 0.8, 1)
    np.testing.assert_allclose(matched.forces_xy, historical.model_b_forces_xy, rtol=4e-13, atol=4e-14)


def test_matched_baseline_equals_global_dimer():
    positions = np.array([[-1.05, 0, 0], [1.05, 0, 0]])
    for lmax in (1, 3, 5):
        matched = matched_multipolar_pairwise_baseline(positions, 0.1, 1, 1, 0, 1.0, lmax)
        global_result = solve_multipolar_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 1.0, lmax)
        np.testing.assert_allclose(matched.forces_xy, global_result.forces_xy, rtol=4e-13, atol=4e-14)


def test_vector_identity_closes():
    positions = compact_cluster(3, 2.5)
    base = compare_nodal_force_models(positions, 0.1, 1, 1, 0, 0.8)
    matched = matched_multipolar_pairwise_baseline(positions, 0.1, 1, 1, 0, 0.8, 5)
    model_d = solve_multipolar_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 0.8, 5)
    np.testing.assert_allclose(
        model_d.forces_xy - base.model_a_forces_xy,
        (matched.forces_xy - base.model_a_forces_xy) + (model_d.forces_xy - matched.forces_xy),
        atol=1e-16,
    )


def test_spectral_radius_uses_balanced_coupling_operator():
    positions = irregular_cluster(6, 2.5)
    result = solve_multipolar_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 0.8, 1)
    coupling = np.eye(len(result.solution.system_matrix)) - result.solution.system_matrix
    expected = np.max(np.abs(np.linalg.eigvals(coupling)))
    assert spectral_radius_l1(result) == expected
    assert not np.isclose(expected, np.max(np.abs(np.linalg.eigvals(result.solution.system_matrix))))


def test_predictor_invariances_and_inverse_cube_scaling():
    positions = irregular_cluster(6, 2.5)
    angle = 0.43
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    transformed = positions.copy(); transformed[:, :2] = positions[:, :2] @ rotation.T + [1.3, -2.2]
    transformed = transformed[[4, 0, 5, 2, 1, 3]]
    first_lambda = maximum_geometric_coupling(positions, 1, 0.8)
    np.testing.assert_allclose(
        maximum_geometric_coupling(transformed, 1, 0.8), first_lambda, rtol=3e-15
    )
    scaled = irregular_cluster(6, 5.0)
    np.testing.assert_allclose(
        maximum_geometric_coupling(scaled, 1, 0.8), first_lambda / 8, rtol=3e-15
    )
    first_rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 0.8, 1))
    second_rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(transformed, 0.1, 1, 1, 0, 0.8, 1))
    np.testing.assert_allclose(second_rho, first_rho, rtol=3e-13, atol=3e-15)


def test_two_step_protocol_and_scale_aware_metric():
    assert two_step_converged([1e-2, 8e-4, 7e-4])
    assert not two_step_converged([8e-4])
    assert not two_step_converged([8e-4, 2e-3])
    value, applicable = normalized_rms_difference(np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2)))
    assert value == 0 and not applicable


def test_matched_baseline_is_covariant_under_rotation_and_permutation():
    positions = irregular_cluster(4, 2.5)
    angle = 0.37
    rotation = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)],
    ])
    order = np.array([3, 1, 0, 2])
    transformed = positions.copy()
    transformed[:, :2] = positions[:, :2] @ rotation.T + [1.2, -0.7]
    transformed = transformed[order]

    original = matched_multipolar_pairwise_baseline(
        positions, 0.1, 1, 1, 0, 0.8, 3
    )
    observed = matched_multipolar_pairwise_baseline(
        transformed, 0.1, 1, 1, 0, 0.8, 3
    )
    expected = original.forces_xy[order] @ rotation.T
    np.testing.assert_allclose(observed.forces_xy, expected, rtol=4e-13, atol=4e-14)
