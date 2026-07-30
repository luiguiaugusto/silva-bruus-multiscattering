"""Numerical-stability regressions for the T11.1 balanced Model-E solver."""

from pathlib import Path

import mpmath as mp
import numpy as np

from acoustic_ms import (
    scalene_trimer,
    solve_mie_multiparticle_nodal,
    solve_model_e_nodal,
)


ROOT = Path(__file__).resolve().parents[1]


def _active(solution, full):
    indices = np.concatenate([
        particle * len(solution.modes) + np.asarray(solution.active_mode_indices)
        for particle in range(len(full))
    ])
    return np.asarray(full).ravel()[indices]


def _relative(first, second):
    scale = max(float(np.linalg.norm(first)), float(np.linalg.norm(second)))
    if scale == 0.0:
        return float(np.linalg.norm(np.asarray(first) - np.asarray(second)))
    return float(np.linalg.norm(np.asarray(first) - np.asarray(second))) / scale


def test_balanced_matrices_and_legacy_aliases_are_explicit():
    positions = np.array([[-1.2, 0.0, 0.0], [1.2, 0.0, 0.0]])
    solution = solve_mie_multiparticle_nodal(
        positions, 0.1, 1.0, 0.0, 0.8, 4
    )
    identity = np.eye(len(solution.balanced_coefficients))
    diagonal = solution.scattering_diagonal
    square_root = solution.square_root_scattering_diagonal
    translation = solution.translation_matrix
    expected_b = identity - translation * diagonal[None, :]
    expected_d = identity - diagonal[:, None] * translation
    expected_q = (
        identity
        - square_root[:, None] * translation * square_root[None, :]
    )
    np.testing.assert_array_equal(
        solution.system_matrix, solution.effective_incident_system_matrix
    )
    np.testing.assert_array_equal(
        solution.right_hand_side, solution.effective_incident_right_hand_side
    )
    np.testing.assert_allclose(solution.system_matrix, expected_b)
    np.testing.assert_allclose(solution.scattered_system_matrix, expected_d)
    np.testing.assert_allclose(solution.balanced_system_matrix, expected_q)
    np.testing.assert_allclose(
        solution.scattered_right_hand_side,
        diagonal * solution.right_hand_side,
    )
    np.testing.assert_allclose(
        solution.balanced_right_hand_side,
        square_root * solution.right_hand_side,
    )
    np.testing.assert_array_equal(square_root, np.sqrt(diagonal))
    assert solution.condition_number == np.linalg.cond(solution.system_matrix)
    effective = _active(solution, solution.effective_incident_coefficients)
    expected_residual = (
        np.linalg.norm(solution.system_matrix @ effective - solution.right_hand_side)
        / max(np.linalg.norm(solution.right_hand_side), np.finfo(float).eps)
    )
    assert solution.residual_relative == expected_residual

    assert solution.production_solver == "balanced_sqrt"


def test_balanced_reconstruction_and_physical_closures():
    solution = solve_mie_multiparticle_nodal(
        scalene_trimer(2.7), 0.1, 1.0, 0.0, 0.8, 5
    )
    q = solution.balanced_coefficients
    d = solution.square_root_scattering_diagonal * q
    b = solution.right_hand_side + solution.translation_matrix @ d
    np.testing.assert_allclose(
        d, _active(solution, solution.scattered_coefficients),
        rtol=2e-13, atol=2e-15,
    )
    np.testing.assert_allclose(
        b, _active(solution, solution.effective_incident_coefficients),
        rtol=2e-13, atol=2e-15,
    )
    assert solution.balanced_backward_error < 1e-12
    assert solution.effective_incident_closure_error < 1e-12
    assert solution.scattering_closure_error < 1e-12


def test_three_formulations_agree_at_safe_order():
    solution = solve_mie_multiparticle_nodal(
        scalene_trimer(2.7), 0.1, 1.0, 0.0, 0.8, 3
    )
    legacy_b = np.linalg.solve(solution.system_matrix, solution.right_hand_side)
    scattered_d = np.linalg.solve(
        solution.scattered_system_matrix, solution.scattered_right_hand_side
    )
    scattered_b = solution.right_hand_side + solution.translation_matrix @ scattered_d
    balanced_b = _active(solution, solution.effective_incident_coefficients)
    balanced_d = _active(solution, solution.scattered_coefficients)
    assert _relative(balanced_b, legacy_b) < 1e-11
    assert _relative(balanced_b, scattered_b) < 1e-12
    assert _relative(balanced_d, scattered_d) < 1e-12


def test_high_order_conditioning_is_fixed_without_changing_legacy_diagnostic():
    solution = solve_mie_multiparticle_nodal(
        scalene_trimer(2.7), 0.1, 1.0, 0.0, 0.8, 9
    )
    assert solution.condition_number > 1e20
    assert solution.balanced_condition_number < 10.0
    assert solution.balanced_backward_error < 1e-12
    assert solution.effective_incident_closure_error < 1e-12
    assert solution.scattering_closure_error < 1e-12


def test_transparent_material_has_exact_zero_balanced_unknowns():
    solution = solve_mie_multiparticle_nodal(
        scalene_trimer(2.7), 0.1, 1.0, 0.0, 0.0, 6
    )
    assert np.all(solution.square_root_scattering_diagonal == 0.0)
    assert np.all(solution.balanced_coefficients == 0.0)
    assert np.all(solution.scattered_coefficients == 0.0)
    np.testing.assert_array_equal(
        solution.effective_incident_coefficients,
        solution.external_incident_coefficients,
    )
    assert solution.balanced_backward_error == 0.0
    assert solution.effective_incident_closure_error == 0.0
    assert solution.scattering_closure_error == 0.0


def test_small_high_precision_balanced_linear_oracle():
    solution = solve_mie_multiparticle_nodal(
        np.array([[-1.25, 0.0, 0.0], [1.25, 0.0, 0.0]]),
        0.1, 1.0, 0.0, 0.8, 3,
    )
    mp.mp.dps = 70
    matrix = mp.matrix([
        [mp.mpc(float(value.real), float(value.imag)) for value in row]
        for row in solution.balanced_system_matrix
    ])
    rhs = mp.matrix([
        mp.mpc(float(value.real), float(value.imag))
        for value in solution.balanced_right_hand_side
    ])
    reference = mp.lu_solve(matrix, rhs)
    q_reference = np.array([complex(reference[index]) for index in range(len(reference))])
    assert _relative(solution.balanced_coefficients, q_reference) < 1e-12


def test_production_source_uses_only_direct_balanced_solve():
    source = (ROOT / "src" / "acoustic_ms" / "mie_multiparticle.py").read_text(
        encoding="utf-8"
    )
    assert "np.linalg.solve(balanced_system, balanced_rhs)" in source
    assert "np.linalg.solve(system, rhs)" not in source
    for forbidden in ("np.linalg.inv", "np.linalg.pinv", "np.linalg.lstsq"):
        assert forbidden not in source


def test_model_e_force_channels_remain_finite_under_balanced_solver():
    result = solve_model_e_nodal(
        scalene_trimer(2.7), 0.1, 1.0, 1.0, 0.0, 0.8, 9
    )
    for values in (
        result.total_forces_xyz,
        result.interaction_forces_xyz,
        result.external_scattered_forces_xyz,
        result.scattered_scattered_forces_xyz,
    ):
        assert np.all(np.isfinite(values))
    assert result.decomposition_residual < 1e-15
