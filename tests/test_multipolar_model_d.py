"""Independent and regression tests for multipolar Model D (T07)."""

from itertools import product

import numpy as np
import pytest

from acoustic_ms import (
    compare_nodal_model_d,
    decompose_multipolar_cluster,
    equilateral_trimer,
    irregular_quartet,
    linear_quartet,
    linear_trimer,
    mode_index,
    modes,
    rayleigh_multipolar_scattering_coefficients,
    rms_vector_magnitude,
    scalene_trimer,
    separation_coefficient,
    solve_multipolar_nodal,
    solve_multipolar_nodal_interaction_forces,
    solve_rayleigh_nodal_interaction_forces,
    square_quartet,
)
from acoustic_ms.corrected_pair import corrected_nodal_pair_force_magnitude
from acoustic_ms.incident import nodal_standing_wave_coefficients


def _pair(distance=2.1):
    return np.array([[-distance / 2, 0.0, 0.0], [distance / 2, 0.0, 0.0]])


def _reduced_dimer_coefficients(ka, distance, f1, lmax):
    """One-particle symmetry reduction independent of the production solver."""
    positions = _pair(distance)
    selected = [
        (ell, m)
        for ell in range(1, lmax + 1)
        for m in range(ell + 1)
        if (ell + m) % 2 == 1
    ]
    scattering = rayleigh_multipolar_scattering_coefficients(ka, 0.0, f1, lmax)
    incident = nodal_standing_wave_coefficients(lmax)
    matrix = np.eye(len(selected), dtype=complex)
    rhs = np.zeros(len(selected), dtype=complex)
    for row, (ell, m) in enumerate(selected):
        rhs[row] = scattering[ell] * incident[mode_index(ell, m)]
        for column, (source_ell, source_m) in enumerate(selected):
            translated = (-1) ** source_m * separation_coefficient(
                ell, m, source_ell, source_m, ka, positions[0], positions[1]
            )
            if source_m:
                translated += separation_coefficient(
                    ell, m, source_ell, -source_m, ka, positions[0], positions[1]
                )
            matrix[row, column] -= scattering[ell] * translated
    return tuple(selected), np.linalg.solve(matrix, rhs)


def test_general_scattering_coefficients_reproduce_appendix_a():
    ka, f1 = 0.07, 0.8
    values = rayleigh_multipolar_scattering_coefficients(ka, 0.3, f1, 5)
    assert values[0] == pytest.approx(-1j * 0.3 * ka**3 / 3)
    assert values[1] == pytest.approx(1j * f1 * ka**3 / 6)
    assert values[3] == pytest.approx(1j * f1 * ka**7 / (350 * (7 - f1)))
    assert values[5] == pytest.approx(1j * f1 * ka**11 / (1309770 * (11 - 2 * f1)))


@pytest.mark.parametrize("ell", [1, 2, 3, 5, 9])
def test_scattering_coefficient_scaling(ell):
    first = rayleigh_multipolar_scattering_coefficients(0.08, 0.0, 0.6, ell)[ell]
    second = rayleigh_multipolar_scattering_coefficients(0.04, 0.0, 0.6, ell)[ell]
    assert second / first == pytest.approx(0.5 ** (2 * ell + 1))


def test_scattering_validation_and_zero_contrast():
    assert np.all(rayleigh_multipolar_scattering_coefficients(0.1, 0.0, 0.0, 5) == 0)
    for args in [(0.0, 0, 0.8, 3), (0.11, 0, 0.8, 3), (0.1, 0, 1.1, 3), (0.1, 0, 0.8, 0)]:
        with pytest.raises(ValueError):
            rayleigh_multipolar_scattering_coefficients(*args)


def test_complete_mode_ordering_and_planar_rule():
    result = solve_multipolar_nodal(_pair(), 0.1, 1.0, 0.0, 0.8, 5)
    assert result.modes == modes(5)
    assert all((ell + m) % 2 == 1 for ell, m in result.active_modes)
    inactive = set(range(36)) - set(result.active_mode_indices)
    assert np.count_nonzero(result.coefficients[:, sorted(inactive)]) == 0


def test_active_and_complete_bases_agree_for_generic_cluster():
    positions = scalene_trimer(2.1)
    active = solve_multipolar_nodal(positions, 0.1, 1, 0, 0.8, 3)
    complete = solve_multipolar_nodal(
        positions, 0.1, 1, 0, 0.8, 3, use_planar_symmetry=False
    )
    np.testing.assert_allclose(active.coefficients, complete.coefficients, rtol=2e-12, atol=2e-16)


def test_balanced_solution_matches_raw_and_has_physical_residual():
    result = solve_multipolar_nodal(scalene_trimer(2.1), 0.1, 1, 0, 0.8, 3)
    raw = np.linalg.solve(result.physical_system_matrix, result.physical_right_hand_side)
    active = result.coefficients[:, np.asarray(result.active_mode_indices)].ravel()
    np.testing.assert_allclose(active, raw, rtol=2e-12, atol=2e-16)
    assert result.residual_relative < 2e-14
    assert result.condition_number < result.physical_condition_number


@pytest.mark.parametrize("lmax", [1, 3, 5])
def test_reduced_dimer_oracle(lmax):
    selected, expected = _reduced_dimer_coefficients(0.1, 2.1, 0.8, lmax)
    actual = solve_multipolar_nodal(_pair(), 0.1, 1, 0, 0.8, lmax)
    observed = np.array([actual.coefficients[0, mode_index(*mode)] for mode in selected])
    np.testing.assert_allclose(observed, expected, rtol=3e-12, atol=3e-16)


def test_model_d_l1_reproduces_model_c_for_canonical_geometries():
    geometries = [
        _pair(), np.array([[-0.8, -1.05, 0], [0.8, 1.05, 0]]),
        linear_trimer(2.1), equilateral_trimer(2.1), scalene_trimer(2.1),
        linear_quartet(2.1), square_quartet(2.1), irregular_quartet(2.1),
    ]
    for positions in geometries:
        new = solve_multipolar_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 0.8, 1)
        old = solve_rayleigh_nodal_interaction_forces(positions, 0.1, 1, 1, 0, 0.8)
        np.testing.assert_allclose(new.forces_xy, old.forces_xy, rtol=3e-13, atol=3e-15)
        np.testing.assert_allclose(new.solution.coefficients, old.solution.coefficients, rtol=3e-13, atol=3e-16)


def test_general_planar_difference_from_eq30_is_separately_resolved():
    differences = []
    for ka in [0.1, 0.05, 0.025]:
        result = solve_multipolar_nodal_interaction_forces(_pair(2.1), ka, 1, 1, 0, 0.8, 5)
        radial = float(result.forces_xy[0, 0])
        analytic = -corrected_nodal_pair_force_magnitude(ka, 1, 2.1, 1, 0.8)
        differences.append(abs(radial - analytic) / abs(radial))
    assert differences[0] > differences[1] > differences[2]
    assert differences[2] < differences[0]


def test_dimer_invariances_energy_and_attraction():
    base = _pair(2.1)
    reference = solve_multipolar_nodal_interaction_forces(base, 0.1, 1, 1, 0, 0.8, 5).forces_xy
    assert reference[0, 0] > 0 and reference[1, 0] < 0
    shifted = solve_multipolar_nodal_interaction_forces(base + [2.3, -1.4, 0], 0.1, 1, 1, 0, 0.8, 5).forces_xy
    np.testing.assert_allclose(shifted, reference, rtol=2e-12, atol=2e-14)
    angle = 0.37
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rotated_positions = base.copy(); rotated_positions[:, :2] = base[:, :2] @ rotation.T
    rotated = solve_multipolar_nodal_interaction_forces(rotated_positions, 0.1, 1, 1, 0, 0.8, 5).forces_xy
    np.testing.assert_allclose(rotated, reference @ rotation.T, rtol=3e-12, atol=3e-14)
    doubled = solve_multipolar_nodal_interaction_forces(base, 0.1, 1, 2, 0, 0.8, 5).forces_xy
    np.testing.assert_allclose(doubled, 2 * reference, rtol=2e-14, atol=2e-14)
    permuted = solve_multipolar_nodal_interaction_forces(base[::-1], 0.1, 1, 1, 0, 0.8, 5).forces_xy
    np.testing.assert_allclose(permuted, reference[::-1], rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(reference.sum(axis=0), 0, atol=3e-14)
    with pytest.raises(ValueError):
        solve_multipolar_nodal_interaction_forces(_pair(1.99), 0.1, 1, 1, 0, 0.8, 3)


def test_zero_contrast_and_model_comparison():
    zero = compare_nodal_model_d(scalene_trimer(2.1), 0.1, 1, 1, 0, 0, 5)
    assert np.count_nonzero(zero.model_d_forces_xy) == 0
    comparison = compare_nodal_model_d(scalene_trimer(2.1), 0.1, 1, 1, 0, 0.8, 3)
    np.testing.assert_allclose(
        comparison.model_d_forces_xy - comparison.model_c_forces_xy,
        comparison.multipolar_correction_xy,
    )


@pytest.mark.parametrize("positions", [scalene_trimer(2.1), irregular_quartet(2.1)])
def test_connected_inclusion_exclusion_and_l1_regression(positions):
    expansion = decompose_multipolar_cluster(positions, 0.1, 1, 1, 0, 0.8, 1)
    reconstructed = expansion.two_body_sum_xy + expansion.irreducible_three_body_sum_xy
    if len(positions) == 4:
        reconstructed += expansion.irreducible_four_body_xy
    np.testing.assert_allclose(reconstructed, expansion.model_d_forces_xy, rtol=2e-13, atol=2e-14)
    if len(positions) == 3:
        old = compare_nodal_model_d(positions, 0.1, 1, 1, 0, 0.8, 1)
        np.testing.assert_allclose(expansion.irreducible_three_body_sum_xy, old.collective_rayleigh_correction_xy, rtol=3e-13, atol=2e-14)


def test_high_order_translation_affects_local_force_finitely():
    low = solve_multipolar_nodal_interaction_forces(_pair(), 0.1, 1, 1, 0, 0.8, 1)
    high = solve_multipolar_nodal_interaction_forces(_pair(), 0.1, 1, 1, 0, 0.8, 5)
    assert np.all(np.isfinite(high.local_scattered_coefficients))
    assert rms_vector_magnitude(high.forces_xy - low.forces_xy) > 0
