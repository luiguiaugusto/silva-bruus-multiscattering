"""Tests for the exact-Mie globally coupled effective-incident solver."""

import numpy as np
import pytest

from acoustic_ms import (
    irregular_quartet,
    linear_trimer,
    nodal_standing_wave_coefficients,
    solve_mie_multiparticle_nodal,
    translation_matrix,
)


def test_mie_solver_effective_and_scattered_systems_are_equivalent():
    positions = irregular_quartet(2.3)
    solution = solve_mie_multiparticle_nodal(
        positions, 0.1, 1.0, 0.2, 0.8, 3, use_planar_symmetry=False
    )
    count = len(solution.modes)
    translation = np.zeros((4 * count, 4 * count), dtype=complex)
    for target in range(4):
        for source in range(4):
            if target != source:
                rows = slice(target * count, (target + 1) * count)
                columns = slice(source * count, (source + 1) * count)
                translation[rows, columns] = translation_matrix(
                    0.1, positions[target], positions[source], 3
                )
    diagonal = np.tile(
        [solution.scattering_coefficients[ell] for ell, _ in solution.modes], 4
    )
    external = np.tile(nodal_standing_wave_coefficients(3), 4)
    scattered_alternate = np.linalg.solve(
        np.eye(4 * count) - diagonal[:, None] * translation,
        diagonal * external,
    )
    np.testing.assert_allclose(
        scattered_alternate,
        solution.scattered_coefficients.ravel(),
        rtol=3e-12,
        atol=3e-13,
    )


def test_single_particle_has_external_effective_field_and_diagonal_scattering():
    solution = solve_mie_multiparticle_nodal(
        np.zeros((1, 3)), 0.1, 1.0, 0.1, 0.8, 4
    )
    np.testing.assert_array_equal(
        solution.effective_incident_coefficients,
        solution.external_incident_coefficients,
    )
    expected = np.array([
        solution.scattering_coefficients[ell] * coefficient
        for (ell, _), coefficient in zip(
            solution.modes, solution.external_incident_coefficients[0]
        )
    ])
    np.testing.assert_allclose(solution.scattered_coefficients[0], expected)


def test_effective_incident_identity_is_satisfied_directly():
    positions = linear_trimer(2.4)
    solution = solve_mie_multiparticle_nodal(
        positions, 0.1, 1.0, 0.0, 1.0, 3, use_planar_symmetry=False
    )
    reconstructed = np.array(solution.external_incident_coefficients, copy=True)
    for target in range(3):
        for source in range(3):
            if target != source:
                reconstructed[target] += translation_matrix(
                    0.1, positions[target], positions[source], 3
                ) @ solution.scattered_coefficients[source]
    np.testing.assert_allclose(
        solution.effective_incident_coefficients,
        reconstructed,
        rtol=3e-12,
        atol=3e-13,
    )


@pytest.mark.parametrize("positions", [linear_trimer(2.2), irregular_quartet(2.2)])
def test_planar_active_base_matches_complete_base(positions):
    reduced = solve_mie_multiparticle_nodal(positions, 0.1, 1.0, 0.0, 0.8, 3)
    complete = solve_mie_multiparticle_nodal(
        positions, 0.1, 1.0, 0.0, 0.8, 3, use_planar_symmetry=False
    )
    np.testing.assert_allclose(
        reduced.effective_incident_coefficients,
        complete.effective_incident_coefficients,
        rtol=2e-12,
        atol=2e-13,
    )
    np.testing.assert_allclose(
        reduced.scattered_coefficients,
        complete.scattered_coefficients,
        rtol=2e-12,
        atol=2e-13,
    )
    for index, (ell, m) in enumerate(reduced.modes):
        if (ell + m) % 2 == 0:
            assert np.all(reduced.effective_incident_coefficients[:, index] == 0.0)
            assert np.all(reduced.scattered_coefficients[:, index] == 0.0)


def test_material_matching_is_handled_exactly():
    solution = solve_mie_multiparticle_nodal(
        linear_trimer(2.1), 0.3, 1.0, 0.0, 0.0, 4
    )
    assert np.all(solution.scattering_coefficients == 0.0)
    assert np.all(solution.scattered_coefficients == 0.0)
    np.testing.assert_array_equal(
        solution.effective_incident_coefficients,
        solution.external_incident_coefficients,
    )
    assert solution.residual_relative == 0.0


@pytest.mark.parametrize(
    "positions,k,radius,lmax,error",
    [
        (np.zeros((2, 2)), 0.1, 1.0, 2, ValueError),
        (np.array([[0.0, 0.0, 0.1]]), 0.1, 1.0, 2, ValueError),
        (np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0]]), 0.1, 1.0, 2, ValueError),
        (np.zeros((1, 3)), 0.0, 1.0, 2, ValueError),
        (np.zeros((1, 3)), 0.1, -1.0, 2, ValueError),
        (np.zeros((1, 3)), 0.1, 1.0, 0, ValueError),
        (np.zeros((1, 3)), 0.1, 1.0, 2.5, TypeError),
    ],
)
def test_mie_solver_rejects_invalid_inputs(positions, k, radius, lmax, error):
    with pytest.raises(error):
        solve_mie_multiparticle_nodal(positions, k, radius, 0.0, 0.8, lmax)


def test_mie_solver_accepts_generic_size_parameter_beyond_rayleigh_domain():
    solution = solve_mie_multiparticle_nodal(
        np.array([[-1.1, 0.0, 0.0], [1.1, 0.0, 0.0]]),
        0.5,
        1.0,
        0.0,
        0.8,
        3,
    )
    assert np.all(np.isfinite(solution.scattered_coefficients))
    assert solution.residual_relative < 1e-11
