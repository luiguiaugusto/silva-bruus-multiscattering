"""Analytical-operator and Neumann-series checks for T09."""

import numpy as np
import pytest

from acoustic_ms.cluster_families import compact_cluster, irregular_cluster
from acoustic_ms.comparison import compare_nodal_force_models
from acoustic_ms.incident import nodal_standing_wave_coefficients
from acoustic_ms.multipoles import mode_index
from acoustic_ms.multipolar_solver import solve_multipolar_nodal
from acoustic_ms.rho_foundation import (
    dipolar_balanced_coupling_matrix,
    dipolar_coupling_diagnostics,
    dipolar_coupling_entry,
    near_field_dipolar_coupling_entry,
    neumann_partial_solutions,
)
from acoustic_ms.scattering import rayleigh_scattering_coefficients
from acoustic_ms.special import spherical_hankel1
from acoustic_ms.translation import translation_matrix


@pytest.mark.parametrize("distance", (2.1, 4.0, 10.0))
@pytest.mark.parametrize("f1", (-0.7, 0.1, 1.0))
def test_closed_dipolar_entry_matches_production_hankel_expression(distance, f1):
    k, radius = 0.1, 1.0
    s1 = rayleigh_scattering_coefficients(k * radius, 0.0, f1)[1]
    expected = s1 * (
        spherical_hankel1(0, k * distance)
        + spherical_hankel1(2, k * distance)
    )
    observed = dipolar_coupling_entry(distance, k, radius, f1)
    np.testing.assert_allclose(observed, expected, rtol=3e-14, atol=3e-16)


def test_analytic_matrix_is_exactly_the_balanced_solver_operator():
    positions = irregular_cluster(6, 2.5)
    solution = solve_multipolar_nodal(
        positions, k=0.1, radius=1.0, f0=0.0, f1=0.8, lmax=1
    )
    production = np.eye(len(solution.system_matrix)) - solution.system_matrix
    analytical = dipolar_balanced_coupling_matrix(
        positions, k=0.1, radius=1.0, f1=0.8
    )
    np.testing.assert_allclose(analytical, production, rtol=4e-14, atol=4e-16)


def test_pair_spectral_radius_has_closed_retarded_form():
    distance, k, radius, f1 = 2.1, 0.1, 1.0, 0.8
    positions = np.array([[-distance / 2, 0, 0], [distance / 2, 0, 0]])
    diagnostics = dipolar_coupling_diagnostics(positions, k, radius, f1)
    expected = 0.5 * abs(f1) * (radius / distance) ** 3
    expected *= np.sqrt(1.0 + (k * distance) ** 2)
    np.testing.assert_allclose(diagnostics.spectral_radius, expected, rtol=3e-15)
    np.testing.assert_allclose(diagnostics.infinity_norm, expected, rtol=3e-15)


def test_near_field_entry_has_quadratic_first_retardation_correction():
    distance, radius, f1 = 2.1, 1.0, 0.8
    k = 1.0e-3 / distance
    exact = dipolar_coupling_entry(distance, k, radius, f1)
    quasistatic = near_field_dipolar_coupling_entry(
        distance, k, radius, f1
    )
    phase = k * distance
    correction = exact / quasistatic - 1.0
    np.testing.assert_allclose(
        correction.real / phase**2, 0.5, rtol=6e-7, atol=0.0
    )
    np.testing.assert_allclose(
        correction.imag / phase**3, 1.0 / 3.0, rtol=2e-6, atol=0.0
    )


def test_row_sum_bounds_spectral_radius_and_sample_is_nearly_normal():
    diagnostics = dipolar_coupling_diagnostics(
        irregular_cluster(10, 2.1), 0.1, 1.0, 1.0
    )
    assert diagnostics.spectral_radius <= diagnostics.infinity_norm
    assert diagnostics.spectral_radius <= diagnostics.spectral_norm
    assert diagnostics.normalized_commutator >= 0.0

    nearly_normal = dipolar_coupling_diagnostics(
        compact_cluster(10, 2.1), 0.1, 1.0, 1.0
    )
    assert nearly_normal.spectral_norm / nearly_normal.spectral_radius < 1.001


def test_neumann_partial_solutions_converge_to_direct_solve():
    coupling = dipolar_balanced_coupling_matrix(
        compact_cluster(10, 2.1), 0.1, 1.0, 1.0
    )
    source = np.ones(len(coupling), dtype=complex)
    partials = neumann_partial_solutions(coupling, source, 12)
    exact = np.linalg.solve(np.eye(len(coupling)) - coupling, source)
    errors = np.linalg.norm(partials - exact, axis=1) / np.linalg.norm(exact)
    assert np.all(errors[1:] < errors[:-1])
    assert errors[-1] < 2e-8


def test_zero_order_neumann_force_is_model_a():
    positions = irregular_cluster(4, 2.5)
    k, radius, energy, f1 = 0.1, 1.0, 1.0, 0.8
    comparison = compare_nodal_force_models(
        positions, k, radius, energy, 0.0, f1
    )
    s1 = rayleigh_scattering_coefficients(k * radius, 0.0, f1)[1]
    coefficients = np.zeros((len(positions), 4), dtype=complex)
    coefficients[:, mode_index(1, 0)] = (
        s1 * nodal_standing_wave_coefficients(1)[mode_index(1, 0)]
    )
    local = np.zeros((len(positions), 9), dtype=complex)
    for target in range(len(positions)):
        for source in range(len(positions)):
            if target != source:
                local[target] += translation_matrix(
                    k, positions[target], positions[source], 2, 1
                ) @ coefficients[source]
    b_minus = local[:, mode_index(2, -1)]
    b_plus = local[:, mode_index(2, 1)]
    prefactor = np.sqrt(30.0 * np.pi) * k * radius**3 * energy / 15.0
    observed = np.column_stack(
        (
            prefactor * np.real(f1 * (b_minus - b_plus)),
            prefactor * np.real(-1j * f1 * (b_plus + b_minus)),
        )
    )
    np.testing.assert_allclose(
        observed, comparison.model_a_forces_xy, rtol=4e-13, atol=4e-14
    )


@pytest.mark.parametrize(
    "arguments",
    [
        (1.9, 0.1, 1.0, 0.8),
        (2.1, 0.0, 1.0, 0.8),
        (2.1, 0.1, -1.0, 0.8),
        (2.1, 0.2, 1.0, 0.8),
        (2.1, 0.1, 1.0, 1.1),
    ],
)
def test_invalid_entry_parameters_are_rejected(arguments):
    with pytest.raises(ValueError):
        dipolar_coupling_entry(*arguments)


def test_invalid_neumann_inputs_are_rejected():
    with pytest.raises(ValueError):
        neumann_partial_solutions(np.ones((2, 3)), np.ones(2), 2)
    with pytest.raises(ValueError):
        neumann_partial_solutions(np.eye(2), np.ones(3), 2)
    with pytest.raises(ValueError):
        neumann_partial_solutions(np.eye(2), np.ones(2), -1)
