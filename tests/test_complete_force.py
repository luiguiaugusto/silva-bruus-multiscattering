"""Tests of the complete multipolar radiation-force functional."""

import numpy as np
import pytest

from acoustic_ms import (
    complete_radiation_force_from_bsc,
    nodal_standing_wave_coefficients,
    rayleigh_multipolar_scattering_coefficients,
    solve_model_e_nodal,
    solve_multipolar_nodal_interaction_forces,
    translation_matrix,
)
from scripts.t11_stress_oracle import stress_tensor_force


def test_complete_force_is_real_finite_and_quadratic():
    rng = np.random.default_rng(413)
    bsc = rng.normal(size=16) + 1j * rng.normal(size=16)
    scattering = np.array([0.02j, 0.01j, 0.003j, 0.0002j])
    force = complete_radiation_force_from_bsc(bsc, scattering, 0.4, 1.7)
    assert force.shape == (3,)
    assert np.all(np.isfinite(force))
    np.testing.assert_allclose(
        complete_radiation_force_from_bsc(2.0 * bsc, scattering, 0.4, 1.7),
        4.0 * force,
        rtol=2e-15,
        atol=2e-15,
    )
    phase = np.exp(0.73j)
    np.testing.assert_allclose(
        complete_radiation_force_from_bsc(phase * bsc, scattering, 0.4, 1.7),
        force,
        rtol=2e-15,
        atol=2e-15,
    )
    np.testing.assert_allclose(
        complete_radiation_force_from_bsc(bsc, scattering, 0.4, 3.4),
        2.0 * force,
        rtol=2e-15,
        atol=2e-15,
    )


def test_complete_force_zero_inputs_and_nodal_external_force():
    scattering = np.array([0.01j, 0.002j, 0.0001j])
    assert np.all(
        complete_radiation_force_from_bsc(np.zeros(9), scattering, 0.1, 1.0)
        == 0.0
    )
    assert np.all(
        complete_radiation_force_from_bsc(
            nodal_standing_wave_coefficients(2), scattering, 0.1, 1.0
        )
        == 0.0
    )
    assert np.all(
        complete_radiation_force_from_bsc(
            nodal_standing_wave_coefficients(2), scattering, 0.1, 0.0
        )
        == 0.0
    )


def test_complete_force_decomposition_is_exact_for_arbitrary_bscs():
    rng = np.random.default_rng(927)
    a = rng.normal(size=25) + 1j * rng.normal(size=25)
    c = rng.normal(size=25) + 1j * rng.normal(size=25)
    scattering = 1j * np.linspace(0.02, 0.001, 5)
    functional = lambda value: complete_radiation_force_from_bsc(
        value, scattering, 0.2, 1.0
    )
    interaction = functional(a + c) - functional(a)
    external_scattered = interaction - functional(c)
    np.testing.assert_allclose(
        functional(a) + external_scattered + functional(c),
        functional(a + c),
        rtol=2e-15,
        atol=2e-15,
    )


def test_complete_force_reduces_to_model_d_external_scattered_force():
    positions = np.array([[-1.05, 0.0, 0.0], [1.05, 0.0, 0.0]])
    model_d = solve_multipolar_nodal_interaction_forces(
        positions, 0.1, 1.0, 1.0, 0.0, 0.8, 1
    )
    external = np.tile(nodal_standing_wave_coefficients(2), (2, 1))
    translated = np.zeros_like(external)
    for target in range(2):
        translated[target] = translation_matrix(
            0.1, positions[target], positions[1 - target], 2, 1
        ) @ model_d.solution.coefficients[1 - target]
    scattering = rayleigh_multipolar_scattering_coefficients(0.1, 0.0, 0.8, 2)
    scattering[[0, 2]] = 0.0
    cross = np.array([
        complete_radiation_force_from_bsc(external[i] + translated[i], scattering, 0.1, 1.0)
        - complete_radiation_force_from_bsc(external[i], scattering, 0.1, 1.0)
        - complete_radiation_force_from_bsc(translated[i], scattering, 0.1, 1.0)
        for i in range(2)
    ])
    np.testing.assert_allclose(cross[:, :2], model_d.forces_xy, rtol=3e-12, atol=3e-13)


@pytest.mark.parametrize("radius", [1.01, 1.04])
@pytest.mark.parametrize("resolution", [(24, 48), (32, 64)])
def test_complete_force_matches_independent_stress_tensor(radius, resolution):
    positions = np.array([[-1.05, 0.0, 0.0], [1.05, 0.0, 0.0]])
    result = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.8, 4)
    theta_order, phi_count = resolution
    oracle = stress_tensor_force(
        result.solution.effective_incident_coefficients[0],
        result.solution.scattering_coefficients,
        0.1,
        1.0,
        radius,
        theta_order,
        phi_count,
    )
    error = np.linalg.norm(oracle - result.total_forces_xyz[0]) / np.linalg.norm(
        result.total_forces_xyz[0]
    )
    assert error < 1e-10


@pytest.mark.parametrize(
    "bsc,scattering,k,energy",
    [
        (np.zeros(8), np.zeros(3), 0.1, 1.0),
        (np.zeros(9), np.zeros((3, 1)), 0.1, 1.0),
        (np.zeros(9), np.zeros(3), 0.0, 1.0),
        (np.zeros(9), np.zeros(3), 0.1, -1.0),
        (np.full(9, np.nan), np.zeros(3), 0.1, 1.0),
    ],
)
def test_complete_force_rejects_invalid_inputs(bsc, scattering, k, energy):
    with pytest.raises((TypeError, ValueError)):
        complete_radiation_force_from_bsc(bsc, scattering, k, energy)
