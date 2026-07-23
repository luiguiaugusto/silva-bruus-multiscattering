import numpy as np
import pytest

from acoustic_ms import (
    corrected_nodal_pair_force_magnitude,
    corrected_nodal_pair_force_on_probe,
    corrected_nodal_pair_forces,
    corrected_pair_coefficients,
    nodal_pair_force_magnitude,
)


KA = 0.1
RADIUS = 1.0
K = KA / RADIUS
ENERGY = 1.0


def relative_error(kd: float, f1: float) -> float:
    distance = kd / K
    corrected = corrected_nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, f1)
    reference = nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, f1)
    return 100.0 * abs(corrected - reference) / abs(corrected)


@pytest.mark.parametrize(
    ("f1", "a0_over_d0", "a2_over_d0"),
    [
        (0.1, -0.506355807146221, 2.537084546937989),
        (0.4, -0.527267939835848, 2.659985720443811),
        (0.8, -0.560671850474030, 2.858882625216557),
        (1.0, -0.580568297988085, 2.978697291068864),
    ],
)
def test_coefficient_ratios_match_independent_contact_values(f1, a0_over_d0, a2_over_d0):
    a0, a2, d0 = corrected_pair_coefficients(f1, 0.5)
    assert np.isclose(a0 / d0, a0_over_d0, rtol=2e-12)
    assert np.isclose(a2 / d0, a2_over_d0, rtol=2e-12)


@pytest.mark.parametrize(
    ("f1", "expected"),
    [(0.1, 1.252519728707), (0.4, 5.160511340274), (0.8, 10.798343941865), (1.0, 13.848266387733)],
)
def test_contact_relative_errors_match_published_regression_values(f1, expected):
    assert np.isclose(relative_error(0.2, f1), expected, rtol=2e-12)


def test_additional_published_relative_error_value():
    assert np.isclose(relative_error(0.3, 1.0), 2.206083641258, rtol=2e-12)


def test_zero_contrast_coefficient_limits_are_nontrivial():
    a0, a2, d0 = corrected_pair_coefficients(0.0, 0.5)
    assert a0 / d0 == -0.5
    assert a2 / d0 == 2.5


def test_corrected_force_converges_to_silva_bruus_for_small_nonzero_contrast():
    distance = 5.0
    corrected = corrected_nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, 1e-7)
    reference = nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, 1e-7)
    assert np.isclose(corrected, reference, rtol=1e-7, atol=1e-25)


def test_corrected_force_converges_to_silva_bruus_at_large_separation():
    distance = 1.0e5
    corrected = corrected_nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, 0.8)
    reference = nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, 0.8)
    assert np.isclose(corrected, reference, rtol=1e-12)


def test_small_kd_limit_and_near_contact_attraction():
    f1 = 0.8
    distance = 2.5
    k = 1e-6
    a0, _, d0 = corrected_pair_coefficients(f1, RADIUS / distance)
    expected = 12.0 * np.pi * RADIUS**2 * f1**2 * ENERGY * (RADIUS / distance) ** 4 * a0 / d0
    actual = corrected_nodal_pair_force_magnitude(k, RADIUS, distance, ENERGY, f1)
    assert actual < 0.0
    assert np.isclose(actual, expected, rtol=1e-11)


def test_vector_symmetries_action_reaction_and_zero_forces():
    p1 = np.array([2.0, -1.0])
    p2 = np.array([-1.0, 0.5])
    force_1, force_2 = corrected_nodal_pair_forces(p1, p2, K, RADIUS, ENERGY, 0.8)
    assert np.allclose(force_1, -force_2)
    shift = np.array([3.0, 9.0])
    assert np.allclose(corrected_nodal_pair_force_on_probe(p1 + shift, p2 + shift, K, RADIUS, ENERGY, 0.8), force_1)
    angle = 0.4
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    assert np.allclose(corrected_nodal_pair_force_on_probe(rotation @ p1, rotation @ p2, K, RADIUS, ENERGY, 0.8), rotation @ force_1)
    assert np.allclose(corrected_nodal_pair_force_on_probe(p1, p2, K, RADIUS, 0.0, 0.8), 0.0)
    assert np.allclose(corrected_nodal_pair_force_on_probe(p1, p2, K, RADIUS, ENERGY, 0.0), 0.0)


def test_invalid_domains_are_rejected():
    with pytest.raises(ValueError):
        corrected_pair_coefficients(1.1, 0.5)
    with pytest.raises(ValueError):
        corrected_pair_coefficients(0.1, 0.0)
    with pytest.raises(ValueError):
        corrected_pair_coefficients(np.nan, 0.5)
    with pytest.raises(ValueError):
        corrected_nodal_pair_force_magnitude(K, RADIUS, 1.9, ENERGY, 0.5)
    with pytest.raises(ValueError):
        corrected_nodal_pair_force_on_probe([0.0, 0.0], [0.0, 0.0], K, RADIUS, ENERGY, 0.5)
    with pytest.raises(ValueError):
        corrected_nodal_pair_force_magnitude(np.inf, RADIUS, 2.0, ENERGY, 0.5)


def test_figure_error_decreases_monotonically_from_contact_to_kd_point_three():
    kd_values = np.linspace(0.2, 0.3, 101)
    for f1 in (0.1, 0.4, 0.8, 1.0):
        errors = np.array([relative_error(kd, f1) for kd in kd_values])
        assert np.all(np.diff(errors) < 0.0)
