import numpy as np
import pytest

from acoustic_ms.silva_bruus import (
    nodal_pair_force_magnitude,
    nodal_pair_force_on_probe,
    nodal_pair_forces,
)


K = 2.0
RADIUS = 0.1
ENERGY = 3.0
F1 = 0.5


def potential(distance: float, k: float = K, radius: float = RADIUS, energy: float = ENERGY, f1: float = F1) -> float:
    """Independent implementation of the specified interaction potential."""
    x = k * distance
    n1 = -np.sin(x) / x - np.cos(x) / x**2
    return 2.0 * np.pi * energy * k**3 * radius**6 * f1**2 * n1 / x


@pytest.mark.parametrize("x", [0.7, 1.3, 4.7, 9.2])
def test_closed_force_matches_negative_numerical_potential_derivative(x):
    distance = x / K
    step = distance * 1.0e-6
    expected = -(potential(distance + step) - potential(distance - step)) / (2.0 * step)
    actual = nodal_pair_force_magnitude(K, RADIUS, distance, ENERGY, F1)
    assert np.isclose(actual, expected, rtol=2.0e-9, atol=1.0e-13)


def test_near_field_limit_and_attraction():
    k = 1.0e-5
    distance = 1.0
    actual = nodal_pair_force_magnitude(k, RADIUS, distance, ENERGY, F1)
    expected = -6.0 * np.pi * RADIUS**2 * ENERGY * F1**2 * (RADIUS / distance) ** 4
    assert actual < 0.0
    assert np.isclose(actual / expected, 1.0, rtol=1.0e-8)


def test_far_field_limit_away_from_cosine_zero():
    k = 1.0
    radius = 1.0e-4
    distance = 100.0 * np.pi
    actual = nodal_pair_force_magnitude(k, radius, distance, ENERGY, F1)
    expected = 2.0 * np.pi * radius**2 * ENERGY * F1**2 * (k * radius) ** 2 * (radius / distance) ** 2 * np.cos(k * distance)
    assert np.isclose(actual / expected, 1.0, rtol=3.0e-4)


def test_action_reaction_translation_and_rotation_covariance():
    p1 = np.array([0.3, -0.4])
    p2 = np.array([-0.5, 0.1])
    force_1, force_2 = nodal_pair_forces(p1, p2, K, RADIUS, ENERGY, F1)
    assert np.allclose(force_1, -force_2)

    shift = np.array([7.0, -2.0])
    shifted = nodal_pair_force_on_probe(p1 + shift, p2 + shift, K, RADIUS, ENERGY, F1)
    assert np.allclose(shifted, force_1)

    theta = 0.71
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    rotated = nodal_pair_force_on_probe(rotation @ p1, rotation @ p2, K, RADIUS, ENERGY, F1)
    assert np.allclose(rotated, rotation @ force_1)


def test_force_is_zero_for_zero_contrast_or_energy():
    p1 = np.array([0.0, 0.0])
    p2 = np.array([1.0, 0.0])
    assert np.allclose(nodal_pair_force_on_probe(p1, p2, K, RADIUS, ENERGY, 0.0), 0.0)
    assert np.allclose(nodal_pair_force_on_probe(p1, p2, K, RADIUS, 0.0, F1), 0.0)


def test_overlap_coincident_and_invalid_inputs_are_rejected():
    with pytest.raises(ValueError, match="distance"):
        nodal_pair_force_magnitude(K, RADIUS, 0.19, ENERGY, F1)
    with pytest.raises(ValueError, match="coincide"):
        nodal_pair_force_on_probe([0.0, 0.0], [0.0, 0.0], K, RADIUS, ENERGY, F1)

    for name, value in [("k", 0.0), ("radius", 0.0), ("energy_density", -1.0), ("f1", 1.1)]:
        arguments = dict(k=K, radius=RADIUS, distance=1.0, energy_density=ENERGY, f1=F1)
        arguments[name] = value
        with pytest.raises(ValueError):
            nodal_pair_force_magnitude(**arguments)
    with pytest.raises(ValueError):
        nodal_pair_force_on_probe([np.nan, 0.0], [1.0, 0.0], K, RADIUS, ENERGY, F1)
