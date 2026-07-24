"""Canonical centered planar trimer and quartet geometries."""

import numpy as np


def _distance(value: float) -> float:
    distance = float(value)
    if not np.isfinite(distance) or distance <= 0.0:
        raise ValueError("distance must be finite and positive")
    return distance


def linear_trimer(distance: float) -> np.ndarray:
    """Return a centered three-particle chain with nearest-neighbor distance."""
    distance = _distance(distance)
    return np.array([[-distance, 0.0, 0.0], [0.0, 0.0, 0.0], [distance, 0.0, 0.0]])


def equilateral_trimer(distance: float) -> np.ndarray:
    """Return a centered equilateral planar trimer of side ``distance``."""
    distance = _distance(distance)
    return np.array([
        [-distance / 2, -np.sqrt(3) * distance / 6, 0.0],
        [distance / 2, -np.sqrt(3) * distance / 6, 0.0],
        [0.0, np.sqrt(3) * distance / 3, 0.0],
    ])


def scalene_trimer(distance: float) -> np.ndarray:
    """Return the fixed-shape centered scalene planar trimer with minimum side."""
    distance = _distance(distance)
    raw = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0], [3 * distance / 11, 12 * distance / 11, 0.0]])
    return raw - raw.mean(axis=0)


def linear_quartet(distance: float) -> np.ndarray:
    """Return a centered four-particle chain with nearest-neighbor distance."""
    distance = _distance(distance)
    return np.array([
        [-1.5 * distance, 0.0, 0.0],
        [-0.5 * distance, 0.0, 0.0],
        [0.5 * distance, 0.0, 0.0],
        [1.5 * distance, 0.0, 0.0],
    ])


def square_quartet(distance: float) -> np.ndarray:
    """Return a centered planar square with side ``distance``."""
    distance = _distance(distance)
    half = distance / 2
    return np.array([
        [-half, -half, 0.0],
        [half, -half, 0.0],
        [half, half, 0.0],
        [-half, half, 0.0],
    ])


def irregular_quartet(distance: float) -> np.ndarray:
    """Return the centered fixed-shape irregular quartet with minimum side."""
    distance = _distance(distance)
    raw = np.array([
        [0.0, 0.0, 0.0],
        [distance, 0.0, 0.0],
        [3 * distance / 10, 6 * distance / 5, 0.0],
        [7 * distance / 5, 11 * distance / 10, 0.0],
    ])
    return raw - raw.mean(axis=0)
