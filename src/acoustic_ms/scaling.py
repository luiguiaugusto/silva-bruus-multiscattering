"""Coupling predictors and unweighted log-space power-law fits."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PowerLawFit:
    """Diagnostics for ``y = prefactor * x**exponent`` fitted in log space."""

    point_count: int
    prefactor: float
    exponent: float
    r_squared_log: float
    rmse_log: float
    max_abs_log_residual: float


def _physical_scalars(radius: object, f1: object) -> tuple[float, float]:
    try:
        radius, f1 = float(radius), float(f1)
    except (TypeError, ValueError) as exc:
        raise ValueError("radius and f1 must be real scalars") from exc
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be finite and positive")
    if not np.isfinite(f1) or not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must be finite and satisfy -2 <= f1 <= 1")
    return radius, f1


def coupling_eta(radius: object, minimum_distance: object, f1: object) -> float:
    """Return ``|f1| * (radius / minimum_distance)**3``."""
    radius, f1 = _physical_scalars(radius, f1)
    try:
        minimum_distance = float(minimum_distance)
    except (TypeError, ValueError) as exc:
        raise ValueError("minimum_distance must be a real scalar") from exc
    if not np.isfinite(minimum_distance) or minimum_distance <= 0.0:
        raise ValueError("minimum_distance must be finite and positive")
    if minimum_distance < 2.0 * radius:
        raise ValueError("minimum_distance must satisfy separation >= 2 * radius")
    return abs(f1) * (radius / minimum_distance) ** 3


def maximum_geometric_coupling(
    positions_xyz: object, radius: object, f1: object
) -> float:
    """Return the maximum per-particle inverse-cube geometric coupling."""
    radius, f1 = _physical_scalars(radius, f1)
    positions = np.asarray(positions_xyz, dtype=float)
    if (
        positions.ndim != 2
        or positions.shape[1:] != (3,)
        or positions.shape[0] < 2
        or not np.all(np.isfinite(positions))
    ):
        raise ValueError("positions_xyz must be finite with shape (N, 3), N >= 2")
    coupling = np.zeros(len(positions))
    for first in range(len(positions)):
        for second in range(first + 1, len(positions)):
            distance = float(np.linalg.norm(positions[first] - positions[second]))
            if distance == 0.0:
                raise ValueError("particle centers must not coincide")
            if distance < 2.0 * radius:
                raise ValueError("particle centers must satisfy separation >= 2 * radius")
            contribution = abs(f1) * (radius / distance) ** 3
            coupling[first] += contribution
            coupling[second] += contribution
    return float(np.max(coupling))


def fit_power_law(x: object, y: object) -> PowerLawFit:
    """Fit ``ln(y) = ln(C) + p ln(x)`` without weighting or filtering."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be one-dimensional with matching shapes")
    if len(x) < 2:
        raise ValueError("at least two points are required")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must be finite")
    if np.any(x <= 0.0) or np.any(y <= 0.0):
        raise ValueError("x and y must be strictly positive")
    if np.all(x == x[0]) or np.all(y == y[0]):
        raise ValueError("x and y must both vary")

    log_x = np.log(x)
    log_y = np.log(y)
    exponent, log_prefactor = np.polyfit(log_x, log_y, 1)
    predicted = log_prefactor + exponent * log_x
    residual = log_y - predicted
    residual_sum = float(np.sum(residual**2))
    total_sum = float(np.sum((log_y - np.mean(log_y)) ** 2))
    return PowerLawFit(
        point_count=len(x),
        prefactor=float(np.exp(log_prefactor)),
        exponent=float(exponent),
        r_squared_log=float(1.0 - residual_sum / total_sum),
        rmse_log=float(np.sqrt(residual_sum / len(x))),
        max_abs_log_residual=float(np.max(np.abs(residual))),
    )
