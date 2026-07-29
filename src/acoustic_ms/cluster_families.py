"""Deterministic planar cluster families used by the T08 transferability study."""

from dataclasses import dataclass

import numpy as np

from .geometries import equilateral_trimer, irregular_quartet, scalene_trimer, square_quartet


@dataclass(frozen=True)
class TransferabilityConfiguration:
    """One physical T08 configuration before force evaluation."""

    case_id: str
    split: str
    particle_count: int
    family: str
    distance_ratio: float
    f1: float
    positions_xyz: np.ndarray


_IRREGULAR = {
    6: np.array([
        [0.00, 0.00], [1.05, 0.08], [2.18, -0.06],
        [0.32, 1.14], [1.49, 1.03], [0.91, 2.22],
    ]),
    10: np.array([
        [0.00, 0.00], [1.10, 0.05], [2.25, -0.10], [3.30, 0.18],
        [0.35, 1.12], [1.48, 0.95], [2.62, 1.30],
        [0.82, 2.15], [1.92, 2.30], [1.28, 3.28],
    ]),
}


def _validated_distance(distance: object) -> float:
    try:
        value = float(distance)
    except (TypeError, ValueError) as exc:
        raise ValueError("distance must be a real scalar") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError("distance must be finite and positive")
    return value


def _normalize_template(points_xy: np.ndarray, distance: float) -> np.ndarray:
    separations = points_xy[:, None, :] - points_xy[None, :, :]
    norms = np.linalg.norm(separations, axis=2)
    minimum = float(np.min(norms[np.triu_indices(len(points_xy), 1)]))
    centered = points_xy / minimum
    centered -= centered.mean(axis=0)
    return np.column_stack((distance * centered, np.zeros(len(points_xy))))


def linear_cluster(particle_count: int, distance: object) -> np.ndarray:
    """Return a centered planar chain for ``particle_count >= 2``."""
    if not isinstance(particle_count, int) or isinstance(particle_count, bool) or particle_count < 2:
        raise ValueError("particle_count must be an integer >= 2")
    distance = _validated_distance(distance)
    x = (np.arange(particle_count) - (particle_count - 1) / 2) * distance
    return np.column_stack((x, np.zeros(particle_count), np.zeros(particle_count)))


def compact_cluster(particle_count: int, distance: object) -> np.ndarray:
    """Return the canonical compact family for N=3,4,6,10."""
    distance = _validated_distance(distance)
    if particle_count == 3:
        return equilateral_trimer(distance)
    if particle_count == 4:
        return square_quartet(distance)
    if particle_count not in (6, 10):
        raise ValueError("compact clusters are defined for N in {3, 4, 6, 10}")
    rows = 3 if particle_count == 6 else 4
    points = []
    for row in range(rows):
        count = rows - row
        for column in range(count):
            points.append((column + row / 2, row * np.sqrt(3) / 2))
    return _normalize_template(np.asarray(points), distance)


def irregular_cluster(particle_count: int, distance: object) -> np.ndarray:
    """Return the canonical irregular family for N=3,4,6,10."""
    distance = _validated_distance(distance)
    if particle_count == 3:
        return scalene_trimer(distance)
    if particle_count == 4:
        return irregular_quartet(distance)
    if particle_count not in _IRREGULAR:
        raise ValueError("irregular clusters are defined for N in {3, 4, 6, 10}")
    return _normalize_template(_IRREGULAR[particle_count], distance)


def cluster_family(particle_count: int, family: str, distance: object) -> np.ndarray:
    """Build one of the thirteen geometries prescribed by T08."""
    if particle_count == 2 and family == "pair":
        return linear_cluster(2, distance)
    if particle_count not in (3, 4, 6, 10):
        raise ValueError("particle_count and family do not identify a T08 geometry")
    builders = {"linear": linear_cluster, "compact": compact_cluster, "irregular": irregular_cluster}
    if family not in builders:
        raise ValueError("family must be pair, linear, compact, or irregular")
    return builders[family](particle_count, distance)


def enumerate_transferability_configurations() -> tuple[TransferabilityConfiguration, ...]:
    """Return the exact deterministic 312-case T08 enumeration."""
    geometries = [(2, "pair")]
    geometries.extend((n, family) for n in (3, 4, 6, 10) for family in ("linear", "compact", "irregular"))
    configurations = []
    for particle_count, family in geometries:
        split = "calibration" if particle_count <= 4 else "holdout"
        for f1 in (0.1, 0.4, 0.8, 1.0):
            for distance in (2.1, 2.5, 3.0, 4.0, 6.0, 10.0):
                case_id = f"n{particle_count}_{family}_f{f1:.1f}_d{distance:.1f}"
                configurations.append(TransferabilityConfiguration(
                    case_id, split, particle_count, family, distance, f1,
                    cluster_family(particle_count, family, distance),
                ))
    return tuple(configurations)
