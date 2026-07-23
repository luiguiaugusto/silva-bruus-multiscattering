"""Silva--Bruus pairwise force for identical spheres in a nodal plane."""

import numpy as np


def _finite_scalar(value: float, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a real finite scalar") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _validate_parameters(k: float, radius: float, energy_density: float, f1: float) -> tuple[float, float, float, float]:
    k = _finite_scalar(k, "k")
    radius = _finite_scalar(radius, "radius")
    energy_density = _finite_scalar(energy_density, "energy_density")
    f1 = _finite_scalar(f1, "f1")
    if k <= 0.0:
        raise ValueError("k must be positive")
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if energy_density < 0.0:
        raise ValueError("energy_density must be non-negative")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")
    return k, radius, energy_density, f1


def _xy_position(position: object, name: str) -> np.ndarray:
    array = np.asarray(position, dtype=float)
    if array.shape != (2,):
        raise ValueError(f"{name} must be a two-component planar position")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def nodal_pair_force_magnitude(
    k: float, radius: float, distance: float, energy_density: float, f1: float
) -> float:
    """Return the signed radial nodal pair force in newtons.

    The radial direction is from source to probe.  Thus a negative result is
    attractive.  ``k`` is in ``m^-1``, ``radius`` and ``distance`` in metres,
    and ``energy_density`` is in ``J m^-3``.
    """
    k, radius, energy_density, f1 = _validate_parameters(k, radius, energy_density, f1)
    distance = _finite_scalar(distance, "distance")
    if distance < 2.0 * radius:
        raise ValueError("distance must satisfy distance >= 2 * radius")

    x = k * distance
    bracket = -1.5 * (np.cos(x) + x * np.sin(x)) + 0.5 * x**2 * np.cos(x)
    return float(4.0 * np.pi * radius**2 * energy_density * f1**2 * (radius / distance) ** 4 * bracket)


def nodal_pair_force_on_probe(
    probe_xy: object, source_xy: object, k: float, radius: float,
    energy_density: float, f1: float,
) -> np.ndarray:
    """Return the 2D force (N) on ``probe_xy`` due to ``source_xy``.

    Both positions are in metres and belong to the pressure nodal plane.
    """
    probe = _xy_position(probe_xy, "probe_xy")
    source = _xy_position(source_xy, "source_xy")
    displacement = probe - source
    distance = float(np.linalg.norm(displacement))
    if distance == 0.0:
        raise ValueError("probe_xy and source_xy must not coincide")
    magnitude = nodal_pair_force_magnitude(k, radius, distance, energy_density, f1)
    return magnitude * displacement / distance


def nodal_pair_forces(
    position_1_xy: object, position_2_xy: object, k: float, radius: float,
    energy_density: float, f1: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return equal-and-opposite forces (N) on particles 1 and 2."""
    force_1 = nodal_pair_force_on_probe(
        position_1_xy, position_2_xy, k, radius, energy_density, f1
    )
    return force_1, -force_1
