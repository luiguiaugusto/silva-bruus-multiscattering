"""Corrected fifth-order two-sphere nodal interaction force.

This module implements Eqs. (30a)--(30d) of the 2026 two-particle benchmark.
It is an analytical truncated result, not a general multiple-scattering solver.
"""

import numpy as np

from .silva_bruus import _finite_scalar, _validate_parameters, _xy_position


def corrected_pair_coefficients(f1: float, separation_ratio: float) -> tuple[float, float, float]:
    """Return the ``(A0, A2, D0)`` coefficients for ``chi = a / d``.

    ``f1`` and ``separation_ratio`` are dimensionless; the latter must satisfy
    ``0 < chi <= 1/2`` for non-overlapping identical spheres.
    """
    f1 = _finite_scalar(f1, "f1")
    chi = _finite_scalar(separation_ratio, "separation_ratio")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")
    if not 0.0 < chi <= 0.5:
        raise ValueError("separation_ratio must satisfy 0 < separation_ratio <= 0.5")

    a0 = (
        4.0 * (77.0 - 25.0 * f1 + 2.0 * f1**2)
        + 90.0 * f1 * (-11.0 + 2.0 * f1) * chi**7
        + 5250.0 * (-7.0 + f1) * f1 * chi**11
        + 33075.0 * f1**2 * chi**18
    )
    a2 = (
        -20.0 * (77.0 - 25.0 * f1 + 2.0 * f1**2)
        - 162.0 * f1 * (-11.0 + 2.0 * f1) * chi**7
        - 24250.0 * (-7.0 + f1) * f1 * chi**11
        + 20025.0 * f1**2 * chi**18
    )
    d0 = (
        -616.0
        + 4.0 * f1 * (50.0 + 77.0 * chi**3 + 1485.0 * chi**7 + 22050.0 * chi**11)
        + f1**3 * chi**3 * (8.0 + 324.0 * chi**7 + 5850.0 * chi**11 + 73575.0 * chi**18)
        - 2.0 * f1**2 * (
            8.0
            + 50.0 * chi**3
            + 9.0 * chi**7 * (60.0 + chi**3 * (99.0 + 700.0 * chi + 2275.0 * chi**4 + 17850.0 * chi**8))
        )
    )
    return float(a0), float(a2), float(d0)


def corrected_nodal_pair_force_magnitude(
    k: float, radius: float, distance: float, energy_density: float, f1: float
) -> float:
    """Return the signed corrected radial pair force in newtons.

    The value is a radial component along the source-to-probe direction, not
    an unsigned magnitude; negative values are attractive.
    """
    k, radius, energy_density, f1 = _validate_parameters(k, radius, energy_density, f1)
    distance = _finite_scalar(distance, "distance")
    if distance < 2.0 * radius:
        raise ValueError("distance must satisfy distance >= 2 * radius")

    chi = radius / distance
    a0, a2, d0 = corrected_pair_coefficients(f1, chi)
    x = k * distance
    bracket = 3.0 * a0 / d0 * (np.cos(x) + x * np.sin(x)) + a2 / (5.0 * d0) * x**2 * np.cos(x)
    return float(4.0 * np.pi * radius**2 * f1**2 * energy_density * chi**4 * bracket)


def corrected_nodal_pair_force_on_probe(
    probe_xy: object, source_xy: object, k: float, radius: float,
    energy_density: float, f1: float,
) -> np.ndarray:
    """Return the corrected 2D force (N) on a probe due to a source."""
    probe = _xy_position(probe_xy, "probe_xy")
    source = _xy_position(source_xy, "source_xy")
    displacement = probe - source
    distance = float(np.linalg.norm(displacement))
    if distance == 0.0:
        raise ValueError("probe_xy and source_xy must not coincide")
    radial_force = corrected_nodal_pair_force_magnitude(k, radius, distance, energy_density, f1)
    return radial_force * displacement / distance


def corrected_nodal_pair_forces(
    position_1_xy: object, position_2_xy: object, k: float, radius: float,
    energy_density: float, f1: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return equal-and-opposite corrected forces (N) on two particles."""
    force_1 = corrected_nodal_pair_force_on_probe(
        position_1_xy, position_2_xy, k, radius, energy_density, f1
    )
    return force_1, -force_1
