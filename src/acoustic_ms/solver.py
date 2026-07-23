"""Dense coupled Rayleigh solver in the full Lmax=1 nodal basis."""

from dataclasses import dataclass

import numpy as np

from .incident import nodal_standing_wave_coefficients
from .multipoles import modes
from .scattering import rayleigh_scattering_coefficients
from .translation import translation_matrix


@dataclass(frozen=True)
class RayleighNodalSolution:
    """Coupled field coefficients and diagnostics; no radiation force is computed."""

    coefficients: np.ndarray
    system_matrix: np.ndarray
    right_hand_side: np.ndarray
    residual_relative: float
    condition_number: float
    modes: tuple[tuple[int, int], ...]


def _validate_positions(positions_xyz: object, radius: float) -> np.ndarray:
    positions = np.asarray(positions_xyz, dtype=float)
    if positions.ndim != 2 or positions.shape[1:] != (3,) or positions.shape[0] < 1 or not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite with shape (N, 3), N >= 1")
    scale = max(1.0, radius, float(np.max(np.abs(positions[:, :2]))))
    if np.any(np.abs(positions[:, 2]) > 1e-12 * scale):
        raise ValueError("all centers must lie in the nodal plane z=0 within tolerance")
    for source in range(len(positions)):
        for target in range(source):
            distance = float(np.linalg.norm(positions[source] - positions[target]))
            if distance == 0.0:
                raise ValueError("particle centers must not coincide")
            if distance < 2.0 * radius:
                raise ValueError("particle centers must satisfy separation >= 2 * radius")
    return positions


def solve_rayleigh_nodal(positions_xyz: object, k: float, radius: float, f0: float, f1: float, lmax: int = 1) -> RayleighNodalSolution:
    """Solve ``(I - D_g U)s = D_g a_ext`` for fixed nodal-plane spheres.

    ``Lmax=1`` means four field modes per particle, while the dense solve
    resums arbitrarily long rescattering chains admitted by that basis.
    """
    if lmax != 1:
        raise ValueError("the T03 Rayleigh solver supports lmax=1 only")
    try:
        k, radius, f0, f1 = (float(value) for value in (k, radius, f0, f1))
    except (TypeError, ValueError) as exc:
        raise ValueError("k, radius, f0, and f1 must be real scalars") from exc
    if not all(np.isfinite((k, radius, f0, f1))) or k <= 0.0 or radius <= 0.0 or not -2.0 <= f1 <= 1.0:
        raise ValueError("require finite k > 0, radius > 0, finite f0, and -2 <= f1 <= 1")
    ka = k * radius
    if ka > 0.1:
        raise ValueError("Rayleigh solver requires ka <= 0.1")
    positions = _validate_positions(positions_xyz, radius)
    count = len(positions)
    local_modes = modes(1)
    local_count = len(local_modes)
    coefficients = rayleigh_scattering_coefficients(ka, f0, f1)
    local_diagonal = np.array([coefficients[0], coefficients[1], coefficients[1], coefficients[1]])
    d_global = np.tile(local_diagonal, count)
    translation = np.zeros((local_count * count, local_count * count), dtype=complex)
    for target in range(count):
        for source in range(count):
            if target != source:
                rows = slice(target * local_count, (target + 1) * local_count)
                columns = slice(source * local_count, (source + 1) * local_count)
                translation[rows, columns] = translation_matrix(k, positions[target], positions[source], 1)
    system = np.eye(local_count * count, dtype=complex) - d_global[:, None] * translation
    rhs = d_global * np.tile(nodal_standing_wave_coefficients(1), count)
    solved = np.linalg.solve(system, rhs)
    residual = float(np.linalg.norm(system @ solved - rhs) / max(np.linalg.norm(rhs), np.finfo(float).eps))
    return RayleighNodalSolution(solved.reshape(count, local_count), system, rhs, residual, float(np.linalg.cond(system)), local_modes)
