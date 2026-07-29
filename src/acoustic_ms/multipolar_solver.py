"""Balanced global multipolar solver for planar nodal clusters."""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .incident import nodal_standing_wave_coefficients
from .multipolar_scattering import rayleigh_multipolar_scattering_coefficients
from .multipoles import mode_count, modes
from .solver import _validate_positions
from .translation import separation_coefficient


@dataclass(frozen=True)
class MultipolarNodalSolution:
    """Full coefficients and balanced/physical diagnostics for Model D."""

    coefficients: np.ndarray
    system_matrix: np.ndarray
    right_hand_side: np.ndarray
    physical_system_matrix: np.ndarray
    physical_right_hand_side: np.ndarray
    residual_relative: float
    condition_number: float
    physical_condition_number: float
    modes: tuple[tuple[int, int], ...]
    active_modes: tuple[tuple[int, int], ...]
    active_mode_indices: tuple[int, ...]
    lmax: int
    used_planar_symmetry: bool


def _validated_scalars(
    k: object, radius: object, f0: object, f1: object
) -> tuple[float, float, float, float]:
    try:
        k, radius, f0, f1 = (float(value) for value in (k, radius, f0, f1))
    except (TypeError, ValueError) as exc:
        raise ValueError("k, radius, f0, and f1 must be real scalars") from exc
    if not np.all(np.isfinite((k, radius, f0, f1))):
        raise ValueError("k, radius, f0, and f1 must be finite")
    if k <= 0.0 or radius <= 0.0:
        raise ValueError("k and radius must be positive")
    if k * radius > 0.1:
        raise ValueError("the Rayleigh regime requires ka <= 0.1")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")
    return k, radius, f0, f1


@lru_cache(maxsize=None)
def _cached_translation_block(
    k: float,
    displacement: tuple[float, float, float],
    selected_modes: tuple[tuple[int, int], ...],
) -> np.ndarray:
    """Return one translation block for a displacement and selected modes."""
    origin = np.zeros(3)
    source = np.asarray(displacement)
    block = np.empty((len(selected_modes), len(selected_modes)), dtype=complex)
    for row, (target_ell, target_m) in enumerate(selected_modes):
        for column, (source_ell, source_m) in enumerate(selected_modes):
            block[row, column] = separation_coefficient(
                target_ell,
                target_m,
                source_ell,
                source_m,
                k,
                origin,
                source,
            )
    block.setflags(write=False)
    return block


def _selected_modes(
    full_modes: tuple[tuple[int, int], ...],
    scattering_by_ell: np.ndarray,
    use_planar_symmetry: bool,
) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]]:
    selected = []
    indices = []
    for index, (ell, m) in enumerate(full_modes):
        if scattering_by_ell[ell] == 0.0:
            continue
        if use_planar_symmetry and (ell + m) % 2 == 0:
            continue
        selected.append((ell, m))
        indices.append(index)
    return tuple(selected), tuple(indices)


def solve_multipolar_nodal(
    positions_xyz: object,
    k: object,
    radius: object,
    f0: object,
    f1: object,
    lmax: int = 1,
    *,
    use_planar_symmetry: bool = True,
) -> MultipolarNodalSolution:
    """Solve the globally coupled nodal Model D system in a balanced basis.

    Planar reflection keeps modes with odd ``ell + m``.  The returned
    coefficient matrix always uses the complete project ordering; omitted
    modes are represented by exact zeros.
    """
    if not isinstance(lmax, int) or isinstance(lmax, bool) or lmax < 1:
        raise ValueError("lmax must be a positive integer")
    if not isinstance(use_planar_symmetry, bool):
        raise ValueError("use_planar_symmetry must be boolean")
    k, radius, f0, f1 = _validated_scalars(k, radius, f0, f1)
    positions = _validate_positions(positions_xyz, radius)
    particle_count = len(positions)
    full_modes = modes(lmax)
    full_count = mode_count(lmax)
    scattering_by_ell = rayleigh_multipolar_scattering_coefficients(
        k * radius, f0, f1, lmax
    )
    active_modes, active_indices = _selected_modes(
        full_modes, scattering_by_ell, use_planar_symmetry
    )
    active_count = len(active_modes)
    full_coefficients = np.zeros((particle_count, full_count), dtype=complex)

    if active_count == 0:
        empty_matrix = np.empty((0, 0), dtype=complex)
        empty_vector = np.empty(0, dtype=complex)
        return MultipolarNodalSolution(
            coefficients=full_coefficients,
            system_matrix=empty_matrix,
            right_hand_side=empty_vector,
            physical_system_matrix=empty_matrix.copy(),
            physical_right_hand_side=empty_vector.copy(),
            residual_relative=0.0,
            condition_number=1.0,
            physical_condition_number=1.0,
            modes=full_modes,
            active_modes=active_modes,
            active_mode_indices=active_indices,
            lmax=lmax,
            used_planar_symmetry=use_planar_symmetry,
        )

    local_scattering = np.array(
        [scattering_by_ell[ell] for ell, _ in active_modes], dtype=complex
    )
    global_scattering = np.tile(local_scattering, particle_count)
    translation = np.zeros(
        (particle_count * active_count, particle_count * active_count),
        dtype=complex,
    )
    for target in range(particle_count):
        for source in range(particle_count):
            if target == source:
                continue
            rows = slice(target * active_count, (target + 1) * active_count)
            columns = slice(source * active_count, (source + 1) * active_count)
            displacement = tuple(float(value) for value in positions[source] - positions[target])
            translation[rows, columns] = _cached_translation_block(
                k, displacement, active_modes
            )

    external_full = nodal_standing_wave_coefficients(lmax)
    external_active = external_full[np.asarray(active_indices)]
    external_global = np.tile(external_active, particle_count)
    physical_system = (
        np.eye(particle_count * active_count, dtype=complex)
        - global_scattering[:, None] * translation
    )
    physical_rhs = global_scattering * external_global

    square_root = np.sqrt(global_scattering)
    balanced_system = (
        np.eye(particle_count * active_count, dtype=complex)
        - square_root[:, None] * translation * square_root[None, :]
    )
    balanced_rhs = square_root * external_global
    balanced_solution = np.linalg.solve(balanced_system, balanced_rhs)
    physical_solution = square_root * balanced_solution

    residual_denominator = max(
        float(np.linalg.norm(physical_rhs)), np.finfo(float).eps
    )
    residual = float(
        np.linalg.norm(physical_system @ physical_solution - physical_rhs)
        / residual_denominator
    )
    active_coefficients = physical_solution.reshape(particle_count, active_count)
    full_coefficients[:, np.asarray(active_indices)] = active_coefficients
    return MultipolarNodalSolution(
        coefficients=full_coefficients,
        system_matrix=balanced_system,
        right_hand_side=balanced_rhs,
        physical_system_matrix=physical_system,
        physical_right_hand_side=physical_rhs,
        residual_relative=residual,
        condition_number=float(np.linalg.cond(balanced_system)),
        physical_condition_number=float(np.linalg.cond(physical_system)),
        modes=full_modes,
        active_modes=active_modes,
        active_mode_indices=active_indices,
        lmax=lmax,
        used_planar_symmetry=use_planar_symmetry,
    )
