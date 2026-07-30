"""Globally coupled exact-Mie solver for planar nodal clusters."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numbers

import numpy as np

from .incident import nodal_standing_wave_coefficients
from .mie_scattering import mie_scattering_coefficients_from_contrasts
from .multipoles import mode_count, modes
from .solver import _validate_positions
from .translation import separation_coefficient


@dataclass(frozen=True)
class MieMultiparticleSolution:
    """Effective-incident and scattered coefficients for exact-Mie spheres."""

    effective_incident_coefficients: np.ndarray
    scattered_coefficients: np.ndarray
    external_incident_coefficients: np.ndarray
    scattering_coefficients: np.ndarray
    system_matrix: np.ndarray
    right_hand_side: np.ndarray
    residual_relative: float
    condition_number: float
    modes: tuple[tuple[int, int], ...]
    active_modes: tuple[tuple[int, int], ...]
    active_mode_indices: tuple[int, ...]
    lmax: int
    used_planar_symmetry: bool


def _positive_real(name: str, value: object) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real scalar")
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _positive_order(lmax: object) -> int:
    if isinstance(lmax, (bool, np.bool_)) or not isinstance(lmax, numbers.Integral):
        raise TypeError("lmax must be an integer")
    result = int(lmax)
    if result < 1:
        raise ValueError("lmax must be a positive integer")
    return result


@lru_cache(maxsize=None)
def _translation_block(
    k: float,
    displacement: tuple[float, float, float],
    selected_modes: tuple[tuple[int, int], ...],
) -> np.ndarray:
    origin = np.zeros(3)
    source = np.asarray(displacement)
    block = np.empty((len(selected_modes), len(selected_modes)), dtype=complex)
    for row, (target_ell, target_m) in enumerate(selected_modes):
        for column, (source_ell, source_m) in enumerate(selected_modes):
            block[row, column] = separation_coefficient(
                target_ell, target_m, source_ell, source_m, k, origin, source
            )
    block.setflags(write=False)
    return block


def solve_mie_multiparticle_nodal(
    positions_xyz: object,
    k: object,
    radius: object,
    f0: object,
    f1: object,
    lmax: object,
    *,
    use_planar_symmetry: bool = True,
) -> MieMultiparticleSolution:
    r"""Solve ``(I-U D)b=a_ext`` and return ``d=D b``.

    Exact isolated-sphere Mie coefficients supply the diagonal operator
    ``D``.  Returned coefficient matrices always use the complete project
    ordering ``ell**2 + ell + m``; modes excluded by planar reflection
    symmetry are exact zeros.
    """

    wave_number = _positive_real("k", k)
    sphere_radius = _positive_real("radius", radius)
    order = _positive_order(lmax)
    if not isinstance(use_planar_symmetry, (bool, np.bool_)):
        raise TypeError("use_planar_symmetry must be boolean")
    positions = _validate_positions(positions_xyz, sphere_radius)
    particle_count = len(positions)
    full_modes = modes(order)
    full_count = mode_count(order)
    if use_planar_symmetry:
        active_indices = tuple(
            index for index, (ell, m) in enumerate(full_modes) if (ell + m) % 2 == 1
        )
    else:
        active_indices = tuple(range(full_count))
    active_modes = tuple(full_modes[index] for index in active_indices)
    active_count = len(active_modes)

    scattering_by_ell = mie_scattering_coefficients_from_contrasts(
        wave_number * sphere_radius, f0, f1, order
    )
    local_scattering = np.asarray(
        [scattering_by_ell[ell] for ell, _ in active_modes], dtype=complex
    )
    global_scattering = np.tile(local_scattering, particle_count)
    translation = np.zeros(
        (particle_count * active_count, particle_count * active_count), dtype=complex
    )
    for target in range(particle_count):
        for source in range(particle_count):
            if target == source:
                continue
            rows = slice(target * active_count, (target + 1) * active_count)
            columns = slice(source * active_count, (source + 1) * active_count)
            displacement = tuple(float(value) for value in positions[source] - positions[target])
            translation[rows, columns] = _translation_block(
                wave_number, displacement, active_modes
            )

    external_local = nodal_standing_wave_coefficients(order)
    external_active = external_local[np.asarray(active_indices)]
    rhs = np.tile(external_active, particle_count)
    system = np.eye(particle_count * active_count, dtype=complex)
    system -= translation * global_scattering[None, :]
    effective_active = np.linalg.solve(system, rhs)
    scattered_active = global_scattering * effective_active
    residual = float(
        np.linalg.norm(system @ effective_active - rhs)
        / max(float(np.linalg.norm(rhs)), np.finfo(float).eps)
    )

    external = np.tile(external_local, (particle_count, 1))
    effective = np.zeros((particle_count, full_count), dtype=complex)
    scattered = np.zeros_like(effective)
    active_array = np.asarray(active_indices)
    effective[:, active_array] = effective_active.reshape(particle_count, active_count)
    scattered[:, active_array] = scattered_active.reshape(particle_count, active_count)
    return MieMultiparticleSolution(
        effective_incident_coefficients=effective,
        scattered_coefficients=scattered,
        external_incident_coefficients=external,
        scattering_coefficients=scattering_by_ell,
        system_matrix=system,
        right_hand_side=rhs,
        residual_relative=residual,
        condition_number=float(np.linalg.cond(system)),
        modes=full_modes,
        active_modes=active_modes,
        active_mode_indices=active_indices,
        lmax=order,
        used_planar_symmetry=bool(use_planar_symmetry),
    )
