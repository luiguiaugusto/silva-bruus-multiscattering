"""Nodal-plane interaction force from the T03 Rayleigh coupled field."""

from dataclasses import dataclass

import numpy as np

from .multipoles import mode_index
from .solver import RayleighNodalSolution, solve_rayleigh_nodal
from .translation import translation_matrix


@dataclass(frozen=True)
class RayleighNodalInteractionResult:
    """Auditable field and 2D interaction-force observables for the nodal plane."""

    solution: RayleighNodalSolution
    local_scattered_coefficients: np.ndarray
    forces_xy: np.ndarray


def solve_rayleigh_nodal_interaction_forces(
    positions_xyz: object,
    k: float,
    radius: float,
    energy_density: float,
    f0: float,
    f1: float,
    lmax: int = 1,
) -> RayleighNodalInteractionResult:
    """Return Eq. (27) nodal interaction forces in newtons.

    The T03 solver remains at ``lmax=1``.  The scattered field is reexpanded
    only for local force evaluation through ``ell=2``; self fields are omitted.
    This is the external--scattered cross-term interaction force, not a total
    off-nodal force and not a scattered--scattered approximation.
    """
    try:
        energy_density = float(energy_density)
    except (TypeError, ValueError) as exc:
        raise ValueError("energy_density must be a real finite scalar") from exc
    if not np.isfinite(energy_density) or energy_density < 0.0:
        raise ValueError("energy_density must be finite and non-negative")

    solution = solve_rayleigh_nodal(positions_xyz, k, radius, f0, f1, lmax=lmax)
    positions = np.asarray(positions_xyz, dtype=float)
    count = len(positions)
    local = np.zeros((count, 9), dtype=complex)
    for target in range(count):
        for source in range(count):
            if source != target:
                local[target] += translation_matrix(
                    k, positions[target], positions[source], target_lmax=2,
                    source_lmax=1,
                ) @ solution.coefficients[source]

    prefactor = np.sqrt(30.0 * np.pi) / 15.0 * k * radius**3 * energy_density
    b_minus = local[:, mode_index(2, -1)]
    b_plus = local[:, mode_index(2, 1)]
    forces = np.empty((count, 2), dtype=float)
    forces[:, 0] = prefactor * np.real(np.conj(f1) * (b_minus - b_plus))
    forces[:, 1] = prefactor * np.real(-1j * np.conj(f1) * (b_plus + b_minus))
    return RayleighNodalInteractionResult(solution, local, forces)
