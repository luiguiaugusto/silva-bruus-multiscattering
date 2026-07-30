"""Exact-Mie global multiple scattering with the complete radiation force."""

from __future__ import annotations

from dataclasses import dataclass
import numbers

import numpy as np

from .complete_force import complete_radiation_force_from_bsc
from .mie_multiparticle import MieMultiparticleSolution, solve_mie_multiparticle_nodal


@dataclass(frozen=True)
class ModelENodalResult:
    """Model-E coefficients, force channels, and decomposition diagnostic."""

    solution: MieMultiparticleSolution
    total_forces_xyz: np.ndarray
    external_forces_xyz: np.ndarray
    interaction_forces_xyz: np.ndarray
    external_scattered_forces_xyz: np.ndarray
    scattered_scattered_forces_xyz: np.ndarray
    scattered_incident_coefficients: np.ndarray
    decomposition_residual: float
    lmax: int


def _energy(value: object) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)) or not isinstance(value, numbers.Real):
        raise TypeError("energy_density must be a real scalar")
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError("energy_density must be finite and non-negative")
    return result


def solve_model_e_nodal(
    positions_xyz: object,
    k: object,
    radius: object,
    energy_density: object,
    f0: object,
    f1: object,
    lmax: object,
    *,
    use_planar_symmetry: bool = True,
) -> ModelENodalResult:
    """Solve exact-Mie Model E and separate complete force contributions."""

    if isinstance(lmax, (bool, np.bool_)) or not isinstance(lmax, numbers.Integral):
        raise TypeError("lmax must be an integer")
    if int(lmax) < 2:
        raise ValueError("Model E requires lmax >= 2 for the complete force")
    energy = _energy(energy_density)
    solution = solve_mie_multiparticle_nodal(
        positions_xyz, k, radius, f0, f1, int(lmax),
        use_planar_symmetry=use_planar_symmetry,
    )
    particle_count = len(solution.effective_incident_coefficients)
    total = np.empty((particle_count, 3), dtype=float)
    external = np.empty_like(total)
    scattered_scattered = np.empty_like(total)
    scattered_incident = solution.effective_incident_coefficients - solution.external_incident_coefficients
    for particle in range(particle_count):
        total[particle] = complete_radiation_force_from_bsc(
            solution.effective_incident_coefficients[particle], solution.scattering_coefficients, k, energy
        )
        external[particle] = complete_radiation_force_from_bsc(
            solution.external_incident_coefficients[particle], solution.scattering_coefficients, k, energy
        )
        scattered_scattered[particle] = complete_radiation_force_from_bsc(
            scattered_incident[particle], solution.scattering_coefficients, k, energy
        )
    interaction = total - external
    external_scattered = interaction - scattered_scattered
    reconstruction = external + external_scattered + scattered_scattered
    residual = float(
        np.linalg.norm(total - reconstruction)
        / max(float(np.linalg.norm(total)), np.finfo(float).eps)
    )
    return ModelENodalResult(
        solution=solution,
        total_forces_xyz=total,
        external_forces_xyz=external,
        interaction_forces_xyz=interaction,
        external_scattered_forces_xyz=external_scattered,
        scattered_scattered_forces_xyz=scattered_scattered,
        scattered_incident_coefficients=scattered_incident,
        decomposition_residual=residual,
        lmax=int(lmax),
    )
