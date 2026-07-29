"""Interaction forces and A--D comparisons for multipolar Model D."""

from dataclasses import dataclass

import numpy as np

from .comparison import NodalForceModelComparison, compare_nodal_force_models
from .multipolar_solver import MultipolarNodalSolution, solve_multipolar_nodal
from .translation import translation_matrix


@dataclass(frozen=True)
class MultipolarNodalInteractionResult:
    """Model-D solution, local reexpansion, and planar interaction forces."""

    solution: MultipolarNodalSolution
    local_scattered_coefficients: np.ndarray
    forces_xy: np.ndarray


@dataclass(frozen=True)
class NodalModelDComparison:
    """Models A--D evaluated without changing the established A--C APIs."""

    model_a_forces_xy: np.ndarray
    model_b_forces_xy: np.ndarray
    model_c_forces_xy: np.ndarray
    model_d_forces_xy: np.ndarray
    two_body_correction_xy: np.ndarray
    collective_rayleigh_correction_xy: np.ndarray
    multipolar_correction_xy: np.ndarray
    base_comparison: NodalForceModelComparison
    model_d_result: MultipolarNodalInteractionResult


def solve_multipolar_nodal_interaction_forces(
    positions_xyz: object,
    k: object,
    radius: object,
    energy_density: object,
    f0: object,
    f1: object,
    lmax: int = 1,
    *,
    use_planar_symmetry: bool = True,
) -> MultipolarNodalInteractionResult:
    """Solve Model D and evaluate the established external--scattered force."""
    try:
        energy_density = float(energy_density)
    except (TypeError, ValueError) as exc:
        raise ValueError("energy_density must be a real scalar") from exc
    if not np.isfinite(energy_density) or energy_density < 0.0:
        raise ValueError("energy_density must be finite and non-negative")

    solution = solve_multipolar_nodal(
        positions_xyz, k, radius, f0, f1, lmax,
        use_planar_symmetry=use_planar_symmetry,
    )
    positions = np.asarray(positions_xyz, dtype=float)
    local = np.zeros((len(positions), 9), dtype=complex)
    for target in range(len(positions)):
        for source in range(len(positions)):
            if source == target:
                continue
            local[target] += translation_matrix(
                float(k), positions[target], positions[source],
                target_lmax=2, source_lmax=lmax,
            ) @ solution.coefficients[source]

    b2m1 = local[:, 5]
    b21 = local[:, 7]
    prefactor = np.sqrt(30.0 * np.pi) * float(k) * float(radius) ** 3
    prefactor *= energy_density / 15.0
    forces = np.empty((len(positions), 2), dtype=float)
    forces[:, 0] = prefactor * np.real(float(f1) * (b2m1 - b21))
    forces[:, 1] = prefactor * np.real(-1j * float(f1) * (b21 + b2m1))
    return MultipolarNodalInteractionResult(solution, local, forces)


def compare_nodal_model_d(
    positions_xyz: object, k: object, radius: object,
    energy_density: object, f0: object, f1: object, lmax: int,
) -> NodalModelDComparison:
    """Compare the approved A--C hierarchy with multipolar Model D."""
    base = compare_nodal_force_models(
        positions_xyz, k, radius, energy_density, f0, f1, lmax=1
    )
    model_d = solve_multipolar_nodal_interaction_forces(
        positions_xyz, k, radius, energy_density, f0, f1, lmax
    )
    return NodalModelDComparison(
        base.model_a_forces_xy, base.model_b_forces_xy,
        base.model_c_forces_xy, model_d.forces_xy,
        base.two_body_correction_xy,
        base.model_c_forces_xy - base.model_b_forces_xy,
        model_d.forces_xy - base.model_c_forces_xy,
        base, model_d,
    )
