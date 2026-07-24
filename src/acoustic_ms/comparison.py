"""Model A/B/C comparisons for planar clusters using approved T01 and T04 APIs."""

from dataclasses import dataclass

import numpy as np

from .force import RayleighNodalInteractionResult, solve_rayleigh_nodal_interaction_forces
from .silva_bruus import nodal_pair_force_on_probe


@dataclass(frozen=True)
class NodalForceModelComparison:
    """Per-particle Model A/B/C forces and their two- and three-body decomposition."""

    model_a_forces_xy: np.ndarray
    model_b_forces_xy: np.ndarray
    model_c_forces_xy: np.ndarray
    two_body_correction_xy: np.ndarray
    irreducible_multibody_xy: np.ndarray
    global_result: RayleighNodalInteractionResult


def compare_nodal_force_models(positions_xyz: object, k: float, radius: float, energy_density: float, f0: float, f1: float, lmax: int = 1) -> NodalForceModelComparison:
    """Compare SB pairwise A, isolated-pair MS B, and global MS C forces.

    Model B uses the same T04 observable as C, solving each unordered pair once.
    Hence C-B isolates collective rescattering within the fixed Lmax=1 basis.
    """
    global_result = solve_rayleigh_nodal_interaction_forces(positions_xyz, k, radius, energy_density, f0, f1, lmax=lmax)
    positions = np.asarray(positions_xyz, dtype=float)
    count = len(positions)
    model_a = np.zeros((count, 2))
    model_b = np.zeros((count, 2))
    for first in range(count):
        for second in range(first + 1, count):
            model_a[first] += nodal_pair_force_on_probe(positions[first, :2], positions[second, :2], k, radius, energy_density, f1)
            model_a[second] += nodal_pair_force_on_probe(positions[second, :2], positions[first, :2], k, radius, energy_density, f1)
            pair_result = solve_rayleigh_nodal_interaction_forces(positions[[first, second]], k, radius, energy_density, f0, f1, lmax=lmax)
            model_b[first] += pair_result.forces_xy[0]
            model_b[second] += pair_result.forces_xy[1]
    model_c = global_result.forces_xy
    return NodalForceModelComparison(model_a, model_b, model_c, model_b - model_a, model_c - model_b, global_result)
