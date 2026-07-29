"""Connected subset expansion with a common multipolar truncation."""

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .model_d import MultipolarNodalInteractionResult, solve_multipolar_nodal_interaction_forces


@dataclass(frozen=True)
class MultipolarConnectedTerm:
    """One connected contribution embedded in global particle ordering."""

    indices: tuple[int, ...]
    forces_xy: np.ndarray
    subset_result: MultipolarNodalInteractionResult


@dataclass(frozen=True)
class MultipolarClusterExpansion:
    """Connected force expansion for two, three, or four particles."""

    model_d_forces_xy: np.ndarray
    two_body_sum_xy: np.ndarray
    irreducible_three_body_by_subset_xy: np.ndarray
    irreducible_three_body_sum_xy: np.ndarray
    irreducible_four_body_xy: np.ndarray
    connected_terms: tuple[MultipolarConnectedTerm, ...]
    subset_results: tuple[tuple[tuple[int, ...], MultipolarNodalInteractionResult], ...]
    full_result: MultipolarNodalInteractionResult
    lmax: int


def decompose_multipolar_cluster(
    positions_xyz: object, k: object, radius: object,
    energy_density: object, f0: object, f1: object, lmax: int,
) -> MultipolarClusterExpansion:
    """Resolve every subset once and apply vector inclusion--exclusion."""
    positions = np.asarray(positions_xyz, dtype=float)
    if positions.ndim != 2 or positions.shape[1:] != (3,):
        raise ValueError("positions_xyz must have shape (N, 3)")
    particle_count = len(positions)
    if particle_count < 2 or particle_count > 4:
        raise ValueError("the connected expansion supports 2 <= N <= 4")

    results: dict[tuple[int, ...], MultipolarNodalInteractionResult] = {}
    connected: dict[tuple[int, ...], np.ndarray] = {}
    terms: list[MultipolarConnectedTerm] = []
    for size in range(2, particle_count + 1):
        for indices in combinations(range(particle_count), size):
            result = solve_multipolar_nodal_interaction_forces(
                positions[np.asarray(indices)], k, radius, energy_density,
                f0, f1, lmax,
            )
            embedded = np.zeros((particle_count, 2), dtype=float)
            embedded[np.asarray(indices)] = result.forces_xy
            value = embedded.copy()
            for proper_indices, proper_value in connected.items():
                if len(proper_indices) < size and set(proper_indices).issubset(indices):
                    value -= proper_value
            results[indices] = result
            connected[indices] = value
            terms.append(MultipolarConnectedTerm(indices, value, result))

    pair_sum = sum(
        (value for indices, value in connected.items() if len(indices) == 2),
        np.zeros((particle_count, 2)),
    )
    triplet_values = [value for indices, value in connected.items() if len(indices) == 3]
    triplet_array = (
        np.stack(triplet_values)
        if triplet_values else np.empty((0, particle_count, 2), dtype=float)
    )
    triplet_sum = sum(triplet_values, np.zeros((particle_count, 2)))
    full_indices = tuple(range(particle_count))
    four_body = (
        connected[full_indices].copy()
        if particle_count == 4 else np.zeros((particle_count, 2))
    )
    return MultipolarClusterExpansion(
        results[full_indices].forces_xy, pair_sum, triplet_array, triplet_sum,
        four_body, tuple(terms), tuple(results.items()), results[full_indices], lmax,
    )
