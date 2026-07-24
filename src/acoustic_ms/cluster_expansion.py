"""Connected body expansion for four-particle nodal Rayleigh clusters."""

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .comparison import NodalForceModelComparison, compare_nodal_force_models


@dataclass(frozen=True)
class NodalQuartetBodyExpansion:
    """Auditable two-, three-, and four-body decomposition for one quartet."""

    model_a_forces_xy: np.ndarray
    model_b_forces_xy: np.ndarray
    model_c_forces_xy: np.ndarray
    two_body_correction_xy: np.ndarray
    collective_correction_xy: np.ndarray
    triplet_indices: tuple[tuple[int, int, int], ...]
    irreducible_three_body_by_triplet_xy: np.ndarray
    irreducible_three_body_sum_xy: np.ndarray
    up_to_three_body_forces_xy: np.ndarray
    irreducible_four_body_xy: np.ndarray
    full_comparison: NodalForceModelComparison
    triplet_comparisons: tuple[NodalForceModelComparison, ...]


def decompose_nodal_quartet(
    positions_xyz: object,
    k: float,
    radius: float,
    energy_density: float,
    f0: float,
    f1: float,
    lmax: int = 1,
) -> NodalQuartetBodyExpansion:
    """Return the connected body expansion of exactly four planar particles."""
    if lmax != 1:
        raise ValueError("the T06 quartet expansion supports lmax=1 only")
    positions = np.asarray(positions_xyz, dtype=float)
    if positions.shape != (4, 3) or not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite with shape (4, 3)")

    full = compare_nodal_force_models(
        positions, k, radius, energy_density, f0, f1, lmax=lmax
    )
    triplet_indices = tuple(combinations(range(4), 3))
    triplet_comparisons = tuple(
        compare_nodal_force_models(
            positions[list(indices)], k, radius, energy_density, f0, f1, lmax=lmax
        )
        for indices in triplet_indices
    )

    by_triplet = np.zeros((4, 4, 2), dtype=float)
    for row, (indices, comparison) in enumerate(
        zip(triplet_indices, triplet_comparisons, strict=True)
    ):
        by_triplet[row, list(indices)] = comparison.irreducible_multibody_xy

    three_body_sum = np.sum(by_triplet, axis=0)
    up_to_three = full.model_b_forces_xy + three_body_sum
    four_body = full.model_c_forces_xy - up_to_three
    return NodalQuartetBodyExpansion(
        model_a_forces_xy=full.model_a_forces_xy,
        model_b_forces_xy=full.model_b_forces_xy,
        model_c_forces_xy=full.model_c_forces_xy,
        two_body_correction_xy=full.two_body_correction_xy,
        collective_correction_xy=full.model_c_forces_xy - full.model_b_forces_xy,
        triplet_indices=triplet_indices,
        irreducible_three_body_by_triplet_xy=by_triplet,
        irreducible_three_body_sum_xy=three_body_sum,
        up_to_three_body_forces_xy=up_to_three,
        irreducible_four_body_xy=four_body,
        full_comparison=full,
        triplet_comparisons=triplet_comparisons,
    )
