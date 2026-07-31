"""Frozen response-blind utilities for the T14 scale-out validation."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from .cluster_families import linear_cluster
from .external_validation import (
    ExternalPredictionMetrics,
    ExternalThresholdAudit,
    LAMBDA_TARGETS,
    TOLERANCES,
)


FloatArray = NDArray[np.float64]
SCALE_OUT_PARTICLE_COUNTS = (15, 28)
SCALE_OUT_FAMILIES = ("linear", "compact", "irregular")
IRREGULAR_AMPLITUDE = 0.15
EXPECTED_SCALE_OUT_CASE_IDS = (
    "t14_n15_linear_level1",
    "t14_n15_linear_level2",
    "t14_n15_linear_level3",
    "t14_n15_linear_level4",
    "t14_n15_compact_level1",
    "t14_n15_compact_level2",
    "t14_n15_compact_level3",
    "t14_n15_compact_level4",
    "t14_n15_irregular_level1",
    "t14_n15_irregular_level2",
    "t14_n15_irregular_level3",
    "t14_n15_irregular_level4",
    "t14_n28_linear_level1",
    "t14_n28_linear_level2",
    "t14_n28_linear_level3",
    "t14_n28_linear_level4",
    "t14_n28_compact_level1",
    "t14_n28_compact_level2",
    "t14_n28_compact_level3",
    "t14_n28_compact_level4",
    "t14_n28_irregular_level1",
    "t14_n28_irregular_level2",
    "t14_n28_irregular_level3",
    "t14_n28_irregular_level4",
)


@dataclass(frozen=True)
class ScaleOutCase:
    """One analytically scaled, response-blind T14 geometry."""

    case_id: str
    scale_order: int
    particle_count: int
    family: str
    target_level: int
    lambda_target: float
    geometry_factor: float
    distance_ratio: float
    minimum_distance: float
    cluster_diameter: float
    positions_xyz: FloatArray


@dataclass(frozen=True)
class ScaleOutGateCriterion:
    """One literal T14 sufficiency or scientific gate criterion."""

    stage: str
    name: str
    observed: float
    threshold: float
    passed: bool
    justification: str


def _minimum_distance(positions: FloatArray) -> float:
    difference = positions[:, None, :] - positions[None, :, :]
    distances = np.linalg.norm(difference, axis=2)
    return float(np.min(distances[np.triu_indices(len(positions), 1)]))


def _cluster_diameter(positions: FloatArray) -> float:
    difference = positions[:, None, :] - positions[None, :, :]
    return float(np.max(np.linalg.norm(difference, axis=2)))


def _normalize(points_xy: FloatArray) -> FloatArray:
    points = np.asarray(points_xy, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2:
        raise ValueError("points_xy must have shape (N, 2), N >= 2")
    if not np.all(np.isfinite(points)):
        raise ValueError("points_xy must be finite")
    centered = points - np.mean(points, axis=0)
    positions = np.column_stack((centered, np.zeros(len(centered))))
    minimum = _minimum_distance(positions)
    if not np.isfinite(minimum) or minimum <= 0.0:
        raise ValueError("template points must be distinct")
    positions[:, :2] /= minimum
    positions[:, :2] -= np.mean(positions[:, :2], axis=0)
    positions[:, 2] = 0.0
    return positions


def triangular_compact_template(particle_count: int) -> FloatArray:
    """Return the centered unit-minimum triangular template."""

    if isinstance(particle_count, bool) or not isinstance(particle_count, int):
        raise TypeError("particle_count must be an integer")
    discriminant = 1 + 8 * particle_count
    q = (math.isqrt(discriminant) - 1) // 2
    if q < 2 or q * (q + 1) // 2 != particle_count:
        raise ValueError("particle_count must be a triangular number with q >= 2")
    points = [
        (column + row / 2.0, np.sqrt(3.0) * row / 2.0)
        for row in range(q)
        for column in range(q - row)
    ]
    return _normalize(np.asarray(points, dtype=float))


def irregular_scale_template(particle_count: int) -> FloatArray:
    """Return the prescribed deterministic perturbed triangular template."""

    compact = triangular_compact_template(particle_count)
    points = compact[:, :2].copy()
    for index in range(particle_count):
        frac_two = (index + 1) * np.sqrt(2.0)
        frac_two -= np.floor(frac_two)
        frac_three = (index + 1) * np.sqrt(3.0)
        frac_three -= np.floor(frac_three)
        points[index, 0] += IRREGULAR_AMPLITUDE * np.cos(2.0 * np.pi * frac_two)
        points[index, 1] += IRREGULAR_AMPLITUDE * np.sin(2.0 * np.pi * frac_three)
    return _normalize(points)


def scale_out_template(particle_count: int, family: str) -> FloatArray:
    """Build one of the deterministic unit-minimum T14 templates."""

    if particle_count not in SCALE_OUT_PARTICLE_COUNTS:
        raise ValueError("T14 particle_count must be 15 or 28")
    if family == "linear":
        result = linear_cluster(particle_count, 1.0)
    elif family == "compact":
        result = triangular_compact_template(particle_count)
    elif family == "irregular":
        result = irregular_scale_template(particle_count)
    else:
        raise ValueError("family must be linear, compact, or irregular")
    result = np.asarray(result, dtype=float)
    result[:, 2] = 0.0
    return result


def geometric_coupling_sum(positions_xyz: object) -> float:
    """Return ``max_i sum_(j!=i) r_ij^-3`` for finite distinct centers."""

    positions = np.asarray(positions_xyz, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2:
        raise ValueError("positions_xyz must have shape (N, 3), N >= 2")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite")
    difference = positions[:, None, :] - positions[None, :, :]
    distances = np.linalg.norm(difference, axis=2)
    off_diagonal = ~np.eye(len(positions), dtype=bool)
    if np.any(distances[off_diagonal] <= 0.0):
        raise ValueError("positions must be distinct")
    inverse_cube = np.zeros_like(distances)
    inverse_cube[off_diagonal] = distances[off_diagonal] ** -3
    return float(np.max(np.sum(inverse_cube, axis=1)))


def build_scale_out_cases(*, f1: float = 0.8) -> tuple[ScaleOutCase, ...]:
    """Construct the exact 24 T14 cases without response information."""

    contrast = float(f1)
    if not np.isfinite(contrast) or contrast != 0.8:
        raise ValueError("T14 freezes f1 at 0.8")
    cases = []
    order = 0
    for particle_count in SCALE_OUT_PARTICLE_COUNTS:
        for family in SCALE_OUT_FAMILIES:
            template = scale_out_template(particle_count, family)
            factor = geometric_coupling_sum(template)
            for level, target in enumerate(LAMBDA_TARGETS, start=1):
                order += 1
                scale = (abs(contrast) * factor / target) ** (1.0 / 3.0)
                positions = scale * template
                positions[:, 2] = 0.0
                cases.append(ScaleOutCase(
                    case_id=f"t14_n{particle_count}_{family}_level{level}",
                    scale_order=order,
                    particle_count=particle_count,
                    family=family,
                    target_level=level,
                    lambda_target=float(target),
                    geometry_factor=factor,
                    distance_ratio=float(scale),
                    minimum_distance=_minimum_distance(positions),
                    cluster_diameter=_cluster_diameter(positions),
                    positions_xyz=positions,
                ))
    if tuple(case.case_id for case in cases) != EXPECTED_SCALE_OUT_CASE_IDS:
        raise RuntimeError("T14 nominal case checksum changed")
    return tuple(cases)


def evaluate_scale_out_gate(
    *,
    eligible_count: int,
    eligible_by_n: Mapping[int, int],
    eligible_by_family: Mapping[str, int],
    eligible_by_level: Mapping[int, int],
    predicted_safe_eligible: Mapping[float, int],
    manifest_intact: bool,
    phase_a_integrity: bool,
    prior_integrity: bool,
    maximum_lmax: int,
    protocol_immutable: bool,
    resource_limit: bool,
    m1_global: ExternalPredictionMetrics | None,
    m1_by_n: Mapping[int, ExternalPredictionMetrics],
    m1_audits: Sequence[ExternalThresholdAudit],
) -> tuple[tuple[ScaleOutGateCriterion, ...], str, str]:
    """Apply sufficiency first and then the immutable scientific M1 gate."""

    sufficiency = [
        ScaleOutGateCriterion("sufficiency", "eligible_all", eligible_count, 20, eligible_count >= 20, "at least 20 of 24 cases eligible"),
    ]
    for particle_count in SCALE_OUT_PARTICLE_COUNTS:
        value = int(eligible_by_n.get(particle_count, 0))
        sufficiency.append(ScaleOutGateCriterion("sufficiency", f"eligible_n{particle_count}", value, 10, value >= 10, "at least 10 eligible cases for each N"))
    for family in SCALE_OUT_FAMILIES:
        value = int(eligible_by_family.get(family, 0))
        sufficiency.append(ScaleOutGateCriterion("sufficiency", f"eligible_{family}", value, 6, value >= 6, "at least six eligible cases for each family"))
    for level in range(1, 5):
        value = int(eligible_by_level.get(level, 0))
        sufficiency.append(ScaleOutGateCriterion("sufficiency", f"eligible_level_{level}", value, 5, value >= 5, "at least five eligible cases for each target level"))
    blind_counts = (6, 12, 18)
    for tolerance, required in zip(TOLERANCES, blind_counts):
        value = int(predicted_safe_eligible.get(tolerance, 0))
        sufficiency.append(ScaleOutGateCriterion("sufficiency", f"predicted_safe_eligible_{tolerance:g}", value, required, value == required, "all blindly predicted-safe cases are eligible"))
    sufficiency.extend((
        ScaleOutGateCriterion("sufficiency", "manifest_intact", float(manifest_intact), 1, manifest_intact, "manifest and frozen predictions are intact"),
        ScaleOutGateCriterion("sufficiency", "phase_a_integrity", float(phase_a_integrity), 1, phase_a_integrity, "four phase-A artifacts are byte-identical"),
        ScaleOutGateCriterion("sufficiency", "prior_integrity", float(prior_integrity), 1, prior_integrity, "all prior artifact hashes are preserved"),
        ScaleOutGateCriterion("sufficiency", "maximum_lmax", maximum_lmax, 13, maximum_lmax <= 13, "no solve exceeds Lmax=13"),
        ScaleOutGateCriterion("sufficiency", "protocol_immutable", float(protocol_immutable), 1, protocol_immutable, "no frozen rule changed after first push"),
    ))
    if not all(item.passed for item in sufficiency):
        decision = (
            "INCONCLUSIVE_T14_RESOURCE_LIMIT"
            if resource_limit
            else "INCONCLUSIVE_T14_INSUFFICIENT_MODEL_E_CONVERGENCE"
        )
        return tuple(sufficiency), decision, "HOLD_T15_T14_INCONCLUSIVE"

    if m1_global is None or any(n not in m1_by_n for n in SCALE_OUT_PARTICLE_COUNTS):
        raise ValueError("M1 metrics are required after sufficiency passes")
    global_audits = {
        round(audit.tolerance, 8): audit
        for audit in m1_audits
        if audit.model == "M1" and audit.scope == "all"
    }
    zero_false = len(global_audits) == 3 and all(
        audit.false_safe_count == 0 for audit in global_audits.values()
    )
    preserved_counts = all(
        round(tolerance, 8) in global_audits
        and global_audits[round(tolerance, 8)].predicted_safe_count == required
        for tolerance, required in zip(TOLERANCES, blind_counts)
    )
    n15 = m1_by_n[15]
    n28 = m1_by_n[28]
    scientific = (
        ScaleOutGateCriterion("scientific", "zero_false_safe", sum(a.false_safe_count for a in global_audits.values()), 0, zero_false, "zero false-safe cases at all tolerances"),
        ScaleOutGateCriterion("scientific", "blind_counts_6_12_18", sum(a.predicted_safe_count for a in global_audits.values()), 36, preserved_counts, "blind 6/12/18 safe counts are preserved"),
        ScaleOutGateCriterion("scientific", "global_rmse_log", m1_global.rmse_log, float(np.log(2.0)), m1_global.rmse_log <= np.log(2.0), "global RMSE log at most ln(2)"),
        ScaleOutGateCriterion("scientific", "n15_rmse_log", n15.rmse_log, float(np.log(2.0)), n15.rmse_log <= np.log(2.0), "N=15 RMSE log at most ln(2)"),
        ScaleOutGateCriterion("scientific", "n28_rmse_log", n28.rmse_log, float(np.log(2.0)), n28.rmse_log <= np.log(2.0), "N=28 RMSE log at most ln(2)"),
        ScaleOutGateCriterion("scientific", "global_factor_2", m1_global.fraction_within_factor_2, 0.80, m1_global.fraction_within_factor_2 >= 0.80, "global factor-two fraction at least 80%"),
        ScaleOutGateCriterion("scientific", "n15_factor_2", n15.fraction_within_factor_2, 0.75, n15.fraction_within_factor_2 >= 0.75, "N=15 factor-two fraction at least 75%"),
        ScaleOutGateCriterion("scientific", "n28_factor_2", n28.fraction_within_factor_2, 0.75, n28.fraction_within_factor_2 >= 0.75, "N=28 factor-two fraction at least 75%"),
        ScaleOutGateCriterion("scientific", "global_spearman", m1_global.spearman, 0.90, m1_global.spearman >= 0.90, "global Spearman at least 0.90"),
        ScaleOutGateCriterion("scientific", "frozen_predictor", float(protocol_immutable), 1, protocol_immutable, "predictor and protocol remain frozen"),
    )
    passed = all(item.passed for item in scientific)
    decision = (
        "PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX"
        if passed else "FAIL_T14_SCALE_OUT_FROZEN_LAMBDA_MAX"
    )
    next_gate = (
        "GO_T15_SYNTHESIS_AND_MANUSCRIPT"
        if passed else "GO_T15_SYNTHESIS_WITH_SCALE_BREAKDOWN"
    )
    return (*sufficiency, *scientific), decision, next_gate
