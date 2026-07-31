"""Frozen, response-blind utilities for the T14.1 large-N validation."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
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
from .scale_out_validation import (
    IRREGULAR_AMPLITUDE,
    geometric_coupling_sum,
    irregular_scale_template,
    triangular_compact_template,
)


FloatArray = NDArray[np.float64]
LARGE_N_PARTICLE_COUNTS = (45, 105)
LARGE_N_FAMILIES = ("linear", "compact", "irregular")
LARGE_N_TRIANGULAR_ROWS = {45: 9, 105: 14}
EXPECTED_LARGE_N_CASE_IDS = tuple(
    f"t14_1_n{particle_count}_{family}_level{level}"
    for particle_count in LARGE_N_PARTICLE_COUNTS
    for family in LARGE_N_FAMILIES
    for level in range(1, 5)
)


@dataclass(frozen=True)
class LargeNCase:
    """One analytically scaled T14.1 geometry selected without Model E."""

    case_id: str
    scale_order: int
    particle_count: int
    family: str
    target_level: int
    triangular_rows: int
    lambda_target: float
    geometry_factor: float
    distance_ratio: float
    minimum_distance: float
    cluster_diameter: float
    positions_xyz: FloatArray
    local_coupling: FloatArray


@dataclass(frozen=True)
class LargeNGateCriterion:
    """One literal sufficiency or scientific T14.1 gate criterion."""

    stage: str
    name: str
    observed: float
    threshold: float
    passed: bool
    justification: str


def _pair_distances(positions: FloatArray) -> FloatArray:
    difference = positions[:, None, :] - positions[None, :, :]
    return np.linalg.norm(difference, axis=2)


def _minimum_distance(positions: FloatArray) -> float:
    distances = _pair_distances(positions)
    return float(np.min(distances[np.triu_indices(len(positions), 1)]))


def _cluster_diameter(positions: FloatArray) -> float:
    return float(np.max(_pair_distances(positions)))


def large_n_template(particle_count: int, family: str) -> FloatArray:
    """Return a centered, planar, unit-minimum T14.1 template."""

    if particle_count not in LARGE_N_PARTICLE_COUNTS:
        raise ValueError("T14.1 particle_count must be 45 or 105")
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


def local_geometric_coupling(
    positions_xyz: object,
    radius: float,
    f1: float,
) -> FloatArray:
    """Return the frozen vector ``Lambda_i`` for every particle."""

    positions = np.asarray(positions_xyz, dtype=float)
    size = float(radius)
    contrast = float(f1)
    if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2:
        raise ValueError("positions_xyz must have shape (N, 3), N >= 2")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite")
    if not np.isfinite(size) or size <= 0.0:
        raise ValueError("radius must be finite and positive")
    if not np.isfinite(contrast):
        raise ValueError("f1 must be finite")
    distances = _pair_distances(positions)
    off_diagonal = ~np.eye(len(positions), dtype=bool)
    if np.any(distances[off_diagonal] <= 0.0):
        raise ValueError("positions must be distinct")
    if np.any(distances[off_diagonal] < 2.0 * size):
        raise ValueError("particles must not overlap")
    inverse_cube = np.zeros_like(distances)
    inverse_cube[off_diagonal] = (size / distances[off_diagonal]) ** 3
    return np.asarray(abs(contrast) * np.sum(inverse_cube, axis=1), dtype=float)


def local_coupling_statistics(values: object) -> dict[str, float | int | str]:
    """Return the immutable descriptive summary and byte checksum of ``Lambda_i``."""

    vector = np.asarray(values, dtype=float)
    if vector.ndim != 1 or vector.size < 2:
        raise ValueError("values must be a one-dimensional vector with N >= 2")
    if not np.all(np.isfinite(vector)) or np.any(vector < 0.0):
        raise ValueError("local coupling values must be finite and nonnegative")
    maximum = float(np.max(vector))
    payload = ";".join(format(float(value), ".17g") for value in vector)
    return {
        "lambda_min": float(np.min(vector)),
        "lambda_mean": float(np.mean(vector)),
        "lambda_median": float(np.median(vector)),
        "lambda_std": float(np.std(vector)),
        "lambda_p10": float(np.percentile(vector, 10.0)),
        "lambda_p90": float(np.percentile(vector, 90.0)),
        "lambda_max": maximum,
        "lambda_mean_over_max": 0.0 if maximum == 0.0 else float(np.mean(vector) / maximum),
        "fraction_ge_0_9_max": float(np.mean(vector >= 0.9 * maximum)),
        "first_argmax": int(np.argmax(vector)),
        "local_coupling_serialized": payload,
        "local_coupling_sha256": sha256(payload.encode("ascii")).hexdigest(),
    }


def build_large_n_cases(*, f1: float = 0.8) -> tuple[LargeNCase, ...]:
    """Construct the exact 24 T14.1 cases without consulting Model E."""

    contrast = float(f1)
    if not np.isfinite(contrast) or contrast != 0.8:
        raise ValueError("T14.1 freezes f1 at 0.8")
    cases: list[LargeNCase] = []
    order = 0
    for particle_count in LARGE_N_PARTICLE_COUNTS:
        for family in LARGE_N_FAMILIES:
            template = large_n_template(particle_count, family)
            factor = geometric_coupling_sum(template)
            for level, target in enumerate(LAMBDA_TARGETS, start=1):
                order += 1
                scale = (abs(contrast) * factor / target) ** (1.0 / 3.0)
                positions = np.asarray(scale * template, dtype=float)
                positions[:, 2] = 0.0
                local = local_geometric_coupling(positions, 1.0, contrast)
                cases.append(LargeNCase(
                    case_id=f"t14_1_n{particle_count}_{family}_level{level}",
                    scale_order=order,
                    particle_count=particle_count,
                    family=family,
                    target_level=level,
                    triangular_rows=LARGE_N_TRIANGULAR_ROWS[particle_count],
                    lambda_target=float(target),
                    geometry_factor=factor,
                    distance_ratio=float(scale),
                    minimum_distance=_minimum_distance(positions),
                    cluster_diameter=_cluster_diameter(positions),
                    positions_xyz=positions,
                    local_coupling=local,
                ))
    if tuple(case.case_id for case in cases) != EXPECTED_LARGE_N_CASE_IDS:
        raise RuntimeError("T14.1 nominal case checksum changed")
    return tuple(cases)


def classify_large_n_trend(ratios: object) -> str:
    """Apply the frozen N=105/N=45 deterioration classification."""

    values = np.asarray(ratios, dtype=float)
    applicable = values[np.isfinite(values) & (values > 0.0)]
    if applicable.size < 10:
        return "INCONCLUSIVE_LARGE_N_TREND"
    median = float(np.median(applicable))
    p90 = float(np.percentile(applicable, 90.0))
    if median <= 1.10 and p90 <= 1.25:
        return "NO_SYSTEMATIC_DETERIORATION"
    if median > 1.25 and int(np.sum(applicable > 1.10)) >= 9:
        return "SYSTEMATIC_DETERIORATION"
    return "MIXED_LARGE_N_TREND"


def evaluate_large_n_gate(
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
) -> tuple[tuple[LargeNGateCriterion, ...], str, str]:
    """Apply the literal preregistered T14.1 sufficiency and M1 gate."""

    criteria: list[LargeNGateCriterion] = [
        LargeNGateCriterion("sufficiency", "eligible_all", eligible_count, 20, eligible_count >= 20, "at least 20 of 24 cases eligible"),
    ]
    for particle_count in LARGE_N_PARTICLE_COUNTS:
        value = int(eligible_by_n.get(particle_count, 0))
        criteria.append(LargeNGateCriterion("sufficiency", f"eligible_n{particle_count}", value, 10, value >= 10, "at least 10 eligible cases for each N"))
    for family in LARGE_N_FAMILIES:
        value = int(eligible_by_family.get(family, 0))
        criteria.append(LargeNGateCriterion("sufficiency", f"eligible_{family}", value, 6, value >= 6, "at least six eligible cases for each family"))
    for level in range(1, 5):
        value = int(eligible_by_level.get(level, 0))
        criteria.append(LargeNGateCriterion("sufficiency", f"eligible_level_{level}", value, 5, value >= 5, "at least five eligible cases for each level"))
    for tolerance, required in zip(TOLERANCES, (6, 12, 18)):
        value = int(predicted_safe_eligible.get(tolerance, 0))
        criteria.append(LargeNGateCriterion("sufficiency", f"predicted_safe_eligible_{tolerance:g}", value, required, value == required, "all blindly predicted-safe cases are eligible"))
    criteria.extend((
        LargeNGateCriterion("sufficiency", "manifest_intact", float(manifest_intact), 1, manifest_intact, "manifest and frozen predictions are intact"),
        LargeNGateCriterion("sufficiency", "phase_a_integrity", float(phase_a_integrity), 1, phase_a_integrity, "five phase-A artifacts are byte-identical"),
        LargeNGateCriterion("sufficiency", "prior_integrity", float(prior_integrity), 1, prior_integrity, "all prior artifact hashes are preserved"),
        LargeNGateCriterion("sufficiency", "maximum_lmax", maximum_lmax, 13, maximum_lmax <= 13, "no solve exceeds Lmax=13"),
        LargeNGateCriterion("sufficiency", "protocol_immutable", float(protocol_immutable), 1, protocol_immutable, "no frozen rule changed after first push"),
    ))
    if not all(item.passed for item in criteria):
        decision = "INCONCLUSIVE_T14_1_RESOURCE_LIMIT" if resource_limit else "INCONCLUSIVE_T14_1_INSUFFICIENT_MODEL_E_CONVERGENCE"
        return tuple(criteria), decision, "HOLD_T15_T14_1_INCONCLUSIVE"
    if m1_global is None or any(n not in m1_by_n for n in LARGE_N_PARTICLE_COUNTS):
        raise ValueError("M1 metrics are required after sufficiency passes")
    global_audits = {
        round(audit.tolerance, 8): audit for audit in m1_audits
        if audit.model == "M1" and audit.scope == "all"
    }
    zero_false = len(global_audits) == 3 and all(a.false_safe_count == 0 for a in global_audits.values())
    preserved = all(
        round(tolerance, 8) in global_audits
        and global_audits[round(tolerance, 8)].predicted_safe_count == required
        for tolerance, required in zip(TOLERANCES, (6, 12, 18))
    )
    scientific = (
        LargeNGateCriterion("scientific", "zero_false_safe", sum(a.false_safe_count for a in global_audits.values()), 0, zero_false, "zero false-safe cases at all tolerances"),
        LargeNGateCriterion("scientific", "blind_counts_6_12_18", sum(a.predicted_safe_count for a in global_audits.values()), 36, preserved, "blind 6/12/18 safe counts are preserved"),
        LargeNGateCriterion("scientific", "global_rmse_log", m1_global.rmse_log, float(np.log(2.0)), m1_global.rmse_log <= np.log(2.0), "global RMSE log at most ln(2)"),
        *(LargeNGateCriterion("scientific", f"n{n}_rmse_log", m1_by_n[n].rmse_log, float(np.log(2.0)), m1_by_n[n].rmse_log <= np.log(2.0), "per-N RMSE log at most ln(2)") for n in LARGE_N_PARTICLE_COUNTS),
        LargeNGateCriterion("scientific", "global_factor_2", m1_global.fraction_within_factor_2, 0.80, m1_global.fraction_within_factor_2 >= 0.80, "global factor-two fraction at least 80%"),
        *(LargeNGateCriterion("scientific", f"n{n}_factor_2", m1_by_n[n].fraction_within_factor_2, 0.75, m1_by_n[n].fraction_within_factor_2 >= 0.75, "per-N factor-two fraction at least 75%") for n in LARGE_N_PARTICLE_COUNTS),
        LargeNGateCriterion("scientific", "global_spearman", m1_global.spearman, 0.90, m1_global.spearman >= 0.90, "global Spearman at least 0.90"),
        LargeNGateCriterion("scientific", "frozen_predictor", float(protocol_immutable), 1, protocol_immutable, "M1 and P3 remain frozen"),
    )
    passed = all(item.passed for item in scientific)
    decision = "PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX" if passed else "FAIL_T14_1_LARGE_N_FROZEN_LAMBDA_MAX"
    next_gate = "GO_T15_SYNTHESIS_AND_MANUSCRIPT" if passed else "GO_T15_SYNTHESIS_WITH_LARGE_N_BREAKDOWN"
    return (*criteria, *scientific), decision, next_gate


__all__ = [
    "EXPECTED_LARGE_N_CASE_IDS",
    "IRREGULAR_AMPLITUDE",
    "LARGE_N_FAMILIES",
    "LARGE_N_PARTICLE_COUNTS",
    "LARGE_N_TRIANGULAR_ROWS",
    "LargeNCase",
    "LargeNGateCriterion",
    "build_large_n_cases",
    "classify_large_n_trend",
    "evaluate_large_n_gate",
    "large_n_template",
    "local_coupling_statistics",
    "local_geometric_coupling",
]
