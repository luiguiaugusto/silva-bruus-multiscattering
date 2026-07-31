"""Frozen protocol utilities for T13 external validation.

The selection and prediction helpers consume metadata and predictors only.
They do not import or call any acoustic force solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import permutations
from typing import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .rho1_model_e_diagnostics import spearman_correlation


M1_PREFACTOR = 4.4964255121671126
M1_EXPONENT = 1.3883601043764593
M1_SAFETY_FACTOR = 2.5699703122019222
P3_PREFACTOR = 14.73950709797405
P3_EXPONENT = 1.4226504975598322
P3_SAFETY_FACTOR = 2.0464420079866286
TOLERANCES = (0.01, 0.05, 0.10)
LAMBDA_TARGETS = (
    0.0031111241226691642,
    0.011108933664494051,
    0.025457132710914911,
    0.065350897425260762,
)
LAMBDA_THRESHOLDS = (
    0.006222248245338328,
    0.019833411059191678,
    0.03267544871263038,
)
RHO_THRESHOLDS = (
    0.00358282366706918,
    0.011105567812764965,
    0.018077516446208253,
)
EXTERNAL_STRATA = (
    "n6_linear",
    "n6_compact",
    "n6_irregular",
    "n10_linear",
    "n10_compact",
    "n10_irregular",
)
EXPECTED_CASE_IDS = (
    "n6_linear_f0.1_d4.0",
    "n6_linear_f1.0_d6.0",
    "n6_linear_f0.1_d2.1",
    "n6_linear_f0.8_d3.0",
    "n6_compact_f0.8_d10.0",
    "n6_compact_f0.1_d3.0",
    "n6_compact_f0.4_d4.0",
    "n6_compact_f1.0_d4.0",
    "n6_irregular_f1.0_d10.0",
    "n6_irregular_f0.1_d3.0",
    "n6_irregular_f0.1_d2.5",
    "n6_irregular_f1.0_d4.0",
    "n10_linear_f0.1_d4.0",
    "n10_linear_f1.0_d6.0",
    "n10_linear_f0.1_d2.1",
    "n10_linear_f0.8_d3.0",
    "n10_compact_f0.1_d6.0",
    "n10_compact_f0.1_d4.0",
    "n10_compact_f0.1_d3.0",
    "n10_compact_f0.1_d2.1",
    "n10_irregular_f0.8_d10.0",
    "n10_irregular_f0.1_d3.0",
    "n10_irregular_f0.1_d2.5",
    "n10_irregular_f1.0_d4.0",
)


@dataclass(frozen=True)
class ExternalValidationCase:
    """One response-blind T08 holdout metadata record."""

    case_id: str
    particle_count: int
    family: str
    stratum: str
    f1: float
    distance_ratio: float
    lambda_max: float
    rho_l1: float
    reference_lmax: int
    target_level: int
    lambda_target: float


@dataclass(frozen=True)
class FrozenExternalPrediction:
    """Frozen point and conservative predictions for one case and model."""

    case_id: str
    model: str
    point_prediction: float
    safety_factor: float
    conservative_prediction: float
    safe_1pct: bool
    safe_5pct: bool
    safe_10pct: bool


@dataclass(frozen=True)
class ExternalPredictionMetrics:
    """Positive-response metrics frozen for the external analysis."""

    point_count: int
    rmse_log: float
    mae_log: float
    median_factor: float
    p90_factor: float
    maximum_factor: float
    fraction_within_factor_2: float
    fraction_within_factor_1_5: float
    spearman: float
    mean_log_bias: float
    median_log_bias: float


@dataclass(frozen=True)
class ExternalThresholdAudit:
    """Strict conservative classification at one external tolerance."""

    model: str
    scope: str
    tolerance: float
    eligible_count: int
    predicted_safe_count: int
    observed_safe_count: int
    false_safe_count: int
    false_unsafe_count: int
    safe_precision: float
    safe_coverage: float
    worst_predicted_safe_error: float
    false_safe_ids: tuple[str, ...]
    false_unsafe_ids: tuple[str, ...]


@dataclass(frozen=True)
class ExternalGateCriterion:
    """One literal T13 sufficiency or scientific gate item."""

    stage: str
    name: str
    observed: float
    threshold: float
    passed: bool
    justification: str


def _finite_positive(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite positive scalar") from exc
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return result


def _metadata_case(row: Mapping[str, object]) -> dict[str, object]:
    """Project a row onto the response-blind fields used by selection."""

    required = (
        "case_id", "split", "particle_count", "family", "f1",
        "distance_ratio", "lambda_max", "rho_l1", "reference_lmax",
        "total_converged",
    )
    if any(field not in row for field in required):
        raise ValueError("a response-blind metadata field is missing")
    case_id = str(row["case_id"])
    family = str(row["family"])
    particle_count = int(row["particle_count"])
    return {
        "case_id": case_id,
        "split": str(row["split"]),
        "particle_count": particle_count,
        "family": family,
        "stratum": f"n{particle_count}_{family}",
        "f1": _finite_positive(row["f1"], "f1"),
        "distance_ratio": _finite_positive(row["distance_ratio"], "distance_ratio"),
        "lambda_max": _finite_positive(row["lambda_max"], "lambda_max"),
        "rho_l1": _finite_positive(row["rho_l1"], "rho_l1"),
        "reference_lmax": int(row["reference_lmax"]),
        "total_converged": str(row["total_converged"]).lower(),
    }


def select_external_validation_cases(
    metadata_rows: Sequence[Mapping[str, object]],
) -> tuple[ExternalValidationCase, ...]:
    """Select four distinct cases per stratum by global log-distance cost."""

    projected = [_metadata_case(row) for row in metadata_rows]
    if len({row["case_id"] for row in projected}) != len(projected):
        raise ValueError("case_id values must be unique")
    selected: list[ExternalValidationCase] = []
    for stratum in EXTERNAL_STRATA:
        candidates = sorted(
            (
                row for row in projected
                if row["split"] == "holdout"
                and row["stratum"] == stratum
                and row["particle_count"] in (6, 10)
                and row["family"] in ("linear", "compact", "irregular")
                and row["total_converged"] == "true"
            ),
            key=lambda row: str(row["case_id"]),
        )
        if len(candidates) < 4:
            raise ValueError(f"stratum {stratum} has fewer than four valid cases")
        best_score = np.inf
        best_ids: tuple[str, ...] | None = None
        best_assignment: tuple[dict[str, object], ...] | None = None
        for assignment in permutations(candidates, len(LAMBDA_TARGETS)):
            score = float(sum(
                abs(np.log(float(row["lambda_max"])) - np.log(target))
                for row, target in zip(assignment, LAMBDA_TARGETS)
            ))
            identifiers = tuple(str(row["case_id"]) for row in assignment)
            if (
                score < best_score - 1.0e-15
                or (abs(score - best_score) <= 1.0e-15 and (best_ids is None or identifiers < best_ids))
            ):
                best_score = score
                best_ids = identifiers
                best_assignment = assignment
        if best_assignment is None:
            raise RuntimeError("external assignment unexpectedly failed")
        for level, (row, target) in enumerate(zip(best_assignment, LAMBDA_TARGETS), start=1):
            selected.append(ExternalValidationCase(
                case_id=str(row["case_id"]),
                particle_count=int(row["particle_count"]),
                family=str(row["family"]),
                stratum=stratum,
                f1=float(row["f1"]),
                distance_ratio=float(row["distance_ratio"]),
                lambda_max=float(row["lambda_max"]),
                rho_l1=float(row["rho_l1"]),
                reference_lmax=int(row["reference_lmax"]),
                target_level=level,
                lambda_target=float(target),
            ))
    if tuple(case.case_id for case in selected) != EXPECTED_CASE_IDS:
        raise RuntimeError("the T13 nominal selection checksum changed")
    return tuple(selected)


def frozen_external_predictions(
    case_id: str,
    lambda_max: object,
    rho_l1: object,
) -> tuple[FrozenExternalPrediction, FrozenExternalPrediction]:
    """Evaluate the immutable M1 and P3 laws and strict safety flags."""

    lam = _finite_positive(lambda_max, "lambda_max")
    rho = _finite_positive(rho_l1, "rho_l1")
    values = (
        ("M1", M1_PREFACTOR * lam**M1_EXPONENT, M1_SAFETY_FACTOR),
        ("P3", P3_PREFACTOR * rho**P3_EXPONENT, P3_SAFETY_FACTOR),
    )
    predictions = []
    for model, point, factor in values:
        conservative = factor * point
        flags = tuple(bool(conservative < tolerance) for tolerance in TOLERANCES)
        predictions.append(FrozenExternalPrediction(
            case_id=str(case_id), model=model, point_prediction=float(point),
            safety_factor=float(factor), conservative_prediction=float(conservative),
            safe_1pct=flags[0], safe_5pct=flags[1], safe_10pct=flags[2],
        ))
    return tuple(predictions)  # type: ignore[return-value]


def canonical_coordinate_hash(positions_xyz: ArrayLike) -> str:
    """Hash finite coordinates using the project's deterministic CSV precision."""

    positions = np.asarray(positions_xyz, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2:
        raise ValueError("positions_xyz must have shape (N, 3), N >= 2")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite")
    payload = ";".join(
        ":".join(format(float(value), ".17g") for value in row)
        for row in positions
    ).encode("ascii")
    return sha256(payload).hexdigest()


def successive_change(
    current: ArrayLike,
    previous: ArrayLike,
) -> tuple[float, bool, float]:
    """Return the scale-aware T11/T12 successive RMS change."""

    current_values = np.asarray(current, dtype=float)
    previous_values = np.asarray(previous, dtype=float)
    if current_values.shape != previous_values.shape or current_values.ndim != 2:
        raise ValueError("successive vector fields must have matching matrix shapes")
    if not np.all(np.isfinite(current_values)) or not np.all(np.isfinite(previous_values)):
        raise ValueError("successive vector fields must be finite")

    def rms(values: NDArray[np.float64]) -> float:
        return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))

    current_rms = rms(current_values)
    previous_rms = rms(previous_values)
    denominator = max(current_rms, previous_rms)
    absolute = rms(current_values - previous_values)
    tolerance = 128.0 * np.finfo(float).eps * denominator
    if denominator <= tolerance:
        return absolute, False, absolute
    return absolute / denominator, True, absolute


def minimum_two_step_confirmation(
    changes: Sequence[float],
    applicable: Sequence[bool],
    orders: Sequence[int],
    *,
    tolerance: float = 1.0e-5,
) -> int:
    """Return the first order closing two successive applicable changes."""

    values = np.asarray(changes, dtype=float)
    flags = np.asarray(applicable, dtype=bool)
    order_values = np.asarray(orders, dtype=int)
    if values.ndim != 1 or not (values.shape == flags.shape == order_values.shape):
        raise ValueError("confirmation inputs must be matching one-dimensional arrays")
    if not np.all(np.isfinite(values)) or not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("confirmation values and tolerance must be finite")
    for index in range(1, len(values)):
        if flags[index - 1] and flags[index] and values[index - 1] <= tolerance and values[index] <= tolerance:
            return int(order_values[index])
    return 0


def external_prediction_metrics(observed: ArrayLike, predicted: ArrayLike) -> ExternalPredictionMetrics:
    """Calculate the immutable positive external metrics without fitting."""

    actual = np.asarray(observed, dtype=float)
    estimate = np.asarray(predicted, dtype=float)
    if actual.ndim != 1 or actual.size < 2 or actual.shape != estimate.shape:
        raise ValueError("observed and predicted must be matching one-dimensional arrays")
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(estimate)) or np.any(actual <= 0.0) or np.any(estimate <= 0.0):
        raise ValueError("observed and predicted must be finite and positive")
    log_ratio = np.log(estimate / actual)
    factor = np.exp(np.abs(log_ratio))
    return ExternalPredictionMetrics(
        point_count=int(actual.size),
        rmse_log=float(np.sqrt(np.mean(log_ratio * log_ratio))),
        mae_log=float(np.mean(np.abs(log_ratio))),
        median_factor=float(np.median(factor)),
        p90_factor=float(np.percentile(factor, 90)),
        maximum_factor=float(np.max(factor)),
        fraction_within_factor_2=float(np.mean(factor <= 2.0)),
        fraction_within_factor_1_5=float(np.mean(factor <= 1.5)),
        spearman=spearman_correlation(actual, estimate),
        mean_log_bias=float(np.mean(log_ratio)),
        median_log_bias=float(np.median(log_ratio)),
    )


def audit_external_threshold(
    case_ids: Sequence[str],
    observed: ArrayLike,
    conservative_predictions: ArrayLike,
    *,
    model: str,
    scope: str,
    tolerance: float,
) -> ExternalThresholdAudit:
    """Apply strict predicted/observed safety comparisons."""

    identifiers = np.asarray(case_ids, dtype=str)
    actual = np.asarray(observed, dtype=float)
    conservative = np.asarray(conservative_predictions, dtype=float)
    if identifiers.ndim != 1 or not (identifiers.shape == actual.shape == conservative.shape):
        raise ValueError("threshold inputs must have matching one-dimensional shapes")
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(conservative)) or np.any(actual <= 0.0) or np.any(conservative <= 0.0):
        raise ValueError("threshold values must be finite and positive")
    threshold = _finite_positive(tolerance, "tolerance")
    predicted_safe = conservative < threshold
    observed_safe = actual < threshold
    false_safe = predicted_safe & ~observed_safe
    false_unsafe = ~predicted_safe & observed_safe
    predicted_count = int(np.sum(predicted_safe))
    observed_count = int(np.sum(observed_safe))
    true_safe = int(np.sum(predicted_safe & observed_safe))
    return ExternalThresholdAudit(
        model=model, scope=scope, tolerance=threshold,
        eligible_count=int(actual.size), predicted_safe_count=predicted_count,
        observed_safe_count=observed_count, false_safe_count=int(np.sum(false_safe)),
        false_unsafe_count=int(np.sum(false_unsafe)),
        safe_precision=float(true_safe / predicted_count) if predicted_count else 0.0,
        safe_coverage=float(true_safe / observed_count) if observed_count else 0.0,
        worst_predicted_safe_error=float(np.max(actual[predicted_safe])) if predicted_count else 0.0,
        false_safe_ids=tuple(sorted(identifiers[false_safe].tolist())),
        false_unsafe_ids=tuple(sorted(identifiers[false_unsafe].tolist())),
    )


def evaluate_external_validation_gate(
    *,
    eligible_count: int,
    eligible_n6: int,
    eligible_n10: int,
    diagnostics_all_passed: bool,
    manifest_intact: bool,
    integrity_passed: bool,
    m1_global: ExternalPredictionMetrics | None,
    m1_n6: ExternalPredictionMetrics | None,
    m1_n10: ExternalPredictionMetrics | None,
    m1_audits: Sequence[ExternalThresholdAudit],
    predicted_safe_n_by_tolerance: Mapping[float, tuple[int, int]],
) -> tuple[tuple[ExternalGateCriterion, ...], str, str]:
    """Evaluate sufficiency first, then the literal scientific M1 gate."""

    sufficiency = (
        ExternalGateCriterion("sufficiency", "eligible_all", float(eligible_count), 20.0, eligible_count >= 20, "at least 20 of 24 cases eligible"),
        ExternalGateCriterion("sufficiency", "eligible_n6", float(eligible_n6), 10.0, eligible_n6 >= 10, "at least 10 of 12 N=6 cases eligible"),
        ExternalGateCriterion("sufficiency", "eligible_n10", float(eligible_n10), 10.0, eligible_n10 >= 10, "at least 10 of 12 N=10 cases eligible"),
        ExternalGateCriterion("sufficiency", "diagnostics", float(diagnostics_all_passed), 1.0, diagnostics_all_passed, "all eligible cases pass numerical diagnostics"),
        ExternalGateCriterion("sufficiency", "manifest", float(manifest_intact), 1.0, manifest_intact, "all 24 frozen IDs present without substitution"),
        ExternalGateCriterion("sufficiency", "frozen_integrity", float(integrity_passed), 1.0, integrity_passed, "phase-A and 70 earlier artifacts remain immutable"),
    )
    if not all(item.passed for item in sufficiency):
        return sufficiency, "INCONCLUSIVE_T13_INSUFFICIENT_MODEL_E_CONVERGENCE", "HOLD_T14_MODEL_E_CONVERGENCE"
    if m1_global is None or m1_n6 is None or m1_n10 is None:
        raise ValueError("scientific metrics are required after sufficiency passes")
    by_tolerance = {round(item.tolerance, 8): item for item in m1_audits if item.scope == "all"}
    zero_false = len(by_tolerance) == 3 and all(item.false_safe_count == 0 for item in by_tolerance.values())
    minimum_safe = all(
        round(tolerance, 8) in by_tolerance and by_tolerance[round(tolerance, 8)].predicted_safe_count >= minimum
        for tolerance, minimum in zip(TOLERANCES, (3, 6, 9))
    )
    both_n = all(
        predicted_safe_n_by_tolerance.get(tolerance, (0, 0))[0] >= 1
        and predicted_safe_n_by_tolerance.get(tolerance, (0, 0))[1] >= 1
        for tolerance in (0.05, 0.10)
    )
    scientific = (
        ExternalGateCriterion("scientific", "zero_false_safe", float(sum(item.false_safe_count for item in by_tolerance.values())), 0.0, zero_false, "zero conservative false-safe cases at 1%, 5%, and 10%"),
        ExternalGateCriterion("scientific", "antivacuity_3_6_9", float(sum(item.predicted_safe_count for item in by_tolerance.values())), 18.0, minimum_safe, "at least 3, 6, and 9 predicted-safe eligible cases"),
        ExternalGateCriterion("scientific", "both_n_5_10pct", float(both_n), 1.0, both_n, "both N values represented among predicted safe at 5% and 10%"),
        ExternalGateCriterion("scientific", "global_rmse_log", m1_global.rmse_log, float(np.log(2.0)), m1_global.rmse_log <= np.log(2.0), "global RMSE log at most ln(2)"),
        ExternalGateCriterion("scientific", "global_factor_2", m1_global.fraction_within_factor_2, 0.80, m1_global.fraction_within_factor_2 >= 0.80, "global factor-two fraction at least 80%"),
        ExternalGateCriterion("scientific", "global_spearman", m1_global.spearman, 0.90, m1_global.spearman >= 0.90, "global Spearman at least 0.90"),
        ExternalGateCriterion("scientific", "n6_rmse_log", m1_n6.rmse_log, float(np.log(2.0)), m1_n6.rmse_log <= np.log(2.0), "N=6 RMSE log at most ln(2)"),
        ExternalGateCriterion("scientific", "n10_rmse_log", m1_n10.rmse_log, float(np.log(2.0)), m1_n10.rmse_log <= np.log(2.0), "N=10 RMSE log at most ln(2)"),
        ExternalGateCriterion("scientific", "n6_factor_2", m1_n6.fraction_within_factor_2, 0.75, m1_n6.fraction_within_factor_2 >= 0.75, "N=6 factor-two fraction at least 75%"),
        ExternalGateCriterion("scientific", "n10_factor_2", m1_n10.fraction_within_factor_2, 0.75, m1_n10.fraction_within_factor_2 >= 0.75, "N=10 factor-two fraction at least 75%"),
        ExternalGateCriterion("scientific", "protocol_immutable", float(integrity_passed), 1.0, integrity_passed, "no frozen coefficient, margin, threshold, case, or rule changed"),
    )
    passed = all(item.passed for item in scientific)
    decision = "PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX" if passed else "FAIL_T13_EXTERNAL_VALIDATION_LAMBDA_MAX"
    t14 = "GO_T14_SCALE_OUT_WITH_FROZEN_LAMBDA_MAX" if passed else "NO_GO_T14_LAMBDA_MAX_NOT_TRANSFERABLE"
    return (*sufficiency, *scientific), decision, t14
