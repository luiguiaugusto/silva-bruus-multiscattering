"""Confirmatory grouped recalibration of rho1 against Model E."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .rho1_model_e_diagnostics import (
    fit_log_linear,
    spearman_correlation,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ConfirmatoryMetrics:
    """Pre-registered multiplicative metrics for positive predictions."""

    point_count: int
    rmse_log: float
    mae_log: float
    median_absolute_log_ratio: float
    fraction_within_factor_2: float
    fraction_within_factor_1_5: float
    spearman: float
    maximum_log_underestimation: float


@dataclass(frozen=True)
class LogoFoldFit:
    """One leave-one-group-out power-law fit."""

    held_out_group: str
    training_count: int
    test_count: int
    prefactor: float
    exponent: float


@dataclass(frozen=True)
class LogoPrediction:
    """One strictly out-of-fold prediction."""

    case_id: str
    held_out_group: str
    observed: float
    predictor: float
    predicted: float


@dataclass(frozen=True)
class SafetyClassification:
    """Safe/unsafe classification for one case and tolerance."""

    case_id: str
    group: str
    tolerance: float
    threshold: float
    predicted_safe: bool
    observed_safe: bool
    false_safe: bool
    false_unsafe: bool


@dataclass(frozen=True)
class SafetyAudit:
    """Aggregated out-of-fold classification at one tolerance."""

    tolerance: float
    predicted_safe_count: int
    predicted_safe_group_count: int
    true_safe_count: int
    false_safe_count: int
    false_unsafe_count: int
    worst_false_safe_excess: float
    coverage_sufficient: bool


@dataclass(frozen=True)
class BootstrapCalibration:
    """Group-bootstrap intervals for the final candidate calibration."""

    seed: int
    valid_samples: int
    attempts: int
    prefactor_interval: tuple[float, float]
    exponent_interval: tuple[float, float]
    threshold_intervals: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class GateCriterion:
    """One immutable pre-registered gate decision item."""

    name: str
    observed: float
    threshold: float
    passed: bool
    justification: str


def _positive_vector(values: ArrayLike, *, name: str) -> FloatArray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError(f"{name} must contain finite positive values")
    return result


def confirmatory_metrics(
    observed: ArrayLike,
    predicted: ArrayLike,
) -> ConfirmatoryMetrics:
    """Calculate the frozen T12.2 metrics for positive paired values."""

    actual = _positive_vector(observed, name="observed")
    estimate = _positive_vector(predicted, name="predicted")
    if actual.shape != estimate.shape or actual.size < 2:
        raise ValueError("observed and predicted must match and contain two points")
    residual = np.log(estimate) - np.log(actual)
    absolute = np.abs(residual)
    return ConfirmatoryMetrics(
        point_count=int(actual.size),
        rmse_log=float(np.sqrt(np.mean(residual * residual))),
        mae_log=float(np.mean(absolute)),
        median_absolute_log_ratio=float(np.median(absolute)),
        fraction_within_factor_2=float(np.mean(absolute <= np.log(2.0))),
        fraction_within_factor_1_5=float(np.mean(absolute <= np.log(1.5))),
        spearman=spearman_correlation(actual, estimate),
        maximum_log_underestimation=float(np.max(np.log(actual) - np.log(estimate))),
    )


def power_law_threshold(tolerance: float, prefactor: float, exponent: float) -> float:
    """Invert ``tolerance = prefactor * rho**exponent``."""

    values = np.asarray([tolerance, prefactor, exponent], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("threshold arguments must be finite")
    if tolerance <= 0.0 or prefactor <= 0.0 or exponent <= 0.0:
        raise ValueError("tolerance, prefactor and exponent must be positive")
    return float(np.exp((np.log(tolerance) - np.log(prefactor)) / exponent))


def logo_power_law_predictions(
    case_ids: Sequence[str],
    predictor: ArrayLike,
    observed: ArrayLike,
    groups: Sequence[str],
) -> tuple[tuple[LogoFoldFit, ...], tuple[LogoPrediction, ...]]:
    """Fit and predict a power law under deterministic LOGO validation."""

    identifiers = np.asarray(case_ids, dtype=str)
    labels = np.asarray(groups, dtype=str)
    x = _positive_vector(predictor, name="predictor")
    y = _positive_vector(observed, name="observed")
    if identifiers.ndim != 1 or labels.ndim != 1:
        raise ValueError("case_ids and groups must be one-dimensional")
    if not (identifiers.shape == labels.shape == x.shape == y.shape):
        raise ValueError("all LOGO inputs must have matching shapes")
    if len(set(identifiers.tolist())) != identifiers.size:
        raise ValueError("case_ids must be unique")
    if np.any(identifiers == "") or np.any(labels == ""):
        raise ValueError("case_ids and groups must be non-empty")
    unique_groups = tuple(sorted(set(labels.tolist())))
    if len(unique_groups) < 2:
        raise ValueError("at least two groups are required")
    fits: list[LogoFoldFit] = []
    predictions: list[LogoPrediction] = []
    for held_out in unique_groups:
        test = np.flatnonzero(labels == held_out)
        train = np.flatnonzero(labels != held_out)
        if test.size == 0 or train.size < 2:
            raise ValueError("each LOGO fold needs test data and two training points")
        train = np.asarray(
            sorted(train.tolist(), key=lambda item: identifiers[item]), dtype=int
        )
        fit = fit_log_linear(x[train], y[train])
        fits.append(LogoFoldFit(
            held_out_group=held_out,
            training_count=int(train.size),
            test_count=int(test.size),
            prefactor=fit.prefactor,
            exponent=fit.coefficient,
        ))
        for index in sorted(test.tolist(), key=lambda item: identifiers[item]):
            predictions.append(LogoPrediction(
                case_id=str(identifiers[index]),
                held_out_group=held_out,
                observed=float(y[index]),
                predictor=float(x[index]),
                predicted=float(fit.prefactor * x[index] ** fit.coefficient),
            ))
    predictions.sort(key=lambda item: item.case_id)
    return tuple(fits), tuple(predictions)


def classify_logo_safety(
    predictions: Sequence[LogoPrediction],
    fits: Sequence[LogoFoldFit],
    tolerances: Sequence[float] = (0.01, 0.05, 0.10),
) -> tuple[tuple[SafetyClassification, ...], tuple[SafetyAudit, ...]]:
    """Classify OOF safety using only the corresponding training-fold fit."""

    if not predictions or not fits:
        raise ValueError("predictions and fits must be non-empty")
    fit_by_group = {fit.held_out_group: fit for fit in fits}
    if len(fit_by_group) != len(fits):
        raise ValueError("held-out groups must be unique")
    tolerance_values = np.asarray(tolerances, dtype=float)
    if (
        tolerance_values.ndim != 1
        or tolerance_values.size == 0
        or not np.all(np.isfinite(tolerance_values))
        or np.any(tolerance_values <= 0.0)
        or len(set(tolerance_values.tolist())) != tolerance_values.size
    ):
        raise ValueError("tolerances must be unique, finite and positive")
    classifications: list[SafetyClassification] = []
    audits: list[SafetyAudit] = []
    for tolerance in tolerance_values:
        current: list[SafetyClassification] = []
        for prediction in predictions:
            if prediction.held_out_group not in fit_by_group:
                raise ValueError("a prediction has no corresponding fold fit")
            fit = fit_by_group[prediction.held_out_group]
            threshold = power_law_threshold(tolerance, fit.prefactor, fit.exponent)
            predicted_safe = prediction.predictor <= threshold
            observed_safe = prediction.observed <= tolerance
            current.append(SafetyClassification(
                case_id=prediction.case_id,
                group=prediction.held_out_group,
                tolerance=float(tolerance),
                threshold=threshold,
                predicted_safe=bool(predicted_safe),
                observed_safe=bool(observed_safe),
                false_safe=bool(predicted_safe and not observed_safe),
                false_unsafe=bool(not predicted_safe and observed_safe),
            ))
        predicted_safe_rows = [item for item in current if item.predicted_safe]
        false_safe_rows = [item for item in current if item.false_safe]
        audits.append(SafetyAudit(
            tolerance=float(tolerance),
            predicted_safe_count=len(predicted_safe_rows),
            predicted_safe_group_count=len({item.group for item in predicted_safe_rows}),
            true_safe_count=sum(item.predicted_safe and item.observed_safe for item in current),
            false_safe_count=len(false_safe_rows),
            false_unsafe_count=sum(item.false_unsafe for item in current),
            worst_false_safe_excess=max(
                (next(p.observed for p in predictions if p.case_id == item.case_id) - tolerance
                 for item in false_safe_rows),
                default=0.0,
            ),
            coverage_sufficient=(
                len(predicted_safe_rows) >= 3
                and len({item.group for item in predicted_safe_rows}) >= 2
            ),
        ))
        classifications.extend(current)
    classifications.sort(key=lambda item: (item.case_id, item.tolerance))
    return tuple(classifications), tuple(audits)


def grouped_bootstrap_calibration(
    predictor: ArrayLike,
    observed: ArrayLike,
    groups: Sequence[str],
    *,
    tolerances: Sequence[float] = (0.01, 0.05, 0.10),
    seed: int = 1202,
    valid_samples: int = 10_000,
    maximum_attempts: int = 100_000,
) -> BootstrapCalibration:
    """Return percentile intervals from whole-group bootstrap resampling."""

    x = _positive_vector(predictor, name="predictor")
    y = _positive_vector(observed, name="observed")
    labels = np.asarray(groups, dtype=str)
    tolerance_values = _positive_vector(tolerances, name="tolerances")
    if x.shape != y.shape or labels.ndim != 1 or labels.shape != x.shape:
        raise ValueError("bootstrap inputs must have matching one-dimensional shapes")
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    if valid_samples < 1 or maximum_attempts < valid_samples:
        raise ValueError("invalid bootstrap sample or attempt count")
    unique_groups = np.asarray(sorted(set(labels.tolist())), dtype=str)
    if unique_groups.size < 2:
        raise ValueError("at least two groups are required")
    indices = {group: np.flatnonzero(labels == group) for group in unique_groups}
    rng = np.random.default_rng(int(seed))
    prefactors: list[float] = []
    exponents: list[float] = []
    thresholds: list[list[float]] = [[] for _ in tolerance_values]
    attempts = 0
    while len(prefactors) < valid_samples and attempts < maximum_attempts:
        attempts += 1
        sampled = rng.choice(unique_groups, size=unique_groups.size, replace=True)
        selected = np.concatenate([indices[str(group)] for group in sampled])
        if np.ptp(np.log(x[selected])) == 0.0:
            continue
        fit = fit_log_linear(x[selected], y[selected])
        if not np.isfinite(fit.prefactor) or not np.isfinite(fit.coefficient):
            continue
        if fit.prefactor <= 0.0 or fit.coefficient <= 0.0:
            continue
        values = [
            power_law_threshold(value, fit.prefactor, fit.coefficient)
            for value in tolerance_values
        ]
        prefactors.append(fit.prefactor)
        exponents.append(fit.coefficient)
        for bucket, value in zip(thresholds, values):
            bucket.append(value)
    if len(prefactors) != valid_samples:
        raise RuntimeError("bootstrap did not produce the requested valid samples")
    interval = lambda values: tuple(float(item) for item in np.quantile(values, [0.025, 0.975]))
    return BootstrapCalibration(
        seed=int(seed),
        valid_samples=valid_samples,
        attempts=attempts,
        prefactor_interval=interval(prefactors),
        exponent_interval=interval(exponents),
        threshold_intervals=tuple(interval(values) for values in thresholds),
    )


def evaluate_recalibration_gate(
    candidate: ConfirmatoryMetrics,
    baseline: ConfirmatoryMetrics,
    fits: Sequence[LogoFoldFit],
    safety_audits: Sequence[SafetyAudit],
    *,
    predictions_finite_positive: bool,
    integrity_passed: bool,
) -> tuple[tuple[GateCriterion, ...], str]:
    """Evaluate the ten immutable T12.2 gate criteria."""

    coefficients_positive = bool(fits) and all(
        np.isfinite(fit.prefactor)
        and np.isfinite(fit.exponent)
        and fit.prefactor > 0.0
        and fit.exponent > 0.0
        for fit in fits
    )
    zero_false_safe = len(safety_audits) == 3 and all(
        audit.false_safe_count == 0 for audit in safety_audits
    )
    coverage = len(safety_audits) == 3 and all(
        audit.coverage_sufficient for audit in safety_audits
    )
    criteria = (
        GateCriterion("finite_positive_oof", float(predictions_finite_positive), 1.0, predictions_finite_positive, "all 28 OOF predictions must be finite and positive"),
        GateCriterion("positive_fold_coefficients", float(coefficients_positive), 1.0, coefficients_positive, "all seven fold prefactors and exponents must be positive"),
        GateCriterion("rmse_log", candidate.rmse_log, float(np.log(2.0)), candidate.rmse_log <= np.log(2.0), "candidate RMSE log must not exceed ln(2)"),
        GateCriterion("fraction_within_factor_2", candidate.fraction_within_factor_2, 0.85, candidate.fraction_within_factor_2 >= 0.85, "candidate factor-two coverage must be at least 85%"),
        GateCriterion("spearman", candidate.spearman, 0.90, candidate.spearman >= 0.90, "candidate Spearman must be at least 0.90"),
        GateCriterion("rmse_improves_p0", candidate.rmse_log, baseline.rmse_log, candidate.rmse_log < baseline.rmse_log, "candidate RMSE log must improve on frozen P0"),
        GateCriterion("factor_2_improves_p0", candidate.fraction_within_factor_2, baseline.fraction_within_factor_2, candidate.fraction_within_factor_2 > baseline.fraction_within_factor_2, "candidate factor-two coverage must improve on P0"),
        GateCriterion("zero_false_safe", float(sum(item.false_safe_count for item in safety_audits)), 0.0, zero_false_safe, "each tolerance must have zero false-safe cases"),
        GateCriterion("minimum_safety_coverage", float(min((item.predicted_safe_count for item in safety_audits), default=0)), 3.0, coverage, "each tolerance needs three predicted-safe cases in two groups"),
        GateCriterion("scientific_numeric_integrity", float(integrity_passed), 1.0, integrity_passed, "all scientific, numerical and integrity checks must pass"),
    )
    decision = (
        "GO_T13_WITH_RECALIBRATED_RHO1"
        if all(item.passed for item in criteria)
        else "NO_GO_T13_RHO1_NOT_QUANTITATIVE"
    )
    return criteria, decision
