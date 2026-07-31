"""Leakage-safe grouped validation of mechanistic Model-A error criteria.

This module contains only statistical post-processing of frozen, positive
error observations.  It never evaluates an acoustic force model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .rho1_model_e_diagnostics import spearman_correlation


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class MechanisticPowerLawFit:
    """Ordinary least-squares fit in logarithmic coordinates."""

    model: str
    point_count: int
    intercept: float
    prefactor: float
    alpha_lambda: float
    alpha_rho: float
    design_rank: int
    standardized_condition_number: float
    log_predictor_correlation: float


@dataclass(frozen=True)
class NestedLogoFold:
    """One outer LOGO fit and its training-only conservative margin."""

    model: str
    held_out_group: str
    training_count: int
    test_count: int
    inner_prediction_count: int
    fit: MechanisticPowerLawFit
    safety_factor: float
    maximum_inner_underprediction_log: float
    valid: bool


@dataclass(frozen=True)
class MechanisticOofPrediction:
    """One point and conservative prediction made strictly out of group."""

    model: str
    case_id: str
    held_out_group: str
    observed: float
    point_prediction: float
    safety_factor: float
    safe_prediction: float


@dataclass(frozen=True)
class MultiplicativeMetrics:
    """Metrics appropriate for positive multiplicative predictions."""

    point_count: int
    rmse_log: float
    mae_log: float
    fraction_within_factor_2: float
    fraction_within_factor_1_5: float
    spearman: float
    worst_multiplicative_ratio: float


@dataclass(frozen=True)
class ThresholdAudit:
    """Strict-threshold safety classification summary."""

    model: str
    rule: str
    tolerance: float
    predicted_safe_count: int
    observed_safe_count: int
    false_safe_count: int
    false_unsafe_count: int
    safe_precision: float
    safe_coverage: float
    false_safe_ids: tuple[str, ...]
    false_unsafe_ids: tuple[str, ...]


@dataclass(frozen=True)
class MechanisticGateCriterion:
    """One literal item in the pre-registered T12.3 decision gate."""

    candidate: str
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


def _validated_inputs(
    case_ids: Sequence[str],
    groups: Sequence[str],
    lambda_max: ArrayLike,
    rho_l1: ArrayLike,
    observed: ArrayLike,
) -> tuple[NDArray[np.str_], NDArray[np.str_], FloatArray, FloatArray, FloatArray]:
    identifiers = np.asarray(case_ids, dtype=str)
    labels = np.asarray(groups, dtype=str)
    lam = _positive_vector(lambda_max, name="lambda_max")
    rho = _positive_vector(rho_l1, name="rho_l1")
    target = _positive_vector(observed, name="observed")
    if identifiers.ndim != 1 or labels.ndim != 1:
        raise ValueError("case_ids and groups must be one-dimensional")
    if not (identifiers.shape == labels.shape == lam.shape == rho.shape == target.shape):
        raise ValueError("all validation inputs must have matching shapes")
    if np.any(identifiers == "") or np.any(labels == ""):
        raise ValueError("case_ids and groups must be non-empty")
    if len(set(identifiers.tolist())) != identifiers.size:
        raise ValueError("case_ids must be unique")
    if len(set(labels.tolist())) < 3:
        raise ValueError("at least three groups are required for nested LOGO")
    return identifiers, labels, lam, rho, target


def fit_mechanistic_power_law(
    lambda_max: ArrayLike,
    observed: ArrayLike,
    rho_l1: ArrayLike | None = None,
) -> MechanisticPowerLawFit:
    """Fit M1 or M2 by unweighted OLS in natural-log coordinates."""

    lam = _positive_vector(lambda_max, name="lambda_max")
    target = _positive_vector(observed, name="observed")
    if lam.shape != target.shape or lam.size < 3:
        raise ValueError("lambda_max and observed must match with at least three points")
    log_lambda = np.log(lam)
    columns = [np.ones(lam.size), log_lambda]
    model = "M1"
    log_rho: FloatArray | None = None
    if rho_l1 is not None:
        rho = _positive_vector(rho_l1, name="rho_l1")
        if rho.shape != lam.shape:
            raise ValueError("rho_l1 must match lambda_max")
        log_rho = np.log(rho)
        columns.append(log_rho)
        model = "M2"
    design = np.column_stack(columns)
    coefficients, _, rank, _ = np.linalg.lstsq(design, np.log(target), rcond=None)
    if int(rank) != design.shape[1]:
        raise ValueError("the logarithmic design matrix is singular")
    standardized = [np.ones(lam.size)]
    for predictor in columns[1:]:
        scale = float(np.std(predictor))
        if scale == 0.0:
            raise ValueError("logarithmic predictors must not be constant")
        standardized.append((predictor - np.mean(predictor)) / scale)
    condition = float(np.linalg.cond(np.column_stack(standardized)))
    correlation = (
        float(np.corrcoef(log_lambda, log_rho)[0, 1])
        if log_rho is not None
        else 0.0
    )
    return MechanisticPowerLawFit(
        model=model,
        point_count=int(lam.size),
        intercept=float(coefficients[0]),
        prefactor=float(np.exp(coefficients[0])),
        alpha_lambda=float(coefficients[1]),
        alpha_rho=float(coefficients[2]) if model == "M2" else 0.0,
        design_rank=int(rank),
        standardized_condition_number=condition,
        log_predictor_correlation=correlation,
    )


def predict_mechanistic_power_law(
    fit: MechanisticPowerLawFit,
    lambda_max: ArrayLike,
    rho_l1: ArrayLike | None = None,
) -> FloatArray:
    """Evaluate an M1 or M2 fit without clipping or a numerical floor."""

    lam = _positive_vector(lambda_max, name="lambda_max")
    log_prediction = fit.intercept + fit.alpha_lambda * np.log(lam)
    if fit.model == "M2":
        if rho_l1 is None:
            raise ValueError("rho_l1 is required for M2")
        rho = _positive_vector(rho_l1, name="rho_l1")
        if rho.shape != lam.shape:
            raise ValueError("rho_l1 must match lambda_max")
        log_prediction = log_prediction + fit.alpha_rho * np.log(rho)
    elif fit.model != "M1":
        raise ValueError("fit model must be M1 or M2")
    return np.exp(log_prediction)


def nested_logo_predictions(
    case_ids: Sequence[str],
    groups: Sequence[str],
    lambda_max: ArrayLike,
    rho_l1: ArrayLike,
    observed: ArrayLike,
    *,
    model: str,
) -> tuple[tuple[NestedLogoFold, ...], tuple[MechanisticOofPrediction, ...]]:
    """Return outer OOF predictions with inner-LOGO safety calibration."""

    if model not in {"M1", "M2"}:
        raise ValueError("model must be M1 or M2")
    identifiers, labels, lam, rho, target = _validated_inputs(
        case_ids, groups, lambda_max, rho_l1, observed
    )
    unique_groups = tuple(sorted(set(labels.tolist())))
    folds: list[NestedLogoFold] = []
    predictions: list[MechanisticOofPrediction] = []
    for outer_group in unique_groups:
        outer_test = np.flatnonzero(labels == outer_group)
        outer_train = np.flatnonzero(labels != outer_group)
        outer_train = np.asarray(
            sorted(outer_train.tolist(), key=lambda index: identifiers[index]), dtype=int
        )
        outer_fit = fit_mechanistic_power_law(
            lam[outer_train], target[outer_train], rho[outer_train] if model == "M2" else None
        )
        inner_predictions: dict[int, float] = {}
        for inner_group in sorted(set(labels[outer_train].tolist())):
            inner_test = np.flatnonzero((labels != outer_group) & (labels == inner_group))
            inner_train = np.flatnonzero((labels != outer_group) & (labels != inner_group))
            inner_train = np.asarray(
                sorted(inner_train.tolist(), key=lambda index: identifiers[index]), dtype=int
            )
            inner_fit = fit_mechanistic_power_law(
                lam[inner_train], target[inner_train],
                rho[inner_train] if model == "M2" else None,
            )
            values = predict_mechanistic_power_law(
                inner_fit, lam[inner_test], rho[inner_test] if model == "M2" else None
            )
            for index, value in zip(inner_test.tolist(), values.tolist()):
                if index in inner_predictions:
                    raise RuntimeError("an inner case received more than one prediction")
                inner_predictions[index] = float(value)
        if set(inner_predictions) != set(outer_train.tolist()):
            raise RuntimeError("inner LOGO did not predict every outer-training case exactly once")
        inner_order = outer_train.tolist()
        inner_values = np.asarray([inner_predictions[index] for index in inner_order])
        maximum_underprediction = float(
            np.max(np.log(target[inner_order]) - np.log(inner_values))
        )
        safety_factor = float(np.exp(maximum_underprediction))
        folds.append(NestedLogoFold(
            model=model,
            held_out_group=outer_group,
            training_count=int(outer_train.size),
            test_count=int(outer_test.size),
            inner_prediction_count=len(inner_predictions),
            fit=outer_fit,
            safety_factor=safety_factor,
            maximum_inner_underprediction_log=maximum_underprediction,
            valid=bool(np.isfinite(safety_factor) and safety_factor > 0.0),
        ))
        point = predict_mechanistic_power_law(
            outer_fit, lam[outer_test], rho[outer_test] if model == "M2" else None
        )
        for index, value in sorted(
            zip(outer_test.tolist(), point.tolist()), key=lambda item: identifiers[item[0]]
        ):
            predictions.append(MechanisticOofPrediction(
                model=model,
                case_id=str(identifiers[index]),
                held_out_group=outer_group,
                observed=float(target[index]),
                point_prediction=float(value),
                safety_factor=safety_factor,
                safe_prediction=float(safety_factor * value),
            ))
    predictions.sort(key=lambda item: item.case_id)
    if len(predictions) != len(identifiers) or len({item.case_id for item in predictions}) != len(identifiers):
        raise RuntimeError("outer LOGO must produce exactly one prediction per case")
    return tuple(folds), tuple(predictions)


def fixed_baseline_nested_predictions(
    case_ids: Sequence[str],
    groups: Sequence[str],
    observed: ArrayLike,
    point_predictions: ArrayLike,
    *,
    model: str,
) -> tuple[tuple[NestedLogoFold, ...], tuple[MechanisticOofPrediction, ...]]:
    """Apply training-only safety margins to a frozen, never-refitted baseline."""

    identifiers = np.asarray(case_ids, dtype=str)
    labels = np.asarray(groups, dtype=str)
    target = _positive_vector(observed, name="observed")
    point = _positive_vector(point_predictions, name="point_predictions")
    if not (identifiers.shape == labels.shape == target.shape == point.shape):
        raise ValueError("fixed baseline inputs must have matching shapes")
    if len(set(identifiers.tolist())) != identifiers.size:
        raise ValueError("case_ids must be unique")
    folds: list[NestedLogoFold] = []
    predictions: list[MechanisticOofPrediction] = []
    dummy = MechanisticPowerLawFit(model, 0, 0.0, 1.0, 0.0, 0.0, 0, 0.0, 0.0)
    for outer_group in sorted(set(labels.tolist())):
        train = np.flatnonzero(labels != outer_group)
        test = np.flatnonzero(labels == outer_group)
        maximum = float(np.max(np.log(target[train]) - np.log(point[train])))
        safety_factor = float(np.exp(maximum))
        folds.append(NestedLogoFold(
            model, outer_group, int(train.size), int(test.size), int(train.size),
            dummy, safety_factor, maximum, True,
        ))
        for index in sorted(test.tolist(), key=lambda item: identifiers[item]):
            predictions.append(MechanisticOofPrediction(
                model, str(identifiers[index]), outer_group, float(target[index]),
                float(point[index]), safety_factor, float(point[index] * safety_factor),
            ))
    predictions.sort(key=lambda item: item.case_id)
    return tuple(folds), tuple(predictions)


def multiplicative_metrics(observed: ArrayLike, predicted: ArrayLike) -> MultiplicativeMetrics:
    """Calculate frozen positive-response metrics in natural-log space."""

    target = _positive_vector(observed, name="observed")
    estimate = _positive_vector(predicted, name="predicted")
    if target.shape != estimate.shape or target.size < 2:
        raise ValueError("observed and predicted must match with at least two points")
    residual = np.log(estimate) - np.log(target)
    absolute = np.abs(residual)
    return MultiplicativeMetrics(
        point_count=int(target.size),
        rmse_log=float(np.sqrt(np.mean(residual * residual))),
        mae_log=float(np.mean(absolute)),
        fraction_within_factor_2=float(np.mean(absolute <= np.log(2.0))),
        fraction_within_factor_1_5=float(np.mean(absolute <= np.log(1.5))),
        spearman=spearman_correlation(target, estimate),
        worst_multiplicative_ratio=float(np.exp(np.max(absolute))),
    )


def audit_safety_thresholds(
    case_ids: Sequence[str],
    observed: ArrayLike,
    predictions: ArrayLike,
    *,
    model: str,
    rule: str,
    tolerances: Sequence[float] = (0.01, 0.05, 0.10),
) -> tuple[ThresholdAudit, ...]:
    """Audit strict ``prediction < tolerance`` safety classifications."""

    identifiers = np.asarray(case_ids, dtype=str)
    target = _positive_vector(observed, name="observed")
    estimate = _positive_vector(predictions, name="predictions")
    if identifiers.ndim != 1 or not (identifiers.shape == target.shape == estimate.shape):
        raise ValueError("safety inputs must have matching one-dimensional shapes")
    results: list[ThresholdAudit] = []
    for tolerance in _positive_vector(tolerances, name="tolerances"):
        predicted_safe = estimate < tolerance
        observed_safe = target < tolerance
        false_safe = predicted_safe & ~observed_safe
        false_unsafe = ~predicted_safe & observed_safe
        predicted_count = int(np.sum(predicted_safe))
        observed_count = int(np.sum(observed_safe))
        true_positive = int(np.sum(predicted_safe & observed_safe))
        results.append(ThresholdAudit(
            model=model,
            rule=rule,
            tolerance=float(tolerance),
            predicted_safe_count=predicted_count,
            observed_safe_count=observed_count,
            false_safe_count=int(np.sum(false_safe)),
            false_unsafe_count=int(np.sum(false_unsafe)),
            safe_precision=float(true_positive / predicted_count) if predicted_count else 0.0,
            safe_coverage=float(true_positive / observed_count) if observed_count else 0.0,
            false_safe_ids=tuple(sorted(identifiers[false_safe].tolist())),
            false_unsafe_ids=tuple(sorted(identifiers[false_unsafe].tolist())),
        ))
    return tuple(results)


def evaluate_mechanistic_gate(
    m1_metrics: MultiplicativeMetrics,
    m2_metrics: MultiplicativeMetrics,
    m1_audits: Sequence[ThresholdAudit],
    m2_audits: Sequence[ThresholdAudit],
    m1_folds: Sequence[NestedLogoFold],
    m2_folds: Sequence[NestedLogoFold],
    m1_full_fit: MechanisticPowerLawFit,
    m2_full_fit: MechanisticPowerLawFit,
    *,
    integrity_passed: bool,
    m2_unstable_collinearity: bool,
) -> tuple[tuple[MechanisticGateCriterion, ...], str, bool, bool]:
    """Apply the literal hierarchical T12.3 gate without adaptive tuning."""

    minima = {0.01: 3, 0.05: 8, 0.10: 12}

    def common(
        candidate: str,
        metrics: MultiplicativeMetrics,
        audits: Sequence[ThresholdAudit],
        folds: Sequence[NestedLogoFold],
        full_fit: MechanisticPowerLawFit,
    ) -> list[MechanisticGateCriterion]:
        audit_by_tolerance = {round(item.tolerance, 8): item for item in audits}
        zero_false = len(audits) == 3 and all(item.false_safe_count == 0 for item in audits)
        antivacuous = all(
            round(tolerance, 8) in audit_by_tolerance
            and audit_by_tolerance[round(tolerance, 8)].predicted_safe_count >= minimum
            for tolerance, minimum in minima.items()
        )
        valid_folds = len(folds) == 7 and all(item.valid for item in folds)
        median_lambda = float(np.median([item.fit.alpha_lambda for item in folds])) if folds else -np.inf
        return [
            MechanisticGateCriterion(candidate, "zero_false_safe", float(sum(item.false_safe_count for item in audits)), 0.0, zero_false, "zero conservative OOF false-safe cases at all tolerances"),
            MechanisticGateCriterion(candidate, "antivacuity_3_8_12", float(sum(item.predicted_safe_count for item in audits)), 23.0, antivacuous, "at least 3, 8 and 12 safe cases at 1%, 5% and 10%"),
            MechanisticGateCriterion(candidate, "rmse_log", metrics.rmse_log, float(np.log(2.0)), metrics.rmse_log <= np.log(2.0), "OOF point RMSE must not exceed ln(2)"),
            MechanisticGateCriterion(candidate, "within_factor_2", metrics.fraction_within_factor_2, 0.85, metrics.fraction_within_factor_2 >= 0.85, "at least 85% of OOF points within factor two"),
            MechanisticGateCriterion(candidate, "spearman", metrics.spearman, 0.90, metrics.spearman >= 0.90, "OOF Spearman must be at least 0.90"),
            MechanisticGateCriterion(candidate, "outer_inner_folds_valid", float(valid_folds), 1.0, valid_folds, "all seven outer and nested inner LOGO folds valid"),
            MechanisticGateCriterion(candidate, "positive_lambda_coefficient", min(full_fit.alpha_lambda, median_lambda), 0.0, full_fit.alpha_lambda > 0.0 and median_lambda > 0.0, "full and median-fold lambda exponents positive"),
            MechanisticGateCriterion(candidate, "frozen_integrity_and_no_holdout", float(integrity_passed), 1.0, integrity_passed, "frozen data preserved and no N=6,10 use"),
        ]

    criteria = common("M1", m1_metrics, m1_audits, m1_folds, m1_full_fit)
    m1_pass = all(item.passed for item in criteria)
    m2_common = common("M2", m2_metrics, m2_audits, m2_folds, m2_full_fit)
    sign_fraction = float(np.mean([
        fold.fit.alpha_lambda > 0.0 and fold.fit.alpha_rho > 0.0 for fold in m2_folds
    ])) if m2_folds else 0.0
    rmse_improvement = (1.0 - m2_metrics.rmse_log / m1_metrics.rmse_log) if m1_metrics.rmse_log > 0.0 else 0.0
    m1_counts = {item.tolerance: item.predicted_safe_count for item in m1_audits}
    coverage_material = (
        all(item.false_safe_count == 0 for item in m2_audits)
        and all(item.predicted_safe_count >= m1_counts.get(item.tolerance, 10**9) for item in m2_audits)
        and any(item.predicted_safe_count > m1_counts.get(item.tolerance, 10**9) for item in m2_audits)
    )
    extras = [
        MechanisticGateCriterion("M2", "identifiable_collinearity", float(not m2_unstable_collinearity), 1.0, not m2_unstable_collinearity, "M2 must not be UNSTABLE_COLLINEARITY"),
        MechanisticGateCriterion("M2", "stable_interpretable_signs", sign_fraction, 0.80, sign_fraction >= 0.80, "both mechanistic exponents positive in at least 80% of outer folds"),
        MechanisticGateCriterion("M2", "incremental_value", rmse_improvement, 0.05, rmse_improvement >= 0.05 or coverage_material, "at least 5% RMSE improvement or strict non-decreasing safe coverage with a gain"),
    ]
    criteria.extend(m2_common)
    criteria.extend(extras)
    m2_pass = all(item.passed for item in (*m2_common, *extras))
    if m1_pass:
        decision = "GO_T13_VALIDATE_LAMBDA_MAX"
    elif m2_pass:
        decision = "GO_T13_VALIDATE_LAMBDA_RHO"
    else:
        decision = "NO_GO_T13_NO_SAFE_LOW_N_CRITERION"
    return tuple(criteria), decision, m1_pass, m2_pass
