"""Reusable diagnostics for the T12.1 Model-E/rho1 failure analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ApplicableScalar:
    """A scalar diagnostic together with its mathematical applicability."""

    value: float
    applicable: bool
    reason: str


@dataclass(frozen=True)
class LogLinearFit:
    """Unweighted fit of ``log(y) = intercept + slope * log(x)``."""

    point_count: int
    intercept: float
    coefficient: float
    prefactor: float


@dataclass(frozen=True)
class OutOfFoldMetrics:
    """Multiplicative prediction diagnostics in logarithmic space."""

    point_count: int
    rmse_log: float
    median_factor: float
    p90_factor: float
    maximum_factor: float
    fraction_within_factor_2: float
    spearman: float


@dataclass(frozen=True)
class ConvergenceTailDiagnostics:
    """Convergence classification and tail-ratio diagnostics for one channel."""

    confirmation_order: int | None
    last_change: float
    classification: str
    q_median: float
    q_minimum: float
    q_maximum: float
    q_count: int
    oscillatory: bool


@dataclass(frozen=True)
class MechanismDiagnostics:
    """Signed-vector diagnostics for ``C_D + C_M + C_S = C``."""

    rms_d: float
    rms_m: float
    rms_s: float
    rms_c: float
    mu_dm: ApplicableScalar
    mu_ds: ApplicableScalar
    mu_ms: ApplicableScalar
    mu_dc: ApplicableScalar
    mu_mc: ApplicableScalar
    mu_sc: ApplicableScalar
    p_d: ApplicableScalar
    p_m: ApplicableScalar
    p_s: ApplicableScalar
    projection_sum: ApplicableScalar
    r_s_over_d: ApplicableScalar
    r_m_over_d: ApplicableScalar
    cancellation_ratio: ApplicableScalar
    closure_rms: float


def _vectors(value: ArrayLike, *, name: str) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 2 or result.shape[1] != 3 or result.shape[0] == 0:
        raise ValueError(f"{name} must have shape (N, 3)")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _rms(value: FloatArray) -> float:
    return float(np.sqrt(np.mean(np.sum(value * value, axis=1))))


def vector_field_inner_product(first: ArrayLike, second: ArrayLike) -> float:
    """Return ``N^-1 sum_i first_i dot second_i``."""

    left = _vectors(first, name="first")
    right = _vectors(second, name="second")
    if left.shape != right.shape:
        raise ValueError("first and second must have matching shapes")
    return float(np.mean(np.sum(left * right, axis=1)))


def _tolerance(*values: FloatArray) -> float:
    scale = max((_rms(value) for value in values), default=0.0)
    return 128.0 * np.finfo(float).eps * scale


def vector_field_cosine(first: ArrayLike, second: ArrayLike) -> ApplicableScalar:
    """Return the signed cosine between two vector fields."""

    left = _vectors(first, name="first")
    right = _vectors(second, name="second")
    if left.shape != right.shape:
        raise ValueError("first and second must have matching shapes")
    left_rms = _rms(left)
    right_rms = _rms(right)
    tolerance = _tolerance(left, right)
    if left_rms <= tolerance or right_rms <= tolerance:
        return ApplicableScalar(0.0, False, "numerically_null_vector_field")
    value = vector_field_inner_product(left, right) / (left_rms * right_rms)
    excess = abs(value) - 1.0
    rounding = 64.0 * np.finfo(float).eps
    if excess > rounding:
        raise RuntimeError("cosine exceeds its physical interval beyond rounding")
    return ApplicableScalar(float(np.clip(value, -1.0, 1.0)), True, "applicable")


def vector_field_projection(component: ArrayLike, total: ArrayLike) -> ApplicableScalar:
    """Return the signed projection ``<component,total>/<total,total>``."""

    part = _vectors(component, name="component")
    whole = _vectors(total, name="total")
    if part.shape != whole.shape:
        raise ValueError("component and total must have matching shapes")
    denominator = vector_field_inner_product(whole, whole)
    tolerance = _tolerance(part, whole)
    if _rms(whole) <= tolerance:
        return ApplicableScalar(0.0, False, "numerically_null_total")
    return ApplicableScalar(
        vector_field_inner_product(part, whole) / denominator,
        True,
        "applicable",
    )


def vector_field_amplitude_ratio(
    numerator: ArrayLike,
    denominator: ArrayLike,
) -> ApplicableScalar:
    """Return the RMS-amplitude ratio without a dimensional absolute floor."""

    top = _vectors(numerator, name="numerator")
    bottom = _vectors(denominator, name="denominator")
    if top.shape != bottom.shape:
        raise ValueError("numerator and denominator must have matching shapes")
    tolerance = _tolerance(top, bottom)
    bottom_rms = _rms(bottom)
    if bottom_rms <= tolerance:
        return ApplicableScalar(0.0, False, "numerically_null_denominator")
    return ApplicableScalar(_rms(top) / bottom_rms, True, "applicable")


def mechanism_diagnostics(
    model_a: ArrayLike,
    model_d: ArrayLike,
    interaction: ArrayLike,
    external_scattered: ArrayLike,
    scattered_scattered: ArrayLike,
) -> MechanismDiagnostics:
    """Diagnose the signed decomposition ``E-A=(D-A)+(E_ext-sc-D)+E_ss``."""

    a = _vectors(model_a, name="model_a")
    d = _vectors(model_d, name="model_d")
    total_e = _vectors(interaction, name="interaction")
    ext = _vectors(external_scattered, name="external_scattered")
    ss = _vectors(scattered_scattered, name="scattered_scattered")
    if len({value.shape for value in (a, d, total_e, ext, ss)}) != 1:
        raise ValueError("all vector fields must have matching shapes")
    c_d = d - a
    c_m = ext - d
    c_s = ss
    correction = total_e - a
    closure = correction - c_d - c_m - c_s
    p_d = vector_field_projection(c_d, correction)
    p_m = vector_field_projection(c_m, correction)
    p_s = vector_field_projection(c_s, correction)
    if p_d.applicable and p_m.applicable and p_s.applicable:
        projection_sum = ApplicableScalar(
            p_d.value + p_m.value + p_s.value, True, "applicable"
        )
    else:
        projection_sum = ApplicableScalar(0.0, False, "numerically_null_total")
    component_sum = _rms(c_d) + _rms(c_m) + _rms(c_s)
    correction_rms = _rms(correction)
    tolerance = _tolerance(c_d, c_m, c_s, correction)
    cancellation = ApplicableScalar(
        component_sum / correction_rms
        if component_sum > tolerance and correction_rms > tolerance
        else 0.0,
        component_sum > tolerance and correction_rms > tolerance,
        "applicable"
        if component_sum > tolerance and correction_rms > tolerance
        else "numerically_null_correction",
    )
    return MechanismDiagnostics(
        rms_d=_rms(c_d),
        rms_m=_rms(c_m),
        rms_s=_rms(c_s),
        rms_c=correction_rms,
        mu_dm=vector_field_cosine(c_d, c_m),
        mu_ds=vector_field_cosine(c_d, c_s),
        mu_ms=vector_field_cosine(c_m, c_s),
        mu_dc=vector_field_cosine(c_d, correction),
        mu_mc=vector_field_cosine(c_m, correction),
        mu_sc=vector_field_cosine(c_s, correction),
        p_d=p_d,
        p_m=p_m,
        p_s=p_s,
        projection_sum=projection_sum,
        r_s_over_d=vector_field_amplitude_ratio(c_s, c_d),
        r_m_over_d=vector_field_amplitude_ratio(c_m, c_d),
        cancellation_ratio=cancellation,
        closure_rms=_rms(closure),
    )


def fit_log_linear(x: ArrayLike, y: ArrayLike) -> LogLinearFit:
    """Fit a deterministic one-feature log-linear power law."""

    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if x_values.ndim != 1 or y_values.ndim != 1:
        raise ValueError("x and y must be one-dimensional")
    if x_values.shape != y_values.shape or x_values.size < 2:
        raise ValueError("x and y must have matching shapes with at least two points")
    if not np.all(np.isfinite(x_values)) or not np.all(np.isfinite(y_values)):
        raise ValueError("x and y must contain only finite values")
    if np.any(x_values <= 0.0) or np.any(y_values <= 0.0):
        raise ValueError("x and y must be strictly positive")
    log_x = np.log(x_values)
    log_y = np.log(y_values)
    if np.ptp(log_x) == 0.0:
        raise ValueError("x must not be constant")
    design = np.column_stack((np.ones(x_values.size), log_x))
    intercept, coefficient = np.linalg.lstsq(design, log_y, rcond=None)[0]
    return LogLinearFit(
        int(x_values.size),
        float(intercept),
        float(coefficient),
        float(np.exp(intercept)),
    )


def leave_group_out_folds(groups: Sequence[str]) -> tuple[tuple[str, NDArray[np.int64], NDArray[np.int64]], ...]:
    """Return deterministic, exhaustive leave-one-group-out folds."""

    labels = np.asarray(groups, dtype=str)
    if labels.ndim != 1 or labels.size < 2:
        raise ValueError("groups must be one-dimensional with at least two entries")
    unique = tuple(dict.fromkeys(labels.tolist()))
    if len(unique) < 2:
        raise ValueError("at least two distinct groups are required")
    result = []
    for group in unique:
        test = np.flatnonzero(labels == group)
        train = np.flatnonzero(labels != group)
        if test.size == 0 or train.size < 2:
            raise ValueError("each fold needs a non-empty test set and two training points")
        result.append((group, train, test))
    return tuple(result)


def _average_ranks(values: FloatArray) -> FloatArray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1.0
        start = end
    return ranks


def spearman_correlation(first: ArrayLike, second: ArrayLike) -> float:
    """Return Spearman's rank correlation with deterministic average ranks."""

    left = np.asarray(first, dtype=float)
    right = np.asarray(second, dtype=float)
    if left.ndim != 1 or right.ndim != 1 or left.shape != right.shape or left.size < 2:
        raise ValueError("rank inputs must be matching one-dimensional arrays")
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError("rank inputs must be finite")
    left_ranks = _average_ranks(left)
    right_ranks = _average_ranks(right)
    if np.ptp(left_ranks) == 0.0 or np.ptp(right_ranks) == 0.0:
        return 0.0
    return float(np.corrcoef(left_ranks, right_ranks)[0, 1])


def out_of_fold_metrics(observed: ArrayLike, predicted: ArrayLike) -> OutOfFoldMetrics:
    """Return multiplicative diagnostics for positive observed/predicted data."""

    actual = np.asarray(observed, dtype=float)
    estimate = np.asarray(predicted, dtype=float)
    if actual.ndim != 1 or estimate.ndim != 1 or actual.shape != estimate.shape:
        raise ValueError("observed and predicted must be matching one-dimensional arrays")
    if actual.size < 2 or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(estimate)):
        raise ValueError("at least two finite predictions are required")
    if np.any(actual <= 0.0) or np.any(estimate <= 0.0):
        raise ValueError("observed and predicted values must be strictly positive")
    residual = np.log(estimate) - np.log(actual)
    factor = np.exp(np.abs(residual))
    return OutOfFoldMetrics(
        int(actual.size),
        float(np.sqrt(np.mean(residual * residual))),
        float(np.median(factor)),
        float(np.quantile(factor, 0.9)),
        float(np.max(factor)),
        float(np.mean(factor <= 2.0)),
        spearman_correlation(actual, estimate),
    )


def convergence_tail_diagnostics(
    orders: Iterable[int],
    changes: Iterable[float],
    applicable: Iterable[bool],
    *,
    tolerance: float = 1.0e-5,
    maximum_order: int = 21,
) -> ConvergenceTailDiagnostics:
    """Classify a convergence channel using two successive applicable changes."""

    order_values = np.asarray(tuple(orders), dtype=int)
    change_values = np.asarray(tuple(changes), dtype=float)
    applicable_values = np.asarray(tuple(applicable), dtype=bool)
    if (
        order_values.ndim != 1
        or order_values.size == 0
        or order_values.shape != change_values.shape
        or order_values.shape != applicable_values.shape
    ):
        raise ValueError("orders, changes and applicable must be non-empty matching vectors")
    if np.any(np.diff(order_values) <= 0) or not np.all(np.isfinite(change_values)):
        raise ValueError("orders must increase and changes must be finite")
    if not np.isfinite(tolerance) or tolerance <= 0.0 or maximum_order < order_values[-1]:
        raise ValueError("invalid tolerance or maximum order")
    confirmation: int | None = None
    for index in range(1, order_values.size):
        if (
            applicable_values[index - 1]
            and applicable_values[index]
            and change_values[index - 1] <= tolerance
            and change_values[index] <= tolerance
        ):
            confirmation = int(order_values[index])
            break
    valid = change_values[applicable_values]
    q_values = [
        current / previous
        for previous, current in zip(valid[:-1], valid[1:])
        if previous > 0.0 and current > 0.0 and np.isfinite(current / previous)
    ][-4:]
    tail = valid[-5:]
    directions = np.sign(np.diff(tail))
    directions = directions[directions != 0.0]
    alternations = int(np.sum(directions[1:] * directions[:-1] < 0.0))
    if confirmation is not None:
        classification = "directly_confirmed"
    elif int(order_values[-1]) == maximum_order:
        classification = "unconfirmed_at_21"
    else:
        classification = "not_applicable"
    return ConvergenceTailDiagnostics(
        confirmation,
        float(change_values[-1]),
        classification,
        float(np.median(q_values)) if q_values else 0.0,
        float(np.min(q_values)) if q_values else 0.0,
        float(np.max(q_values)) if q_values else 0.0,
        len(q_values),
        alternations >= 2,
    )
