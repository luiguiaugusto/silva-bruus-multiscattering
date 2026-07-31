"""Independent unit tests for the confirmatory T12.2 protocol."""

from dataclasses import replace

import numpy as np
import pytest

from acoustic_ms import (
    ConfirmatoryMetrics,
    LogoFoldFit,
    LogoPrediction,
    SafetyAudit,
    classify_logo_safety,
    confirmatory_metrics,
    evaluate_recalibration_gate,
    grouped_bootstrap_calibration,
    logo_power_law_predictions,
    power_law_threshold,
)


def _synthetic_logo():
    case_ids = [f"case_{index:02d}" for index in range(12)]
    groups = ["a"] * 4 + ["b"] * 4 + ["c"] * 4
    predictor = np.geomspace(0.002, 0.08, 12)
    observed = 7.5 * predictor**1.35
    return case_ids, groups, predictor, observed


def test_logo_recovers_exact_power_law_and_strict_holdout():
    case_ids, groups, predictor, observed = _synthetic_logo()
    fits, predictions = logo_power_law_predictions(case_ids, predictor, observed, groups)
    assert len(fits) == 3
    assert len(predictions) == len(case_ids)
    assert {item.case_id for item in predictions} == set(case_ids)
    for fit in fits:
        assert fit.prefactor == pytest.approx(7.5, rel=3e-14)
        assert fit.exponent == pytest.approx(1.35, rel=3e-14)
        assert fit.training_count == 8
        assert fit.test_count == 4
    for prediction in predictions:
        index = case_ids.index(prediction.case_id)
        assert prediction.held_out_group == groups[index]
        assert prediction.predicted == pytest.approx(observed[index], rel=5e-14)


def test_logo_is_invariant_to_input_row_order():
    case_ids, groups, predictor, observed = _synthetic_logo()
    baseline = logo_power_law_predictions(case_ids, predictor, observed, groups)
    order = np.array([8, 2, 10, 0, 6, 4, 11, 1, 9, 3, 7, 5])
    shuffled = logo_power_law_predictions(
        np.asarray(case_ids)[order], predictor[order], observed[order], np.asarray(groups)[order]
    )
    assert baseline == shuffled


def test_logo_rejects_duplicates_nonfinite_and_nonpositive_values():
    case_ids, groups, predictor, observed = _synthetic_logo()
    with pytest.raises(ValueError, match="unique"):
        logo_power_law_predictions([case_ids[0]] * len(case_ids), predictor, observed, groups)
    for invalid in (0.0, -1.0, np.nan, np.inf):
        changed = predictor.copy()
        changed[0] = invalid
        with pytest.raises(ValueError):
            logo_power_law_predictions(case_ids, changed, observed, groups)


def test_power_law_threshold_is_the_independent_analytic_inverse():
    threshold = power_law_threshold(0.05, 8.0, 1.5)
    assert 8.0 * threshold**1.5 == pytest.approx(0.05, rel=2e-15)
    for arguments in ((0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0), (np.inf, 1.0, 1.0)):
        with pytest.raises(ValueError):
            power_law_threshold(*arguments)


def test_confirmatory_metrics_match_manual_example():
    observed = np.array([1.0, 2.0, 4.0, 8.0])
    predicted = np.array([2.0, 1.0, 8.0, 4.0])
    metrics = confirmatory_metrics(observed, predicted)
    assert metrics.rmse_log == pytest.approx(np.log(2.0))
    assert metrics.mae_log == pytest.approx(np.log(2.0))
    assert metrics.median_absolute_log_ratio == pytest.approx(np.log(2.0))
    assert metrics.fraction_within_factor_2 == 1.0
    assert metrics.fraction_within_factor_1_5 == 0.0
    assert metrics.maximum_log_underestimation == pytest.approx(np.log(2.0))


def test_safety_classification_counts_false_safe_and_minimum_coverage():
    fits = (
        LogoFoldFit("a", 4, 2, 1.0, 1.0),
        LogoFoldFit("b", 4, 2, 1.0, 1.0),
    )
    predictions = (
        LogoPrediction("a1", "a", 0.02, 0.005, 0.005),
        LogoPrediction("a2", "a", 0.005, 0.009, 0.009),
        LogoPrediction("b1", "b", 0.004, 0.004, 0.004),
        LogoPrediction("b2", "b", 0.02, 0.02, 0.02),
    )
    classifications, audits = classify_logo_safety(predictions, fits, (0.01,))
    assert len(classifications) == 4
    audit = audits[0]
    assert audit.predicted_safe_count == 3
    assert audit.predicted_safe_group_count == 2
    assert audit.true_safe_count == 2
    assert audit.false_safe_count == 1
    assert audit.false_unsafe_count == 0
    assert audit.worst_false_safe_excess == pytest.approx(0.01)
    assert audit.coverage_sufficient


def test_group_bootstrap_is_reproducible_and_uses_requested_valid_count():
    _, groups, predictor, observed = _synthetic_logo()
    first = grouped_bootstrap_calibration(
        predictor, observed, groups, seed=42, valid_samples=250
    )
    second = grouped_bootstrap_calibration(
        predictor, observed, groups, seed=42, valid_samples=250
    )
    assert first == second
    assert first.valid_samples == 250
    assert first.attempts >= first.valid_samples
    assert first.prefactor_interval == pytest.approx((7.5, 7.5), rel=2e-13)
    assert first.exponent_interval == pytest.approx((1.35, 1.35), rel=2e-13)


def _passing_metrics() -> ConfirmatoryMetrics:
    return ConfirmatoryMetrics(28, 0.2, 0.15, 0.1, 0.95, 0.8, 0.95, 0.3)


def test_gate_positive_and_false_safe_negative_fixtures():
    candidate = _passing_metrics()
    baseline = replace(candidate, rmse_log=0.4, fraction_within_factor_2=0.7)
    fits = tuple(LogoFoldFit(str(index), 24, 4, 10.0, 1.4) for index in range(7))
    audits = tuple(SafetyAudit(value, 4, 2, 4, 0, 0, 0.0, True) for value in (0.01, 0.05, 0.10))
    criteria, decision = evaluate_recalibration_gate(
        candidate, baseline, fits, audits,
        predictions_finite_positive=True, integrity_passed=True,
    )
    assert len(criteria) == 10
    assert all(item.passed for item in criteria)
    assert decision == "GO_T13_WITH_RECALIBRATED_RHO1"
    failing = list(audits)
    failing[-1] = replace(failing[-1], false_safe_count=1, worst_false_safe_excess=0.02)
    criteria, decision = evaluate_recalibration_gate(
        candidate, baseline, fits, failing,
        predictions_finite_positive=True, integrity_passed=True,
    )
    assert not next(item for item in criteria if item.name == "zero_false_safe").passed
    assert decision == "NO_GO_T13_RHO1_NOT_QUANTITATIVE"
