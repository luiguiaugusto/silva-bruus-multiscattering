"""Unit and leakage tests for the T12.3 mechanistic criterion."""

from __future__ import annotations

import numpy as np
import pytest

from acoustic_ms import (
    MultiplicativeMetrics,
    audit_safety_thresholds,
    evaluate_mechanistic_gate,
    fit_mechanistic_power_law,
    nested_logo_predictions,
    predict_mechanistic_power_law,
)


def test_m1_and_m2_recover_exact_prespecified_power_laws():
    lam = np.array([0.01, 0.015, 0.03, 0.06, 0.11, 0.2])
    rho = np.array([0.009, 0.021, 0.025, 0.072, 0.09, 0.24])
    m1 = fit_mechanistic_power_law(lam, 2.5 * lam**1.3)
    assert m1.prefactor == pytest.approx(2.5, rel=1e-13)
    assert m1.alpha_lambda == pytest.approx(1.3, rel=1e-13)
    assert m1.alpha_rho == 0.0
    m2 = fit_mechanistic_power_law(lam, 3.2 * lam**0.8 * rho**0.35, rho)
    assert m2.prefactor == pytest.approx(3.2, rel=1e-13)
    assert m2.alpha_lambda == pytest.approx(0.8, rel=1e-13)
    assert m2.alpha_rho == pytest.approx(0.35, rel=1e-13)
    np.testing.assert_allclose(
        predict_mechanistic_power_law(m2, lam, rho),
        3.2 * lam**0.8 * rho**0.35,
        rtol=2e-14,
    )


@pytest.mark.parametrize(
    "lam, observed, rho",
    [
        ([0.1, 0.2], [0.1, 0.2], None),
        ([0.1, 0.2, np.nan], [0.1, 0.2, 0.3], None),
        ([0.1, 0.2, 0.3], [0.1, 0.0, 0.3], None),
        ([0.1, 0.2, 0.3], [0.1, 0.2, 0.3], [0.2, 0.3]),
        ([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], None),
    ],
)
def test_power_law_input_validation(lam, observed, rho):
    with pytest.raises(ValueError):
        fit_mechanistic_power_law(lam, observed, rho)


def _synthetic_logo(observed_shift: float = 1.0, permutation=None):
    groups = np.repeat(np.array(["g0", "g1", "g2", "g3"]), 3)
    case_ids = np.array([f"case_{index:02d}" for index in range(12)])
    lam = np.geomspace(0.01, 0.2, 12)
    rho = np.array([0.012, 0.014, 0.019, 0.021, 0.033, 0.039, 0.052, 0.071, 0.085, 0.11, 0.16, 0.25])
    observed = 4.0 * lam**1.2
    observed[groups == "g0"] *= observed_shift
    if permutation is None:
        permutation = np.arange(12)
    return nested_logo_predictions(
        case_ids[permutation], groups[permutation], lam[permutation], rho[permutation],
        observed[permutation], model="M1",
    )


def test_outer_group_never_enters_fit_inner_logo_or_safety_factor():
    folds_a, predictions_a = _synthetic_logo(1.0)
    folds_b, predictions_b = _synthetic_logo(1e12)
    fold_a = next(fold for fold in folds_a if fold.held_out_group == "g0")
    fold_b = next(fold for fold in folds_b if fold.held_out_group == "g0")
    assert fold_a.fit == fold_b.fit
    assert fold_a.safety_factor == fold_b.safety_factor
    assert fold_a.inner_prediction_count == fold_a.training_count == 9
    point_a = {row.case_id: row.point_prediction for row in predictions_a if row.held_out_group == "g0"}
    point_b = {row.case_id: row.point_prediction for row in predictions_b if row.held_out_group == "g0"}
    assert point_a == point_b


def test_nested_logo_is_invariant_to_input_row_order():
    folds_a, predictions_a = _synthetic_logo()
    folds_b, predictions_b = _synthetic_logo(permutation=np.array([8, 1, 10, 3, 6, 0, 11, 5, 2, 9, 4, 7]))
    assert folds_a == folds_b
    assert predictions_a == predictions_b
    assert len(predictions_a) == len({row.case_id for row in predictions_a}) == 12


def test_strict_safety_threshold_and_antivacuity_counts():
    ids = ["equal", "safe", "false_safe", "false_unsafe"]
    observed = [0.01, 0.005, 0.02, 0.005]
    predicted = [0.009, 0.004, 0.009, 0.02]
    audit = audit_safety_thresholds(ids, observed, predicted, model="M1", rule="conservative", tolerances=[0.01])[0]
    assert audit.predicted_safe_count == 3
    assert audit.observed_safe_count == 2
    assert audit.false_safe_ids == ("equal", "false_safe")
    assert audit.false_unsafe_ids == ("false_unsafe",)


def test_literal_gate_selects_m1_before_m2():
    ids = [f"c{i}" for i in range(12)]
    groups = np.repeat(["g0", "g1", "g2"], 4)
    lam = np.geomspace(0.001, 0.03, 12)
    rho = np.geomspace(0.0013, 0.04, 12) * (1 + 0.02 * np.sin(np.arange(12)))
    observed = 2 * lam
    m1_folds, m1_predictions = nested_logo_predictions(ids, groups, lam, rho, observed, model="M1")
    m2_folds, m2_predictions = nested_logo_predictions(ids, groups, lam, rho, observed, model="M2")
    m1_folds = m1_folds + m1_folds + m1_folds[:1]
    m2_folds = m2_folds + m2_folds + m2_folds[:1]
    m1_metrics = MultiplicativeMetrics(12, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0)
    m2_metrics = MultiplicativeMetrics(12, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0)
    m1_audits = audit_safety_thresholds(ids, observed, [row.safe_prediction for row in m1_predictions], model="M1", rule="conservative")
    m2_audits = audit_safety_thresholds(ids, observed, [row.safe_prediction for row in m2_predictions], model="M2", rule="conservative")
    # Supply the pre-registered anti-vacuity counts explicitly to isolate gate logic.
    m1_audits = tuple(type(item)(item.model, item.rule, item.tolerance, count, max(count, item.observed_safe_count), 0, 0, 1.0, 1.0, (), ()) for item, count in zip(m1_audits, (3, 8, 12)))
    m2_audits = tuple(type(item)(item.model, item.rule, item.tolerance, count, max(count, item.observed_safe_count), 0, 0, 1.0, 1.0, (), ()) for item, count in zip(m2_audits, (3, 8, 12)))
    _, decision, m1_pass, _ = evaluate_mechanistic_gate(
        m1_metrics, m2_metrics, m1_audits, m2_audits, m1_folds, m2_folds,
        fit_mechanistic_power_law(lam, observed),
        fit_mechanistic_power_law(lam, observed * rho**0.01, rho),
        integrity_passed=True, m2_unstable_collinearity=False,
    )
    assert m1_pass
    assert decision == "GO_T13_VALIDATE_LAMBDA_MAX"
