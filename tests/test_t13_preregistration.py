"""Response-blind selection and frozen-protocol tests for T13 phase A."""

from __future__ import annotations

import csv
from pathlib import Path
import random

import numpy as np
import pytest

from acoustic_ms import (
    EXPECTED_CASE_IDS,
    ExternalPredictionMetrics,
    ExternalThresholdAudit,
    audit_external_threshold,
    evaluate_external_validation_gate,
    frozen_external_predictions,
    minimum_two_step_confirmation,
    select_external_validation_cases,
    successive_change,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"


def _t08_rows() -> list[dict[str, str]]:
    with (DATA / "t08_cases.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_exact_response_blind_selection_and_distribution():
    selected = select_external_validation_cases(_t08_rows())
    assert tuple(case.case_id for case in selected) == EXPECTED_CASE_IDS
    assert len(selected) == 24
    assert [case.target_level for case in selected] == [1, 2, 3, 4] * 6
    assert {case.particle_count for case in selected} == {6, 10}
    for stratum in {case.stratum for case in selected}:
        assert sum(case.stratum == stratum for case in selected) == 4


def test_selection_ignores_responses_and_is_row_order_invariant():
    rows = _t08_rows()
    expected = tuple(case.case_id for case in select_external_validation_cases(rows))
    protected = {
        "case_id", "split", "particle_count", "family", "f1",
        "distance_ratio", "lambda_max", "rho_l1", "reference_lmax",
        "total_converged",
    }
    mutated = []
    for row in rows:
        copy = dict(row)
        for field in set(copy) - protected:
            copy[field] = "response_column_replaced_by_test"
        mutated.append(copy)
    random.Random(1300).shuffle(mutated)
    assert tuple(case.case_id for case in select_external_validation_cases(mutated)) == expected


def test_frozen_m1_p3_formulas_and_blind_safe_counts():
    selected = select_external_validation_cases(_t08_rows())
    predictions = [
        prediction
        for case in selected
        for prediction in frozen_external_predictions(case.case_id, case.lambda_max, case.rho_l1)
    ]
    assert len(predictions) == 48
    m1 = [prediction for prediction in predictions if prediction.model == "M1"]
    assert [sum(getattr(prediction, field) for prediction in m1) for field in ("safe_1pct", "safe_5pct", "safe_10pct")] == [6, 12, 18]
    first_m1, first_p3 = frozen_external_predictions("case", 0.01, 0.005)
    assert first_m1.point_prediction == pytest.approx(4.4964255121671126 * 0.01**1.3883601043764593)
    assert first_m1.safety_factor == 2.5699703122019222
    assert first_p3.point_prediction == pytest.approx(14.73950709797405 * 0.005**1.4226504975598322)
    assert first_p3.safety_factor == 2.0464420079866286


def test_strict_threshold_equality_is_unsafe_and_false_safe_is_correct():
    audit = audit_external_threshold(
        ["equal", "safe", "false_safe", "false_unsafe"],
        [0.01, 0.005, 0.02, 0.005],
        [0.009, 0.004, 0.009, 0.02],
        model="M1", scope="all", tolerance=0.01,
    )
    assert audit.predicted_safe_count == 3
    assert audit.observed_safe_count == 2
    assert audit.false_safe_ids == ("equal", "false_safe")
    assert audit.false_unsafe_ids == ("false_unsafe",)


def test_two_successive_changes_and_scale_aware_zero_handling():
    orders = [2, 3, 4, 5]
    assert minimum_two_step_confirmation([0, 1e-4, 1e-6, 2e-6], [False, True, True, True], orders) == 5
    assert minimum_two_step_confirmation([0, 1e-4, 1e-6, 2e-5], [False, True, True, True], orders) == 0
    change, applicable, absolute = successive_change(np.zeros((2, 3)), np.zeros((2, 3)))
    assert change == absolute == 0.0
    assert not applicable


def _metric() -> ExternalPredictionMetrics:
    return ExternalPredictionMetrics(24, 0.2, 0.1, 1.1, 1.3, 1.5, 0.9, 0.85, 0.95, 0.0, 0.0)


def _audits(false_safe: int = 0) -> tuple[ExternalThresholdAudit, ...]:
    return tuple(
        ExternalThresholdAudit(
            "M1", "all", tolerance, 24, count, count,
            false_safe if tolerance == 0.10 else 0, 0, 1.0, 1.0,
            tolerance / 2, ("bad",) if false_safe and tolerance == 0.10 else (), (),
        )
        for tolerance, count in zip((0.01, 0.05, 0.10), (3, 6, 9))
    )


def test_gate_precedence_and_literal_pass_fail_outcomes():
    common = dict(
        diagnostics_all_passed=True, manifest_intact=True, integrity_passed=True,
        m1_global=_metric(), m1_n6=_metric(), m1_n10=_metric(),
        m1_audits=_audits(),
        predicted_safe_n_by_tolerance={0.05: (3, 3), 0.10: (4, 5)},
    )
    criteria, decision, t14 = evaluate_external_validation_gate(
        eligible_count=19, eligible_n6=10, eligible_n10=9, **common
    )
    assert decision == "INCONCLUSIVE_T13_INSUFFICIENT_MODEL_E_CONVERGENCE"
    assert t14 == "HOLD_T14_MODEL_E_CONVERGENCE"
    assert all(item.stage == "sufficiency" for item in criteria)
    _, decision, t14 = evaluate_external_validation_gate(
        eligible_count=24, eligible_n6=12, eligible_n10=12, **common
    )
    assert decision == "PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX"
    assert t14 == "GO_T14_SCALE_OUT_WITH_FROZEN_LAMBDA_MAX"
    failed = dict(common)
    failed["m1_audits"] = _audits(false_safe=1)
    _, decision, t14 = evaluate_external_validation_gate(
        eligible_count=24, eligible_n6=12, eligible_n10=12, **failed
    )
    assert decision == "FAIL_T13_EXTERNAL_VALIDATION_LAMBDA_MAX"
    assert t14 == "NO_GO_T14_LAMBDA_MAX_NOT_TRANSFERABLE"


def test_p3_cannot_change_the_m1_gate_api():
    parameters = evaluate_external_validation_gate.__annotations__
    assert all("p3" not in name.lower() for name in parameters)
