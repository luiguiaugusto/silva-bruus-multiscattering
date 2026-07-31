"""Response-blind tests for the frozen T14 preregistration."""

from __future__ import annotations

import csv
from hashlib import sha256
import inspect
import subprocess
from pathlib import Path

import numpy as np
import pytest

from acoustic_ms import (
    EXPECTED_SCALE_OUT_CASE_IDS,
    ExternalPredictionMetrics,
    ExternalThresholdAudit,
    build_scale_out_cases,
    compact_cluster,
    evaluate_scale_out_gate,
    frozen_external_predictions,
    geometric_coupling_sum,
    irregular_scale_template,
    maximum_geometric_coupling,
    scale_out_template,
    triangular_compact_template,
)
from acoustic_ms.external_validation import LAMBDA_TARGETS, TOLERANCES


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
PHASE_A = (
    DATA / "t14_scale_manifest.csv",
    DATA / "t14_frozen_predictions.csv",
    DATA / "t14_frozen_protocol.csv",
    DATA / "t14_prior_artifact_hashes.csv",
)
EXPECTED_IDS = (
    "t14_n15_linear_level1", "t14_n15_linear_level2",
    "t14_n15_linear_level3", "t14_n15_linear_level4",
    "t14_n15_compact_level1", "t14_n15_compact_level2",
    "t14_n15_compact_level3", "t14_n15_compact_level4",
    "t14_n15_irregular_level1", "t14_n15_irregular_level2",
    "t14_n15_irregular_level3", "t14_n15_irregular_level4",
    "t14_n28_linear_level1", "t14_n28_linear_level2",
    "t14_n28_linear_level3", "t14_n28_linear_level4",
    "t14_n28_compact_level1", "t14_n28_compact_level2",
    "t14_n28_compact_level3", "t14_n28_compact_level4",
    "t14_n28_irregular_level1", "t14_n28_irregular_level2",
    "t14_n28_irregular_level3", "t14_n28_irregular_level4",
)


def _minimum_distance(positions: np.ndarray) -> float:
    distances = np.linalg.norm(
        positions[:, None, :] - positions[None, :, :], axis=2
    )
    return float(np.min(distances[np.triu_indices(len(positions), 1)]))


def _metrics(**changes: float) -> ExternalPredictionMetrics:
    values = {
        "point_count": 24, "rmse_log": 0.1, "mae_log": 0.08,
        "median_factor": 1.1, "p90_factor": 1.3, "maximum_factor": 1.5,
        "fraction_within_factor_2": 1.0,
        "fraction_within_factor_1_5": 0.9, "spearman": 0.99,
        "mean_log_bias": 0.0, "median_log_bias": 0.0,
    }
    values.update(changes)
    return ExternalPredictionMetrics(**values)


def _audits(false_safe: int = 0) -> tuple[ExternalThresholdAudit, ...]:
    result = []
    for tolerance, count in zip(TOLERANCES, (6, 12, 18)):
        result.append(ExternalThresholdAudit(
            model="M1", scope="all", tolerance=tolerance,
            eligible_count=24, predicted_safe_count=count,
            observed_safe_count=count - false_safe,
            false_safe_count=false_safe, false_unsafe_count=0,
            safe_precision=1.0 if not false_safe else 0.0,
            safe_coverage=1.0, worst_predicted_safe_error=tolerance / 2,
            false_safe_ids=() if not false_safe else ("synthetic",),
            false_unsafe_ids=(),
        ))
    return tuple(result)


def _gate_kwargs() -> dict[str, object]:
    return {
        "eligible_count": 24,
        "eligible_by_n": {15: 12, 28: 12},
        "eligible_by_family": {"linear": 8, "compact": 8, "irregular": 8},
        "eligible_by_level": {1: 6, 2: 6, 3: 6, 4: 6},
        "predicted_safe_eligible": {0.01: 6, 0.05: 12, 0.1: 18},
        "manifest_intact": True, "phase_a_integrity": True,
        "prior_integrity": True, "maximum_lmax": 13,
        "protocol_immutable": True, "resource_limit": False,
        "m1_global": _metrics(),
        "m1_by_n": {15: _metrics(point_count=12), 28: _metrics(point_count=12)},
        "m1_audits": _audits(),
    }


def test_exact_ids_and_distribution() -> None:
    cases = build_scale_out_cases()
    assert EXPECTED_SCALE_OUT_CASE_IDS == EXPECTED_IDS
    assert tuple(case.case_id for case in cases) == EXPECTED_IDS
    assert sum(case.particle_count == 15 for case in cases) == 12
    assert sum(case.particle_count == 28 for case in cases) == 12
    for family in ("linear", "compact", "irregular"):
        assert sum(case.family == family for case in cases) == 8
    for level in range(1, 5):
        assert sum(case.target_level == level for case in cases) == 6


@pytest.mark.parametrize("particle_count", (6, 10))
def test_triangular_template_reproduces_historical_compact(
    particle_count: int,
) -> None:
    assert np.allclose(
        triangular_compact_template(particle_count),
        compact_cluster(particle_count, 1.0), rtol=0.0, atol=3e-16,
    )


@pytest.mark.parametrize("particle_count", (15, 28))
@pytest.mark.parametrize("family", ("linear", "compact", "irregular"))
def test_templates_are_deterministic_centered_planar_and_normalized(
    particle_count: int, family: str,
) -> None:
    first = scale_out_template(particle_count, family)
    second = scale_out_template(particle_count, family)
    assert np.array_equal(first, second)
    assert first.shape == (particle_count, 3)
    assert np.all(np.isfinite(first))
    assert np.all(first[:, 2] == 0.0)
    assert np.allclose(np.mean(first, axis=0), 0.0, rtol=0.0, atol=2e-15)
    assert np.isclose(_minimum_distance(first), 1.0, rtol=2e-15, atol=2e-15)


def test_irregular_algorithm_is_deterministic_and_nontrivial() -> None:
    assert np.array_equal(irregular_scale_template(15), irregular_scale_template(15))
    assert not np.allclose(irregular_scale_template(15), triangular_compact_template(15))
    assert geometric_coupling_sum(irregular_scale_template(15)) > 0.0


@pytest.mark.parametrize("particle_count", (15, 28))
def test_all_exact_lambda_targets_and_nonoverlap(particle_count: int) -> None:
    for case in build_scale_out_cases():
        if case.particle_count != particle_count:
            continue
        assert case.minimum_distance > 2.0
        assert np.isclose(
            maximum_geometric_coupling(case.positions_xyz, 1.0, 0.8),
            LAMBDA_TARGETS[case.target_level - 1], rtol=5e-13, atol=5e-15,
        )


def test_blind_m1_predictions_and_exact_safe_counts() -> None:
    expected_points = (
        0.0014863937048093015, 0.0087007541172098232,
        0.027514200371470676, 0.10186146774330709,
    )
    expected_safe = (
        0.0038199876936037326, 0.02236067977499789,
        0.070710678118654738, 0.26178094805761298,
    )
    predictions = []
    for case in build_scale_out_cases():
        m1 = frozen_external_predictions(case.case_id, case.lambda_target, 0.01)[0]
        assert np.isclose(m1.point_prediction, expected_points[case.target_level - 1], rtol=3e-15, atol=3e-17)
        assert np.isclose(m1.conservative_prediction, expected_safe[case.target_level - 1], rtol=3e-15, atol=3e-17)
        predictions.append(m1)
    assert tuple(sum(getattr(item, field) for item in predictions) for field in ("safe_1pct", "safe_5pct", "safe_10pct")) == (6, 12, 18)


def test_gate_pass_fail_inconclusive_and_p3_cannot_intervene() -> None:
    _, decision, next_gate = evaluate_scale_out_gate(**_gate_kwargs())
    assert (decision, next_gate) == (
        "PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX",
        "GO_T15_SYNTHESIS_AND_MANUSCRIPT",
    )
    failing = _gate_kwargs()
    failing["m1_audits"] = _audits(false_safe=1)
    _, decision, next_gate = evaluate_scale_out_gate(**failing)
    assert (decision, next_gate) == (
        "FAIL_T14_SCALE_OUT_FROZEN_LAMBDA_MAX",
        "GO_T15_SYNTHESIS_WITH_SCALE_BREAKDOWN",
    )
    incomplete = _gate_kwargs()
    incomplete.update(eligible_count=19, m1_global=None, m1_by_n={}, m1_audits=())
    _, decision, next_gate = evaluate_scale_out_gate(**incomplete)
    assert (decision, next_gate) == (
        "INCONCLUSIVE_T14_INSUFFICIENT_MODEL_E_CONVERGENCE",
        "HOLD_T15_T14_INCONCLUSIVE",
    )
    assert "p3" not in inspect.signature(evaluate_scale_out_gate).parameters


def test_preregistration_is_blind_and_phase_a_is_immutable() -> None:
    script = ROOT / "scripts" / "preregister_t14_scale_out.py"
    source = script.read_text(encoding="utf-8")
    assert "solve_model_e_nodal" not in source
    for path in PHASE_A:
        expected = subprocess.run(
            ["git", "show", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT,
            check=True, capture_output=True,
        ).stdout
        assert sha256(path.read_bytes()).digest() == sha256(expected).digest()
    with (DATA / "t14_prior_artifact_hashes.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 81
    for row in rows:
        artifact = ROOT / row["path"]
        assert artifact.stat().st_size == int(row["size_bytes"])
        assert sha256(artifact.read_bytes()).hexdigest() == row["sha256"]

@pytest.mark.parametrize("particle_count", (0, 6, 10, 14, 29))
def test_scale_out_template_rejects_unregistered_sizes(particle_count: int) -> None:
    with pytest.raises(ValueError):
        scale_out_template(particle_count, "linear")
