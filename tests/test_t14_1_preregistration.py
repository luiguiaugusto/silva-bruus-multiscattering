"""Response-blind tests for the frozen T14.1 preregistration."""

from __future__ import annotations

import csv
from hashlib import sha256
import inspect
from pathlib import Path
import subprocess

import numpy as np
import pytest

from acoustic_ms import (
    EXPECTED_LARGE_N_CASE_IDS,
    ExternalPredictionMetrics,
    ExternalThresholdAudit,
    build_large_n_cases,
    classify_large_n_trend,
    evaluate_large_n_gate,
    frozen_external_predictions,
    large_n_template,
    local_coupling_statistics,
    local_geometric_coupling,
    maximum_geometric_coupling,
    scale_out_template,
    triangular_compact_template,
)
from acoustic_ms.external_validation import LAMBDA_TARGETS, TOLERANCES


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
PHASE_A = (
    DATA / "t14_1_large_n_manifest.csv",
    DATA / "t14_1_local_coupling.csv",
    DATA / "t14_1_frozen_predictions.csv",
    DATA / "t14_1_frozen_protocol.csv",
    DATA / "t14_1_prior_artifact_hashes.csv",
)


def _minimum_distance(positions: np.ndarray) -> float:
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)
    return float(np.min(distances[np.triu_indices(len(positions), 1)]))


def _metrics(**changes: float) -> ExternalPredictionMetrics:
    values = {
        "point_count": 24, "rmse_log": 0.1, "mae_log": 0.08,
        "median_factor": 1.1, "p90_factor": 1.3, "maximum_factor": 1.5,
        "fraction_within_factor_2": 1.0, "fraction_within_factor_1_5": 0.9,
        "spearman": 0.99, "mean_log_bias": 0.0, "median_log_bias": 0.0,
    }
    values.update(changes)
    return ExternalPredictionMetrics(**values)


def _audits(false_safe: int = 0) -> tuple[ExternalThresholdAudit, ...]:
    return tuple(
        ExternalThresholdAudit(
            model="M1", scope="all", tolerance=tolerance,
            eligible_count=24, predicted_safe_count=count,
            observed_safe_count=count - false_safe,
            false_safe_count=false_safe, false_unsafe_count=0,
            safe_precision=1.0 if not false_safe else 0.0, safe_coverage=1.0,
            worst_predicted_safe_error=tolerance / 2,
            false_safe_ids=() if not false_safe else ("synthetic",),
            false_unsafe_ids=(),
        )
        for tolerance, count in zip(TOLERANCES, (6, 12, 18))
    )


def _gate_kwargs() -> dict[str, object]:
    return {
        "eligible_count": 24, "eligible_by_n": {45: 12, 105: 12},
        "eligible_by_family": {"linear": 8, "compact": 8, "irregular": 8},
        "eligible_by_level": {1: 6, 2: 6, 3: 6, 4: 6},
        "predicted_safe_eligible": {0.01: 6, 0.05: 12, 0.1: 18},
        "manifest_intact": True, "phase_a_integrity": True,
        "prior_integrity": True, "maximum_lmax": 13,
        "protocol_immutable": True, "resource_limit": False,
        "m1_global": _metrics(),
        "m1_by_n": {45: _metrics(point_count=12), 105: _metrics(point_count=12)},
        "m1_audits": _audits(),
    }


def test_exact_ids_and_strata() -> None:
    cases = build_large_n_cases()
    assert tuple(case.case_id for case in cases) == EXPECTED_LARGE_N_CASE_IDS
    assert len(cases) == 24
    assert [case.scale_order for case in cases] == list(range(1, 25))
    assert sum(case.particle_count == 45 for case in cases) == 12
    assert sum(case.particle_count == 105 for case in cases) == 12
    for family in ("linear", "compact", "irregular"):
        assert sum(case.family == family for case in cases) == 8
    for level in range(1, 5):
        assert sum(case.target_level == level for case in cases) == 6


@pytest.mark.parametrize(("particle_count", "rows"), ((45, 9), (105, 14)))
def test_triangular_identity_and_template_properties(particle_count: int, rows: int) -> None:
    assert rows * (rows + 1) // 2 == particle_count
    for family in ("linear", "compact", "irregular"):
        first = large_n_template(particle_count, family)
        second = large_n_template(particle_count, family)
        assert np.array_equal(first, second)
        assert first.shape == (particle_count, 3)
        assert np.all(np.isfinite(first))
        assert np.all(first[:, 2] == 0.0)
        assert np.allclose(np.mean(first, axis=0), 0.0, rtol=0.0, atol=2e-14)
        assert np.isclose(_minimum_distance(first), 1.0, rtol=3e-14, atol=3e-14)


@pytest.mark.parametrize("particle_count", (15, 28))
def test_large_n_algorithms_reproduce_t14_templates(particle_count: int) -> None:
    assert np.array_equal(triangular_compact_template(particle_count), scale_out_template(particle_count, "compact"))
    from acoustic_ms import irregular_scale_template
    assert np.array_equal(irregular_scale_template(particle_count), scale_out_template(particle_count, "irregular"))


def test_exact_lambda_targets_local_vectors_and_hashes() -> None:
    for case in build_large_n_cases():
        assert case.minimum_distance > 2.0
        assert np.isclose(maximum_geometric_coupling(case.positions_xyz, 1.0, 0.8), LAMBDA_TARGETS[case.target_level - 1], rtol=5e-13, atol=5e-15)
        local = local_geometric_coupling(case.positions_xyz, 1.0, 0.8)
        assert np.array_equal(local, case.local_coupling)
        first = local_coupling_statistics(local)
        second = local_coupling_statistics(local.copy())
        assert first == second
        assert np.isclose(first["lambda_max"], case.lambda_target, rtol=5e-13, atol=5e-15)
        assert len(str(first["local_coupling_sha256"])) == 64


def test_local_coupling_invariant_to_rigid_motion_and_permutation() -> None:
    positions = build_large_n_cases()[4].positions_xyz
    angle = 0.37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
    order = np.arange(len(positions))[::-1]
    reference = local_geometric_coupling(positions, 1.0, 0.8)
    moved = (positions @ rotation.T + np.array([2.1, -4.2, 0.0]))[order]
    assert np.allclose(local_geometric_coupling(moved, 1.0, 0.8), reference[order], rtol=2e-14, atol=2e-16)


def test_blind_m1_predictions_and_safe_counts() -> None:
    predictions = [frozen_external_predictions(case.case_id, case.lambda_target, 0.01)[0] for case in build_large_n_cases()]
    assert tuple(sum(getattr(item, field) for item in predictions) for field in ("safe_1pct", "safe_5pct", "safe_10pct")) == (6, 12, 18)


def test_gate_pass_fail_resource_and_p3_cannot_intervene() -> None:
    _, decision, next_gate = evaluate_large_n_gate(**_gate_kwargs())
    assert (decision, next_gate) == ("PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX", "GO_T15_SYNTHESIS_AND_MANUSCRIPT")
    failing = _gate_kwargs()
    failing["m1_audits"] = _audits(false_safe=1)
    assert evaluate_large_n_gate(**failing)[1:] == ("FAIL_T14_1_LARGE_N_FROZEN_LAMBDA_MAX", "GO_T15_SYNTHESIS_WITH_LARGE_N_BREAKDOWN")
    incomplete = _gate_kwargs()
    incomplete.update(eligible_count=19, m1_global=None, m1_by_n={}, m1_audits=())
    assert evaluate_large_n_gate(**incomplete)[1:] == ("INCONCLUSIVE_T14_1_INSUFFICIENT_MODEL_E_CONVERGENCE", "HOLD_T15_T14_1_INCONCLUSIVE")
    incomplete["resource_limit"] = True
    assert evaluate_large_n_gate(**incomplete)[1:] == ("INCONCLUSIVE_T14_1_RESOURCE_LIMIT", "HOLD_T15_T14_1_INCONCLUSIVE")
    assert "p3" not in inspect.signature(evaluate_large_n_gate).parameters


@pytest.mark.parametrize(
    ("values", "expected"),
    (([1.0] * 12, "NO_SYSTEMATIC_DETERIORATION"), ([1.3] * 9 + [1.0] * 3, "SYSTEMATIC_DETERIORATION"), ([1.2] * 12, "MIXED_LARGE_N_TREND"), ([1.0] * 9, "INCONCLUSIVE_LARGE_N_TREND")),
)
def test_frozen_large_n_trend_classification(values: list[float], expected: str) -> None:
    assert classify_large_n_trend(values) == expected


def test_preregistration_is_blind_idempotent_and_preserves_prior_hashes() -> None:
    source = (ROOT / "scripts" / "preregister_t14_1_large_n.py").read_text(encoding="utf-8")
    assert "solve_model_e_nodal" not in source
    assert "model_e_response_consulted" in source
    for path in PHASE_A:
        assert path.is_file() and path.stat().st_size > 0
        tracked = subprocess.run(["git", "cat-file", "-e", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT).returncode == 0
        if tracked:
            expected = subprocess.run(["git", "show", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT, check=True, capture_output=True).stdout
            assert sha256(path.read_bytes()).digest() == sha256(expected).digest()
    with (DATA / "t14_1_prior_artifact_hashes.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 95
    for row in rows:
        artifact = ROOT / row["path"]
        assert artifact.stat().st_size == int(row["size_bytes"])
        assert sha256(artifact.read_bytes()).hexdigest() == row["sha256"]


@pytest.mark.parametrize("particle_count", (0, 15, 28, 44, 106))
def test_large_n_template_rejects_unregistered_sizes(particle_count: int) -> None:
    with pytest.raises(ValueError):
        large_n_template(particle_count, "linear")


def test_local_coupling_validation() -> None:
    with pytest.raises(ValueError):
        local_geometric_coupling(np.zeros((1, 3)), 1.0, 0.8)
    with pytest.raises(ValueError):
        local_geometric_coupling(np.zeros((2, 3)), 1.0, 0.8)
    with pytest.raises(ValueError):
        local_geometric_coupling(np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0]]), 1.0, 0.8)
