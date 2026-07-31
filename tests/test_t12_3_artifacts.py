"""Integration and integrity tests for deterministic T12.3 products."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURE = ROOT / "results" / "figures" / "t12_3_mechanistic_validity.png"
SCRIPT = ROOT / "scripts" / "analyze_t12_3_mechanistic_validity.py"
ARTIFACTS = (
    DATA / "t12_3_oof_predictions.csv",
    DATA / "t12_3_logo_coefficients.csv",
    DATA / "t12_3_nested_safety_factors.csv",
    DATA / "t12_3_threshold_audit.csv",
    DATA / "t12_3_metrics.csv",
    DATA / "t12_3_group_bootstrap.csv",
    DATA / "t12_3_case_influence.csv",
    DATA / "t12_3_gate.csv",
    FIGURE,
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_exact_frozen_ids_and_no_external_holdout_are_used():
    source = _read(DATA / "t12_1_resolved_comparison.csv")
    rows = _read(DATA / "t12_3_oof_predictions.csv")
    assert len(source) == 28
    assert len(rows) == 112
    assert [row["case_id"] for row in rows[::4]] == [row["case_id"] for row in source]
    assert {int(row["particle_count"]) for row in rows} == {2, 3, 4}
    assert not ({6, 10} & {int(row["particle_count"]) for row in rows})
    for case_id in {row["case_id"] for row in source}:
        assert sum(row["case_id"] == case_id for row in rows) == 4
        assert {row["model"] for row in rows if row["case_id"] == case_id} == {"P0", "P3", "M1", "M2"}


def test_outer_logo_and_nested_margin_metadata_are_complete():
    coefficients = _read(DATA / "t12_3_logo_coefficients.csv")
    outer = [row for row in coefficients if row["scope"] == "outer_fold"]
    assert len(outer) == 14
    assert all(int(row["training_count"]) == 24 for row in outer)
    assert all(int(row["test_count"]) == 4 for row in outer)
    margins = _read(DATA / "t12_3_nested_safety_factors.csv")
    assert len(margins) == 28
    assert all(int(row["inner_prediction_count"]) == 24 for row in margins)
    assert all(row["valid"] == "true" for row in margins)
    assert all(row["outer_group_excluded_from_fit_and_margin"] == "true" for row in margins)


def test_m1_regression_and_conservative_gate_are_exact():
    metrics = _read(DATA / "t12_3_metrics.csv")
    m1 = next(row for row in metrics if row["model"] == "M1" and row["scope_type"] == "global")
    assert float(m1["rmse_log"]) == pytest.approx(0.6293890920247307, rel=5e-14)
    assert float(m1["fraction_within_factor_2"]) == pytest.approx(25 / 28)
    audit = [row for row in _read(DATA / "t12_3_threshold_audit.csv") if row["model"] == "M1" and row["rule"] == "conservative"]
    assert [int(row["predicted_safe_count"]) for row in audit] == [4, 9, 14]
    assert [int(row["false_safe_count"]) for row in audit] == [0, 0, 0]
    gate = _read(DATA / "t12_3_gate.csv")
    assert {row["final_decision"] for row in gate} == {"GO_T13_VALIDATE_LAMBDA_MAX"}
    assert all(row["candidate_pass"] == "true" for row in gate if row["candidate"] == "M1")


def test_m2_collinearity_is_reported_without_regularization():
    fits = _read(DATA / "t12_3_logo_coefficients.csv")
    m2 = [row for row in fits if row["model"] == "M2" and row["scope"] == "outer_fold"]
    assert len(m2) == 7
    assert sum(float(row["alpha_rho"]) > 0 for row in m2) == 4
    assert all(float(row["standardized_condition_number"]) < 1e3 for row in m2)
    gate = _read(DATA / "t12_3_gate.csv")
    assert {row["m2_collinearity_status"] for row in gate} == {"UNSTABLE_COLLINEARITY"}


def test_all_eight_csvs_are_finite_ordered_and_newline_terminated():
    expected = {
        "t12_3_oof_predictions.csv": 112,
        "t12_3_logo_coefficients.csv": 16,
        "t12_3_nested_safety_factors.csv": 28,
        "t12_3_threshold_audit.csv": 24,
        "t12_3_metrics.csv": 96,
        "t12_3_group_bootstrap.csv": 18,
        "t12_3_case_influence.csv": 56,
        "t12_3_gate.csv": 19,
    }
    for name, count in expected.items():
        path = DATA / name
        text = path.read_text(encoding="utf-8")
        assert text.endswith("\n")
        assert len(_read(path)) == count
        assert "nan" not in text.lower()
        assert "inf" not in text.lower()
        header = text.splitlines()[0].split(",")
        assert len(header) == len(set(header))


def test_all_pre_t12_3_result_artifacts_match_git_head():
    paths = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", "HEAD", "results/data", "results/figures"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    protected = [path for path in paths if "/t12_3_" not in path]
    assert len(protected) == 61
    for relative in protected:
        expected = subprocess.run(
            ["git", "show", f"HEAD:{relative}"], cwd=ROOT, check=True,
            capture_output=True,
        ).stdout
        assert hashlib.sha256((ROOT / relative).read_bytes()).digest() == hashlib.sha256(expected).digest()


def test_t12_3_regeneration_is_byte_deterministic():
    before = {path.name: _hash(path) for path in ARTIFACTS}
    subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, capture_output=True, text=True)
    assert {path.name: _hash(path) for path in ARTIFACTS} == before
