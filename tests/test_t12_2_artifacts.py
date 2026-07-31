"""Integrity tests for the deterministic T12.2 calibration artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURE = ROOT / "results" / "figures" / "t12_2_rho1_recalibration.png"
SCRIPT = ROOT / "scripts" / "analyze_t12_2_rho1_recalibration.py"
ARTIFACTS = (
    DATA / "t12_2_logo_predictions.csv",
    DATA / "t12_2_logo_fits.csv",
    DATA / "t12_2_metrics.csv",
    DATA / "t12_2_safety_audit.csv",
    DATA / "t12_2_final_calibration.csv",
    DATA / "t12_2_gate.csv",
    FIGURE,
)
PROTECTED_HASHES = {
    DATA / "t12_sentinel_manifest.csv": "a46cf99bee4802ce42e29d5a9970c9fb8da7ae63940fa772ee2e9dc5f77befe1",
    DATA / "t12_model_e_convergence.csv": "a1cf482541c5d95fc5145d8a01a0e69c39fae737070094e0e03071175aaf8524",
    DATA / "t12_model_comparison.csv": "3fd672da68099a264c36497ca7b6ee548f5e6430a4d9ea0f493b1f06fdd5cf91",
    DATA / "t12_threshold_audit.csv": "a33544c8040693f6be7607c8ecce28a33c9bfd0ab58034dd98bad09b3e88a516",
    DATA / "t12_1_extended_convergence.csv": "d41a956e9c58e5d49ab06f94b7574f8c9f987610223f9692e2c1b67297019e23",
    DATA / "t12_1_convergence_summary.csv": "89addfd05f6d1c33160bd9bc1b3cbb6ff05f93a99dc1ce80c4c7d96a2184cbf0",
    DATA / "t12_1_resolved_comparison.csv": "5097cd7014bac635e09179e5bd4f49a0308dc4f6f02eb7bd76d60f18c2e89f39",
    DATA / "t12_1_mechanism_diagnostics.csv": "d3a968f322feddb7da9fd9f3fb470564c691b4b4db49ddca65248627e30e4334",
    DATA / "t12_1_predictor_diagnostics.csv": "d63c16d216a6235a41d44346e8c1f981e2ffd1ba3977c728713c8c57bd71842f",
    DATA / "t12_1_out_of_fold_predictions.csv": "7686c091a0323d011e79e50ebfe9dc096dc10240ca21c60b0749c9af488ac4d6",
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_t12_2_manifest_is_exact_and_excludes_external_holdout():
    rows = _read(DATA / "t12_2_logo_predictions.csv")
    source = _read(DATA / "t12_1_resolved_comparison.csv")
    assert len(rows) == 28
    assert [int(row["sentinel_order"]) for row in rows] == list(range(1, 29))
    assert [row["case_id"] for row in rows] == [row["case_id"] for row in source]
    assert {int(row["particle_count"]) for row in rows} == {2, 3, 4}
    assert {row["held_out_group"] for row in rows} == {
        "n2_pair", "n3_compact", "n3_irregular", "n3_linear",
        "n4_compact", "n4_irregular", "n4_linear",
    }


def test_all_total_references_and_exact_two_internal_statuses_are_preserved():
    rows = _read(DATA / "t12_2_logo_predictions.csv")
    assert all(row["total_confirmed"] == "true" for row in rows)
    assert all(row["interaction_confirmed"] == "true" for row in rows)
    assert {
        row["case_id"] for row in rows
        if row["scattered_scattered_confirmed"] != "true"
    } == {"n2_pair_f1.0_d2.1", "n3_irregular_f1.0_d2.1"}


def test_logo_folds_cover_each_group_once_without_leakage_metadata():
    fits = _read(DATA / "t12_2_logo_fits.csv")
    predictions = _read(DATA / "t12_2_logo_predictions.csv")
    assert len(fits) == 7
    assert all(int(row["training_count"]) == 24 for row in fits)
    assert all(int(row["test_count"]) == 4 for row in fits)
    assert all(row["coefficients_positive"] == "true" for row in fits)
    assert all(sum(item["held_out_group"] == row["held_out_group"] for item in predictions) == 4 for row in fits)


def test_safety_audit_identifies_the_single_ten_percent_false_safe():
    audit = _read(DATA / "t12_2_safety_audit.csv")
    assert len(audit) == 3
    assert [int(row["predicted_safe_count"]) for row in audit] == [7, 14, 20]
    assert [int(row["predicted_safe_group_count"]) for row in audit] == [7, 7, 7]
    assert [int(row["false_safe_count"]) for row in audit] == [0, 0, 1]
    assert all(row["coverage_sufficient"] == "true" for row in audit)
    predictions = _read(DATA / "t12_2_logo_predictions.csv")
    false_safe = [row["case_id"] for row in predictions if row["false_safe_10pct"] == "true"]
    assert false_safe == ["n2_pair_f0.8_d2.5"]


def test_gate_has_ten_criteria_and_reports_the_preregistered_no_go():
    rows = _read(DATA / "t12_2_gate.csv")
    criteria = [row for row in rows if row["record_type"] == "criterion"]
    assert len(criteria) == 10
    failed = [row["criterion"] for row in criteria if row["passed"] != "true"]
    assert failed == ["zero_false_safe"]
    decision = next(row for row in rows if row["record_type"] == "decision")
    assert decision["decision"] == "NO_GO_T13_RHO1_NOT_QUANTITATIVE"


def test_final_calibration_uses_group_bootstrap_metadata():
    row = _read(DATA / "t12_2_final_calibration.csv")[0]
    assert int(row["point_count"]) == 28
    assert int(row["group_count"]) == 7
    assert int(row["bootstrap_seed"]) == 1202
    assert int(row["bootstrap_valid_samples"]) == 10_000
    assert int(row["bootstrap_attempts"]) >= 10_000
    assert float(row["prefactor"]) > 0.0
    assert float(row["exponent"]) > 0.0
    assert float(row["epsilon_floor"]) == 0.0


def test_all_t12_2_csvs_are_finite_ordered_and_newline_terminated():
    expected_counts = {
        "t12_2_logo_predictions.csv": 28,
        "t12_2_logo_fits.csv": 7,
        "t12_2_metrics.csv": 16,
        "t12_2_safety_audit.csv": 3,
        "t12_2_final_calibration.csv": 1,
        "t12_2_gate.csv": 11,
    }
    for name, count in expected_counts.items():
        path = DATA / name
        text = path.read_text(encoding="utf-8")
        assert len(_read(path)) == count
        assert text.endswith("\n")
        assert "nan" not in text.lower()
        assert "inf" not in text.lower()
        header = text.splitlines()[0].split(",")
        assert len(header) == len(set(header))


def test_all_protected_t12_and_t12_1_hashes_are_unchanged():
    for path, expected in PROTECTED_HASHES.items():
        assert _hash(path) == expected


def test_t12_2_regeneration_is_byte_deterministic():
    before = {path.name: _hash(path) for path in ARTIFACTS}
    subprocess.run(
        [sys.executable, str(SCRIPT)], cwd=ROOT, check=True,
        capture_output=True, text=True,
    )
    assert {path.name: _hash(path) for path in ARTIFACTS} == before
