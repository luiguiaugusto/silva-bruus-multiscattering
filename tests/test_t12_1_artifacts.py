"""Integrity and determinism tests for the T12.1 published artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPT = ROOT / "scripts" / "analyze_t12_1_rho1_failure.py"
EXTENDED_CASES = (
    "n2_pair_f1.0_d2.1",
    "n3_compact_f0.8_d2.1",
    "n3_irregular_f1.0_d2.1",
    "n3_linear_f1.0_d2.1",
    "n4_irregular_f0.8_d2.1",
    "n4_linear_f0.8_d2.1",
    "n3_compact_f0.1_d2.1",
    "n4_compact_f0.1_d2.1",
    "n4_irregular_f0.1_d2.1",
    "n4_linear_f0.1_d2.1",
)
DERIVED = (
    DATA / "t12_1_convergence_summary.csv",
    DATA / "t12_1_resolved_comparison.csv",
    DATA / "t12_1_mechanism_diagnostics.csv",
    DATA / "t12_1_predictor_diagnostics.csv",
    DATA / "t12_1_out_of_fold_predictions.csv",
    ROOT / "results" / "figures" / "t12_1_rho1_failure_diagnostics.png",
)
T12_HASHES = {
    "t12_sentinel_manifest.csv": "a46cf99bee4802ce42e29d5a9970c9fb8da7ae63940fa772ee2e9dc5f77befe1",
    "t12_model_e_convergence.csv": "a1cf482541c5d95fc5145d8a01a0e69c39fae737070094e0e03071175aaf8524",
    "t12_model_comparison.csv": "3fd672da68099a264c36497ca7b6ee548f5e6430a4d9ea0f493b1f06fdd5cf91",
    "t12_threshold_audit.csv": "a33544c8040693f6be7607c8ecce28a33c9bfd0ab58034dd98bad09b3e88a516",
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_extended_campaign_has_exact_cases_provenance_and_order():
    rows = _read(DATA / "t12_1_extended_convergence.csv")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
    assert tuple(grouped) == EXTENDED_CASES
    for case_rows in grouped.values():
        orders = [int(row["lmax"]) for row in case_rows]
        assert orders[:12] == list(range(2, 14))
        assert orders == list(range(2, orders[-1] + 1))
        assert all(row["source"] == "t12" for row in case_rows[:12])
        assert all(row["source"] == "t12_1" for row in case_rows[12:])


def test_resolved_comparison_is_the_frozen_28_case_manifest_without_n6_or_n10():
    rows = _read(DATA / "t12_1_resolved_comparison.csv")
    assert len(rows) == 28
    assert [int(row["sentinel_order"]) for row in rows] == list(range(1, 29))
    assert {int(row["particle_count"]) for row in rows} == {2, 3, 4}
    assert len({row["case_id"] for row in rows}) == 28


def test_convergence_summary_has_four_channels_per_extended_case():
    rows = _read(DATA / "t12_1_convergence_summary.csv")
    assert len(rows) == 40
    assert {
        (row["case_id"], row["channel"]) for row in rows
    } == {
        (case_id, channel)
        for case_id in EXTENDED_CASES
        for channel in ("total", "interaction", "external_scattered", "scattered_scattered")
    }
    assert all(row["classification"] != "not_applicable" for row in rows)


def test_interaction_channels_are_confirmed_and_two_ss_channels_reach_l21_unconfirmed():
    rows = _read(DATA / "t12_1_convergence_summary.csv")
    interaction = [row for row in rows if row["channel"] == "interaction"]
    unconfirmed = [
        row for row in rows if row["classification"] == "unconfirmed_at_21"
    ]
    assert all(row["classification"] == "directly_confirmed" for row in interaction)
    assert {
        (row["case_id"], row["channel"]) for row in unconfirmed
    } == {
        ("n2_pair_f1.0_d2.1", "scattered_scattered"),
        ("n3_irregular_f1.0_d2.1", "scattered_scattered"),
    }


def test_oof_predictions_have_exact_coverage_without_leakage():
    comparison = _read(DATA / "t12_1_resolved_comparison.csv")
    rows = _read(DATA / "t12_1_out_of_fold_predictions.csv")
    assert len(rows) == 28 * 5
    expected = {
        (row["case_id"], candidate)
        for row in comparison
        for candidate in (
            "P0_frozen_rho_l1",
            "P1_eta",
            "P2_lambda_max",
            "P3_rho_l1",
            "P4_epsilon_a_d",
        )
    }
    assert {(row["case_id"], row["candidate"]) for row in rows} == expected
    assert all(row["stratum"] not in ("", "all") for row in rows)


def test_recommendation_is_coherent_with_preregistered_numeric_rule():
    rows = _read(DATA / "t12_1_predictor_diagnostics.csv")
    global_rows = {
        row["candidate"]: row for row in rows if row["record_type"] == "global_oof"
    }
    recommendation = next(
        row["recommendation"] for row in rows if row["record_type"] == "recommendation"
    )
    p3 = global_rows["P3_rho_l1"]
    assert float(p3["rmse_log"]) <= np.log(2.0)
    assert float(p3["fraction_within_factor_2"]) >= 0.8
    best = min(float(global_rows[name]["rmse_log"]) for name in ("P1_eta", "P2_lambda_max", "P3_rho_l1"))
    assert float(p3["rmse_log"]) <= best + 0.05
    assert recommendation == "READY_T12_2_RHO1_RECALIBRATION_STUDY"


def test_all_six_csvs_are_finite_and_have_unique_headers():
    paths = (
        DATA / "t12_1_extended_convergence.csv",
        DATA / "t12_1_convergence_summary.csv",
        DATA / "t12_1_resolved_comparison.csv",
        DATA / "t12_1_mechanism_diagnostics.csv",
        DATA / "t12_1_predictor_diagnostics.csv",
        DATA / "t12_1_out_of_fold_predictions.csv",
    )
    for path in paths:
        text = path.read_text(encoding="utf-8")
        header = text.splitlines()[0].split(",")
        assert len(header) == len(set(header))
        assert "\x00" not in text
        assert "nan" not in text.lower()
        assert "inf" not in text.lower()
        assert text.endswith("\n")


def test_t12_artifact_hashes_are_preserved():
    for name, expected in T12_HASHES.items():
        assert _hash(DATA / name) == expected
    figure = ROOT / "results" / "figures" / "t12_model_e_sentinel_audit.png"
    assert _hash(figure) == "3ad7ab8569640e98e86bf405e3740630ccf043941aa3c83c7d66145ff1124059"


def test_analyze_only_is_deterministic_and_does_not_change_raw_extension():
    raw = DATA / "t12_1_extended_convergence.csv"
    raw_before = _hash(raw)
    before = {_path.name: _hash(_path) for _path in DERIVED}
    subprocess.run(
        [sys.executable, str(SCRIPT), "--analyze-only"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert _hash(raw) == raw_before
    assert {_path.name: _hash(_path) for _path in DERIVED} == before
