"""Integrity checks for the deterministic T12 sentinel artifacts."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"


def _read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _assert_finite_csv(name: str) -> None:
    rows = _read(name)
    assert rows
    for row in rows:
        for value in row.values():
            lowered = value.lower()
            assert "nan" not in lowered
            assert "inf" not in lowered


def test_t12_artifact_counts_and_order() -> None:
    manifest = _read("t12_sentinel_manifest.csv")
    comparisons = _read("t12_model_comparison.csv")
    convergence = _read("t12_model_e_convergence.csv")
    audit = _read("t12_threshold_audit.csv")
    assert len(manifest) == len(comparisons) == 28
    assert [int(row["sentinel_order"]) for row in manifest] == list(range(1, 29))
    assert [row["case_id"] for row in comparisons] == [row["case_id"] for row in manifest]
    assert set(row["case_id"] for row in convergence) == set(row["case_id"] for row in manifest)
    assert len(audit) == 26
    assert sum(row["record_type"] == "threshold" for row in audit) == 24
    assert sum(row["record_type"] == "frozen_prediction_performance" for row in audit) == 1
    assert sum(row["record_type"] == "gate" for row in audit) == 1


def test_t12_manifest_matches_preregistration_constraints() -> None:
    manifest = _read("t12_sentinel_manifest.csv")
    assert all(row["source"] == "pre_registered_t08_calibration" for row in manifest)
    assert all(row["validated"] == "true" for row in manifest)
    assert all(row["split"] == "calibration" for row in manifest)
    assert {int(row["particle_count"]) for row in manifest} == {2, 3, 4}
    assert len({(row["particle_count"], row["family"], row["rho_band"]) for row in manifest}) == 28


def test_t12_convergence_sequences_and_two_step_confirmation() -> None:
    rows = _read("t12_model_e_convergence.csv")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
    for case_rows in grouped.values():
        final = int(case_rows[-1]["final_lmax"])
        assert [int(row["lmax"]) for row in case_rows] == list(range(2, final + 1))
        assert 5 <= final <= 13
        assert all(row["campaign_complete"] == "true" for row in case_rows)
        assert all(row["production_solver"] == "balanced_sqrt" for row in case_rows)
        for channel in ("total", "interaction", "external_scattered", "scattered_scattered"):
            minimum = int(case_rows[-1][f"{channel}_minimum_confirmed_lmax"])
            confirmed = case_rows[-1][f"{channel}_confirmed"] == "true"
            assert confirmed == (minimum > 0)
            if confirmed:
                index = next(i for i, row in enumerate(case_rows) if int(row["lmax"]) == minimum)
                assert index >= 2
                for offset in (-1, 0):
                    row = case_rows[index + offset]
                    assert row[f"{channel}_change_applicable"] == "true"
                    assert float(row[f"{channel}_successive_change"]) <= 1e-5


def test_t12_accepted_numerical_diagnostics() -> None:
    rows = _read("t12_model_e_convergence.csv")
    assert all(float(row["balanced_condition_number"]) < 10.0 for row in rows)
    for field in (
        "balanced_backward_error",
        "effective_incident_closure_error",
        "scattering_closure_error",
        "force_decomposition_residual",
    ):
        assert all(float(row[field]) < 1e-12 for row in rows)
    comparisons = _read("t12_model_comparison.csv")
    assert all(float(row["decomposition_relative_error"]) < 1e-12 for row in comparisons)
    assert all(row["diagnostics_pass"] == "true" for row in comparisons)


def test_t12_applicability_and_false_safe_rules() -> None:
    rows = _read("t12_model_comparison.csv")
    thresholds = {1: 0.0053990295322641655, 5: 0.02000077753569526, 10: 0.03914887870730305}
    for row in rows:
        applicable = row["threshold_metric_applicable"] == "true"
        for percent, threshold in thresholds.items():
            predicted = float(row["rho_l1"]) <= threshold
            observed = applicable and float(row["epsilon_a_e"]) <= percent / 100.0
            assert (row[f"predicted_safe_{percent}pct"] == "true") == predicted
            assert (row[f"observed_safe_{percent}pct"] == "true") == observed
            assert (row[f"false_safe_{percent}pct"] == "true") == (applicable and predicted and not observed)
        if row["prediction_metric_applicable"] == "false":
            assert float(row["prediction_factor"]) == 0.0


def test_t12_csvs_contain_no_nan_or_inf() -> None:
    for name in (
        "t12_sentinel_manifest.csv",
        "t12_model_e_convergence.csv",
        "t12_model_comparison.csv",
        "t12_threshold_audit.csv",
    ):
        _assert_finite_csv(name)


def test_t12_figure_exists_and_is_nonempty() -> None:
    path = ROOT / "results" / "figures" / "t12_model_e_sentinel_audit.png"
    assert path.is_file()
    assert path.stat().st_size > 100_000
