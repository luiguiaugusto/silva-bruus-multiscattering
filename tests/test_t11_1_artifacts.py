"""Artifact integrity checks for the T11.1 stability audit."""

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
STABILITY = DATA / "t11_1_solver_stability.csv"
HIGH_PRECISION = DATA / "t11_1_high_precision_oracle.csv"


def _rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_t11_1_artifact_counts_order_and_finiteness():
    stability = _rows(STABILITY)
    high_precision = _rows(HIGH_PRECISION)
    assert len(stability) == 48
    assert len(high_precision) == 14
    cases = [
        "dimer_axis", "dimer_diagonal", "dimer_rigid",
        "trimer_equilateral", "trimer_scalene", "quartet_irregular",
    ]
    assert list(dict.fromkeys(row["case_id"] for row in stability)) == cases
    for case_id in cases:
        assert [
            int(row["lmax"]) for row in stability if row["case_id"] == case_id
        ] == list(range(2, 10))
    assert [row["case_id"] for row in high_precision] == (
        ["dimer_axis"] * 7 + ["trimer_scalene"] * 7
    )
    assert all(row["production_solver"] == "balanced_sqrt" for row in stability)
    assert max(float(row["condition_number_balanced"]) for row in stability) < 10.0
    assert max(float(row["balanced_backward_error"]) for row in stability) < 1e-12
    assert max(float(row["effective_incident_closure_error"]) for row in stability) < 1e-12
    assert max(float(row["scattering_closure_error"]) for row in stability) < 1e-12
    for channel in (
        "total", "interaction", "external_scattered", "scattered_scattered"
    ):
        resolved = [
            float(row[f"{channel}_balanced_vs_scattered_relative"])
            for row in stability
            if row[f"{channel}_balanced_vs_scattered_relative_applicable"] == "true"
        ]
        assert resolved and max(resolved) < 1e-9
    for rows in (stability, high_precision):
        for row in rows:
            assert not any(
                value.lower() in {"nan", "inf", "+inf", "-inf"}
                for value in row.values()
            )
    assert (FIGURES / "t11_1_model_e_stability.png").stat().st_size > 0


def test_t11_1_high_precision_acceptance():
    rows = _rows(HIGH_PRECISION)
    coefficient_errors = [
        float(row["relative_error_or_absolute_if_unresolved"])
        for row in rows if row["quantity"] in {"q", "d", "b"}
    ]
    force_errors = [
        float(row["relative_error_or_absolute_if_unresolved"])
        for row in rows
        if row["quantity"] not in {"q", "d", "b"}
        and row["quantity_resolved"] == "true"
    ]
    assert max(coefficient_errors) < 1e-11
    assert force_errors
    assert max(force_errors) < 1e-10


def test_t11_and_t11_1_five_csvs_have_no_nan_or_inf():
    paths = (
        DATA / "t11_model_e_convergence.csv",
        DATA / "t11_force_oracle.csv",
        DATA / "t11_force_decomposition.csv",
        STABILITY,
        HIGH_PRECISION,
    )
    for path in paths:
        rows = _rows(path)
        assert rows
        for row in rows:
            for value in row.values():
                assert value.lower() not in {"nan", "inf", "+inf", "-inf"}
                try:
                    numeric = float(value)
                except ValueError:
                    continue
                assert np.isfinite(numeric)


def test_t11_convergence_csv_preserves_legacy_and_adds_balanced_diagnostics():
    rows = _rows(DATA / "t11_model_e_convergence.csv")
    required = {
        "system_residual",
        "condition_number",
        "scattered_condition_number",
        "balanced_condition_number",
        "balanced_backward_error",
        "effective_incident_closure_error",
        "scattering_closure_error",
        "production_solver",
    }
    assert required <= set(rows[0])
    assert len(rows) == 48
    assert all(row["production_solver"] == "balanced_sqrt" for row in rows)
