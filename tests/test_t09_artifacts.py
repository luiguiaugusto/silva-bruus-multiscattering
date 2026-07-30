"""Artifact integrity checks for the T09 analytical audit."""

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"


def _rows(name):
    return list(csv.DictReader((DATA / name).read_text().splitlines()))


def test_operator_audit_reproduces_all_frozen_t08_radii():
    rows = _rows("t09_operator_audit.csv")
    assert len(rows) == 312
    assert len({row["case_id"] for row in rows}) == 312
    assert max(float(row["rho_absolute_difference"]) for row in rows) < 1e-14
    assert all(
        float(row["rho_analytic"]) <= float(row["spectral_norm"]) + 1e-14
        for row in rows
    )
    assert all(
        float(row["rho_analytic"]) <= float(row["infinity_norm"]) + 1e-14
        for row in rows
    )


def test_neumann_artifact_has_three_complete_monotone_sequences():
    rows = _rows("t09_neumann_convergence.csv")
    assert len(rows) == 39
    case_ids = sorted({row["case_id"] for row in rows})
    assert len(case_ids) == 3
    for case_id in case_ids:
        selected = [row for row in rows if row["case_id"] == case_id]
        orders = [int(row["partial_order"]) for row in selected]
        errors = [float(row["relative_solution_error"]) for row in selected]
        assert orders == list(range(13))
        resolved = [error for error in errors if error > 1e-14]
        assert all(
            later < earlier for earlier, later in zip(resolved, resolved[1:])
        )
        assert errors[-1] < errors[0]


def test_t09_csvs_are_finite_and_summary_checks_pass():
    names = (
        "t09_operator_audit.csv",
        "t09_neumann_convergence.csv",
        "t09_analytic_summary.csv",
    )
    for name in names:
        for row in _rows(name):
            for value in row.values():
                assert value.lower() not in {"nan", "inf", "+inf", "-inf"}
                try:
                    numeric = float(value)
                except ValueError:
                    continue
                assert np.isfinite(numeric)
    summary = {
        row["metric"]: float(row["value"])
        for row in _rows("t09_analytic_summary.csv")
    }
    assert summary["symbolic_hankel_identity"] == 1.0
    assert summary["symbolic_near_field_series"] == 1.0
    assert summary["maximum_rho"] < 1.0
