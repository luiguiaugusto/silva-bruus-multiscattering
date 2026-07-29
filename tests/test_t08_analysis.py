"""Leakage, threshold, raw-data, and T07 regression tests for T08."""

import csv
from pathlib import Path

import numpy as np
import pytest

from acoustic_ms.transferability import conservative_threshold, select_predictor_by_group_cv


ROOT = Path(__file__).resolve().parents[1]


def _row(n, family, eta, epsilon):
    return {"particle_count": n, "family": family, "eta": eta, "lambda_max": 1.1 * eta, "rho_l1": 0.9 * eta, "epsilon_a": epsilon}


def test_selection_and_threshold_reject_holdout_leakage():
    calibration = [_row(n, family, 0.01 * index, 0.02 * index) for index, (n, family) in enumerate(
        [(2, "pair"), (3, "linear"), (3, "compact"), (3, "irregular"), (4, "linear"), (4, "compact"), (4, "irregular")], 1
    ) for _ in (0, 1)]
    selected, scores = select_predictor_by_group_cv(calibration)
    assert selected in scores and set(scores) == {"eta", "lambda_max", "rho_l1"}
    threshold, available, count = conservative_threshold(calibration, selected, 0.2, 8)
    assert available and count >= 8 and threshold > 0
    contaminated = calibration + [_row(6, "linear", 0.001, 0.9)]
    with pytest.raises(ValueError): select_predictor_by_group_cv(contaminated)
    with pytest.raises(ValueError): conservative_threshold(contaminated, selected, 0.2)



def test_threshold_never_splits_equal_predictor_values():
    rows = [_row(3, "linear", 0.01, 0.005) for _ in range(7)]
    rows.extend([
        _row(4, "compact", 0.02, 0.005),
        _row(4, "irregular", 0.02, 0.2),
    ])
    threshold, available, count = conservative_threshold(rows, "eta", 0.01, 8)
    assert threshold == 0.0
    assert not available
    assert count == 0
def test_raw_artifact_counts_finiteness_and_split():
    cases = list(csv.DictReader((ROOT / "results/data/t08_cases.csv").read_text().splitlines()))
    assert len(cases) == 312 and len({row["case_id"] for row in cases}) == 312
    assert sum(row["split"] == "calibration" for row in cases) == 168
    assert sum(row["split"] == "holdout" for row in cases) == 144
    numeric = ["eta", "lambda_max", "rho_l1", "epsilon_a", "epsilon_b", "maximum_physical_residual"]
    assert all(np.isfinite(float(row[field])) for row in cases for field in numeric)
    assert max(float(row["maximum_physical_residual"]) for row in cases) < 1e-11


def test_t07_canonical_force_regression():
    old = list(csv.DictReader((ROOT / "results/data/t07_cluster_convergence.csv").read_text().splitlines()))
    new = list(csv.DictReader((ROOT / "results/data/t08_convergence.csv").read_text().splitlines()))
    mapping = {
        "trimer_linear": (3, "linear"), "trimer_equilateral": (3, "compact"),
        "trimer_scalene": (3, "irregular"), "quartet_linear": (4, "linear"),
        "quartet_square": (4, "compact"), "quartet_irregular": (4, "irregular"),
    }
    for geometry, (n, family) in mapping.items():
        expected = [row for row in old if row["geometry"] == geometry]
        prefix = f"n{n}_{family}_f0.8_d2.1"
        observed = {int(row["lmax"]): row for row in new if row["case_id"] == prefix}
        for row in expected:
            lmax = int(row["lmax"])
            if lmax in observed:
                np.testing.assert_allclose(float(observed[lmax]["rms_d"]), float(row["total_force_rms"]), rtol=3e-12, atol=3e-13)
