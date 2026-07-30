"""Integrity, scientific regression, and determinism checks for T11 outputs."""

import csv
import hashlib
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
ARTIFACTS = (
    DATA / "t11_model_e_convergence.csv",
    DATA / "t11_force_oracle.csv",
    DATA / "t11_force_decomposition.csv",
    FIGURES / "t11_model_e_validation.png",
)


def _rows(name):
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _hashes():
    return tuple(hashlib.sha256(path.read_bytes()).hexdigest() for path in ARTIFACTS)


def test_t11_artifact_counts_order_and_finiteness():
    convergence = _rows("t11_model_e_convergence.csv")
    oracle = _rows("t11_force_oracle.csv")
    decomposition = _rows("t11_force_decomposition.csv")
    assert len(convergence) == 48
    assert len(oracle) == 72
    assert len(decomposition) == 16
    expected_cases = [
        "dimer_axis", "dimer_diagonal", "dimer_rigid",
        "trimer_equilateral", "trimer_scalene", "quartet_irregular",
    ]
    assert list(dict.fromkeys(row["case_id"] for row in convergence)) == expected_cases
    for case_id in expected_cases:
        assert [int(row["lmax"]) for row in convergence if row["case_id"] == case_id] == list(range(2, 10))
    for rows in (convergence, oracle, decomposition):
        for row in rows:
            assert not any(value.lower() in {"nan", "inf", "+inf", "-inf"} for value in row.values())


def test_t11_oracle_and_decomposition_regressions():
    oracle = _rows("t11_force_oracle.csv")
    resolved_errors = [
        float(row["relative_error_or_absolute_if_unresolved"])
        for row in oracle if row["force_resolved"] == "true"
    ]
    assert resolved_errors
    assert max(resolved_errors) < 1e-10
    decomposition = _rows("t11_force_decomposition.csv")
    assert max(float(row["decomposition_residual"]) for row in decomposition) < 1e-15
    assert all(float(row["scattered_scattered_over_interaction_rms"]) > 0.0 for row in decomposition)


def test_t11_convergence_flags_require_two_successive_changes():
    rows = _rows("t11_model_e_convergence.csv")
    channels = ("total", "interaction", "external_scattered", "scattered_scattered")
    for case_id in dict.fromkeys(row["case_id"] for row in rows):
        case = [row for row in rows if row["case_id"] == case_id]
        for channel in channels:
            reported = int(case[-1][f"{channel}_minimum_confirmed_lmax"])
            qualifying = []
            for index in range(2, len(case)):
                if (
                    case[index - 1][f"{channel}_change_applicable"] == "true"
                    and case[index][f"{channel}_change_applicable"] == "true"
                    and float(case[index - 1][f"{channel}_successive_change"]) <= 1e-5
                    and float(case[index][f"{channel}_successive_change"]) <= 1e-5
                ):
                    qualifying.append(int(case[index]["lmax"]))
            assert reported == (qualifying[0] if qualifying else 0)


def test_t11_analysis_is_byte_deterministic():
    command = [sys.executable, str(ROOT / "scripts" / "analyze_t11_model_e.py")]
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
    first = _hashes()
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
    assert _hashes() == first
    assert all(path.stat().st_size > 0 for path in ARTIFACTS)
