"""Validation, persistence, and deterministic regressions for T07."""

import csv
from pathlib import Path

import numpy as np
import pytest

from acoustic_ms import solve_multipolar_nodal, solve_multipolar_nodal_interaction_forces


ROOT = Path(__file__).resolve().parents[1]


def _pair():
    return np.array([[-1.05, 0.0, 0.0], [1.05, 0.0, 0.0]])


@pytest.mark.parametrize(
    "positions",
    [
        np.zeros((2, 2)),
        np.zeros((2, 4)),
        np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0], [2.1, 0.0, 0.1]]),
        np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0]]),
    ],
)
def test_multipolar_position_validation(positions):
    with pytest.raises(ValueError):
        solve_multipolar_nodal(positions, 0.1, 1.0, 0.0, 0.8, 3)


@pytest.mark.parametrize("lmax", [0, -1, 1.5, True])
def test_multipolar_order_validation(lmax):
    with pytest.raises(ValueError):
        solve_multipolar_nodal(_pair(), 0.1, 1.0, 0.0, 0.8, lmax)


@pytest.mark.parametrize(
    "parameters",
    [(0.0, 1.0, 0.8), (0.1, 0.0, 0.8), (0.1, 1.0, np.inf), (0.11, 1.0, 0.8)],
)
def test_multipolar_scalar_validation(parameters):
    k, radius, f1 = parameters
    with pytest.raises(ValueError):
        solve_multipolar_nodal(_pair(), k, radius, 0.0, f1, 3)


def test_repeated_solution_is_deterministic():
    first = solve_multipolar_nodal_interaction_forces(_pair(), 0.1, 1, 1, 0, 0.8, 5)
    second = solve_multipolar_nodal_interaction_forces(_pair(), 0.1, 1, 1, 0, 0.8, 5)
    assert np.array_equal(first.forces_xy, second.forces_xy)
    assert np.array_equal(first.solution.coefficients, second.solution.coefficients)


def test_t07_csv_integrity_and_counts():
    paths = {
        "t07_pair_analytic_validation.csv": 3,
        "t07_cluster_convergence.csv": 35,
    }
    for name, count in paths.items():
        rows = list(csv.DictReader((ROOT / "results" / "data" / name).read_text(encoding="utf-8").splitlines()))
        assert len(rows) == count
        assert all("nan" not in str(row).lower() and "inf" not in str(row).lower() for row in rows)
    dimer_rows = list(csv.DictReader((ROOT / "results" / "data" / "t07_dimer_convergence.csv").read_text(encoding="utf-8").splitlines()))
    base = [row for row in dimer_rows if int(row["lmax"]) <= 9]
    assert len(base) == 5 * 4 * 5
    keys = [(row["distance_ratio"], row["f1"], row["lmax"]) for row in dimer_rows]
    assert len(keys) == len(set(keys))
    assert all("nan" not in str(row).lower() and "inf" not in str(row).lower() for row in dimer_rows)
