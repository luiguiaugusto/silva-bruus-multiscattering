"""Produce deterministic validation data for the T03 Rayleigh field solver."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from acoustic_ms.scattering import rayleigh_scattering_coefficients
from acoustic_ms.special import spherical_hankel1
from acoustic_ms.solver import solve_rayleigh_nodal


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "data" / "t03_solver_validation.csv"


def pair_error(solution, ka: float, separation: float, f0: float, f1: float) -> float:
    s1 = rayleigh_scattering_coefficients(ka, f0, f1)[1]
    expected = s1 * np.sqrt(12.0 * np.pi) / (1.0 - s1 * (spherical_hankel1(0, ka * separation) + spherical_hankel1(2, ka * separation)))
    return float(np.max(abs(solution.coefficients[:, 2] - expected)) / abs(expected))


def main() -> None:
    rows: list[dict[str, object]] = []
    ka, f0, f1 = .08, -.4, .7
    single = solve_rayleigh_nodal([[0.0, 0.0, 0.0]], ka, 1.0, f0, f1)
    rows.append({"case": "single", "n_particles": 1, "ka": ka, "d_over_a": "", "f0": f0, "f1": f1, "residual_relative": single.residual_relative, "condition_number": single.condition_number, "benchmark_relative_error": float(abs(single.coefficients[0, 2] - 1j * f1 * np.sqrt(np.pi / 3.0) * ka**3) / abs(1j * f1 * np.sqrt(np.pi / 3.0) * ka**3))})
    for separation in (20.0, 4.0, 2.0):
        pair = solve_rayleigh_nodal([[-separation / 2, 0.0, 0.0], [separation / 2, 0.0, 0.0]], ka, 1.0, f0, f1)
        rows.append({"case": "pair", "n_particles": 2, "ka": ka, "d_over_a": separation, "f0": f0, "f1": f1, "residual_relative": pair.residual_relative, "condition_number": pair.condition_number, "benchmark_relative_error": pair_error(pair, ka, separation, f0, f1)})
    triangle = solve_rayleigh_nodal([[0.0, 0.0, 0.0], [3.0, .4, 0.0], [-.8, 2.6, 0.0]], ka, 1.0, f0, f1)
    rows.append({"case": "scalene_triangle", "n_particles": 3, "ka": ka, "d_over_a": "", "f0": f0, "f1": f1, "residual_relative": triangle.residual_relative, "condition_number": triangle.condition_number, "benchmark_relative_error": ""})
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(f"{row['case']}: N={row['n_particles']}, residual={float(row['residual_relative']):.3e}, cond={float(row['condition_number']):.3e}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
