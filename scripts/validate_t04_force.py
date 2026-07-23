"""Generate deterministic T04 two-particle interaction-force validation data."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from acoustic_ms import (
    corrected_nodal_pair_force_magnitude,
    nodal_pair_force_magnitude,
    rayleigh_scattering_coefficients,
    solve_rayleigh_nodal_interaction_forces,
    spherical_hankel1,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "data" / "t04_pair_force_validation.csv"
FIELDS = (
    "case", "ka", "d_over_a", "kd", "f0", "f1", "force_1_x", "force_1_y",
    "force_2_x", "force_2_y", "radial_force_ms_l1", "scalar_reference",
    "relative_error", "silva_bruus_force", "corrected_l5_force",
    "residual_relative", "condition_number",
)


def scalar_pair_reference(ka: float, separation: float, f1: float, energy: float) -> float:
    """Independent scalar L=1 pair expression from the T04 specification."""
    s1 = rayleigh_scattering_coefficients(ka, 0.0, f1)[1]
    kd = ka * separation
    s10 = s1 * np.sqrt(12.0 * np.pi) / (1.0 - s1 * (spherical_hankel1(0, kd) + spherical_hankel1(2, kd)))
    q = spherical_hankel1(1, kd, derivative=True) / kd - spherical_hankel1(1, kd) / kd**2
    return float(4.0 * np.pi * ka * energy * np.sqrt(3.0 / (4.0 * np.pi)) * np.real(np.conj(f1) * s10 * q))


def pair_row(ka: float, separation: float, f0: float, f1: float, energy: float = 1.0) -> dict[str, object]:
    positions = np.array([[-separation / 2, 0.0, 0.0], [separation / 2, 0.0, 0.0]])
    result = solve_rayleigh_nodal_interaction_forces(positions, ka, 1.0, energy, f0, f1)
    radial = -float(result.forces_xy[0, 0])
    reference = scalar_pair_reference(ka, separation, f1, energy)
    return {
        "case": "pair", "ka": ka, "d_over_a": separation, "kd": ka * separation,
        "f0": f0, "f1": f1, "force_1_x": result.forces_xy[0, 0], "force_1_y": result.forces_xy[0, 1],
        "force_2_x": result.forces_xy[1, 0], "force_2_y": result.forces_xy[1, 1],
        "radial_force_ms_l1": radial, "scalar_reference": reference,
        "relative_error": abs(radial - reference) / abs(reference),
        "silva_bruus_force": nodal_pair_force_magnitude(ka, 1.0, separation, energy, f1),
        "corrected_l5_force": corrected_nodal_pair_force_magnitude(ka, 1.0, separation, energy, f1),
        "residual_relative": result.solution.residual_relative,
        "condition_number": result.solution.condition_number,
    }


def main() -> None:
    rows: list[dict[str, object]] = []
    single = solve_rayleigh_nodal_interaction_forces([[0.0, 0.0, 0.0]], .1, 1.0, 1.0, 0.0, .8)
    rows.append({field: "" for field in FIELDS} | {"case": "single", "ka": .1, "f0": 0.0, "f1": .8, "force_1_x": 0.0, "force_1_y": 0.0, "residual_relative": single.solution.residual_relative, "condition_number": single.solution.condition_number})
    rows.extend(pair_row(.1, 2.0, .2, f1) for f1 in (.1, .4, .8, 1.0))
    rows.extend((pair_row(.08, 4.0, -.5, .8), pair_row(.05, 8.0, .3, .4)))
    if max(float(row["relative_error"]) for row in rows[1:]) > 2e-12:
        raise RuntimeError("T04 scalar pair validation exceeded tolerance")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"T04: {len(rows) - 1} pair benchmarks; max relative error={max(float(row['relative_error']) for row in rows[1:]):.3e}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
