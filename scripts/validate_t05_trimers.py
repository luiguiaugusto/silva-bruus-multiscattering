"""Generate deterministic T05 Model A/B/C trimer regression and sweep artifacts."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    angular_errors_degrees, compare_nodal_force_models, equilateral_trimer,
    linear_trimer, rms_relative_error, rms_vector_magnitude, scalene_trimer, symmetric_particle_errors,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
GEOMETRIES = (("linear_chain", linear_trimer), ("equilateral", equilateral_trimer), ("scalene", scalene_trimer))
F1_VALUES = (.1, .4, .8, 1.0)
DISTANCES = np.linspace(2.1, 10.0, 160)


def _summary_rows(geometry_name, geometry, distance, f1):
    result = compare_nodal_force_models(geometry(distance), .1, 1., 1., 0., f1)
    errors_a = symmetric_particle_errors(result.model_c_forces_xy, result.model_a_forces_xy)
    errors_b = symmetric_particle_errors(result.model_c_forces_xy, result.model_b_forces_xy)
    angles_a = angular_errors_degrees(result.model_c_forces_xy, result.model_a_forces_xy)
    angles_b = angular_errors_degrees(result.model_c_forces_xy, result.model_b_forces_xy)
    return {
        "geometry": geometry_name, "ka": .1, "f0": 0., "f1": f1, "d_min_over_a": distance,
        "rms_a_vs_c": rms_relative_error(result.model_c_forces_xy, result.model_a_forces_xy),
        "rms_b_vs_c": rms_relative_error(result.model_c_forces_xy, result.model_b_forces_xy),
        "max_symmetric_a_vs_c": float(np.max(errors_a)), "max_symmetric_b_vs_c": float(np.max(errors_b)),
        "max_angle_a_vs_c_deg": float(np.nanmax(angles_a)) if np.any(~np.isnan(angles_a)) else float("nan"),
        "max_angle_b_vs_c_deg": float(np.nanmax(angles_b)) if np.any(~np.isnan(angles_b)) else float("nan"),
        "rms_two_body_correction": rms_vector_magnitude(result.two_body_correction_xy),
        "rms_irreducible_multibody": rms_vector_magnitude(result.irreducible_multibody_xy),
        "residual_relative": result.global_result.solution.residual_relative,
        "condition_number": result.global_result.solution.condition_number,
        "sum_c_x": float(np.sum(result.model_c_forces_xy[:, 0])), "sum_c_y": float(np.sum(result.model_c_forces_xy[:, 1])),
    }, result


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    regression_rows = []
    for name, geometry, distance in (("linear_chain", linear_trimer, 2.1), ("equilateral", equilateral_trimer, 2.1), ("scalene", scalene_trimer, 2.2)):
        _, result = _summary_rows(name, geometry, distance, .8)
        for index in range(3):
            regression_rows.append({"geometry": name, "particle": index, "ka": .1, "f0": 0., "f1": .8, "d_min_over_a": distance,
                "a_x": result.model_a_forces_xy[index,0], "a_y": result.model_a_forces_xy[index,1], "b_x": result.model_b_forces_xy[index,0], "b_y": result.model_b_forces_xy[index,1], "c_x": result.model_c_forces_xy[index,0], "c_y": result.model_c_forces_xy[index,1], "delta2_x": result.two_body_correction_xy[index,0], "delta2_y": result.two_body_correction_xy[index,1], "delta3_x": result.irreducible_multibody_xy[index,0], "delta3_y": result.irreducible_multibody_xy[index,1]})
    sweep_rows = []
    for name, geometry in GEOMETRIES:
        for f1 in F1_VALUES:
            for distance in DISTANCES:
                row, _ = _summary_rows(name, geometry, float(distance), f1)
                if not np.all(np.isfinite([row["rms_a_vs_c"], row["rms_b_vs_c"], row["residual_relative"], row["condition_number"]])):
                    raise RuntimeError("non-finite T05 production metric")
                sweep_rows.append(row)
    for path, rows in ((DATA / "t05_trimer_regression.csv", regression_rows), (DATA / "t05_trimer_sweep.csv", sweep_rows)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    figure, axes = plt.subplots(3, 2, figsize=(10, 11), sharex=True, constrained_layout=True)
    colors = {value: color for value, color in zip(F1_VALUES, ("C0", "C1", "C2", "C3"))}
    for row_index, (name, _) in enumerate(GEOMETRIES):
        for column, field in enumerate(("rms_a_vs_c", "rms_b_vs_c")):
            axis = axes[row_index, column]
            for f1 in F1_VALUES:
                rows = [row for row in sweep_rows if row["geometry"] == name and row["f1"] == f1]
                axis.plot([row["d_min_over_a"] for row in rows], [100 * row[field] for row in rows], color=colors[f1], label=fr"$f_1={f1}$")
            axis.set_title(f"{name}: {'A (Silva--Bruus)' if column == 0 else 'B (isolated-pair MS)'} vs C")
            axis.set_ylabel("RMS error (%)"); axis.grid(True, alpha=.3)
            if all(row[field] > 0 for row in sweep_rows if row["geometry"] == name): axis.set_yscale("log")
            if row_index == 2: axis.set_xlabel(r"$d_{\min}/a$")
            if row_index == 0 and column == 1: axis.legend()
    figure.savefig(FIGURES / "t05_trimer_model_errors.png", dpi=200)
    plt.close(figure)
    print(f"T05: wrote {len(regression_rows)} regression rows and {len(sweep_rows)} sweep rows")


if __name__ == "__main__":
    main()
