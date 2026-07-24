"""Generate deterministic T06 quartet body-expansion artifacts."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    angular_errors_degrees,
    decompose_nodal_quartet,
    irregular_quartet,
    linear_quartet,
    rms_relative_error,
    rms_vector_magnitude,
    square_quartet,
    symmetric_particle_errors,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
GEOMETRIES = (
    ("linear_chain", linear_quartet),
    ("square", square_quartet),
    ("irregular", irregular_quartet),
)
F1_VALUES = (.1, .4, .8, 1.0)
DISTANCES = np.linspace(2.1, 10.0, 160)


def _maximum_defined(values):
    return float(np.nanmax(values)) if np.any(~np.isnan(values)) else float("nan")


def _closed_four_body(result):
    closed = result.model_c_forces_xy.copy()
    for indices, comparison in zip(
        result.triplet_indices, result.triplet_comparisons, strict=True
    ):
        closed[list(indices)] -= comparison.model_c_forces_xy
    return closed + result.model_b_forces_xy


def _summary_row(geometry_name, geometry, distance, f1):
    result = decompose_nodal_quartet(geometry(distance), .1, 1., 1., 0., f1)
    comparisons = (
        ("a_vs_c", result.model_a_forces_xy),
        ("b_vs_c", result.model_b_forces_xy),
        ("up_to_three_vs_c", result.up_to_three_body_forces_xy),
    )
    row = {
        "geometry": geometry_name,
        "ka": .1,
        "f0": 0.,
        "f1": f1,
        "d_min_over_a": distance,
    }
    for label, model in comparisons:
        row[f"rms_{label}"] = rms_relative_error(result.model_c_forces_xy, model)
        row[f"max_symmetric_{label}"] = float(
            np.max(symmetric_particle_errors(result.model_c_forces_xy, model))
        )
        row[f"max_angle_{label}_deg"] = _maximum_defined(
            angular_errors_degrees(result.model_c_forces_xy, model)
        )
    full_rms = rms_vector_magnitude(result.model_c_forces_xy)
    row.update({
        "rms_two_body_correction": rms_vector_magnitude(result.two_body_correction_xy),
        "rms_collective_correction": rms_vector_magnitude(result.collective_correction_xy),
        "rms_three_body_sum": rms_vector_magnitude(result.irreducible_three_body_sum_xy),
        "rms_four_body": rms_vector_magnitude(result.irreducible_four_body_xy),
        "relative_three_body_sum_amplitude": rms_vector_magnitude(result.irreducible_three_body_sum_xy) / full_rms,
        "relative_four_body_amplitude": rms_vector_magnitude(result.irreducible_four_body_xy) / full_rms,
        "reconstruction_max_abs": float(np.max(np.abs(
            result.model_c_forces_xy
            - result.model_b_forces_xy
            - result.irreducible_three_body_sum_xy
            - result.irreducible_four_body_xy
        ))),
        "closed_form_max_abs": float(np.max(np.abs(
            _closed_four_body(result) - result.irreducible_four_body_xy
        ))),
        "quartet_residual_relative": result.full_comparison.global_result.solution.residual_relative,
        "quartet_condition_number": result.full_comparison.global_result.solution.condition_number,
        "max_triplet_residual_relative": max(
            comparison.global_result.solution.residual_relative
            for comparison in result.triplet_comparisons
        ),
        "max_triplet_condition_number": max(
            comparison.global_result.solution.condition_number
            for comparison in result.triplet_comparisons
        ),
        "sum_c_x": float(np.sum(result.model_c_forces_xy[:, 0])),
        "sum_c_y": float(np.sum(result.model_c_forces_xy[:, 1])),
    })
    return row, result


def _regression_rows():
    rows = []
    for geometry_name, geometry in GEOMETRIES:
        _, result = _summary_row(geometry_name, geometry, 2.1, .8)
        for particle in range(4):
            row = {
                "geometry": geometry_name,
                "particle": particle,
                "ka": .1,
                "f0": 0.,
                "f1": .8,
                "d_min_over_a": 2.1,
            }
            fields = (
                ("a", result.model_a_forces_xy),
                ("b", result.model_b_forces_xy),
                ("c", result.model_c_forces_xy),
                ("delta2", result.two_body_correction_xy),
                ("collective", result.collective_correction_xy),
                ("phi3_sum", result.irreducible_three_body_sum_xy),
                ("up_to_three", result.up_to_three_body_forces_xy),
                ("phi4", result.irreducible_four_body_xy),
            )
            for label, values in fields:
                row[f"{label}_x"], row[f"{label}_y"] = values[particle]
            for triplet_row, indices in enumerate(result.triplet_indices):
                label = "".join(str(index) for index in indices)
                values = result.irreducible_three_body_by_triplet_xy[
                    triplet_row, particle
                ]
                row[f"phi3_{label}_x"], row[f"phi3_{label}_y"] = values
            rows.append(row)
    return rows


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _plot(sweep_rows, fields, titles, output, suptitle=None):
    figure, axes = plt.subplots(
        3, 2, figsize=(10, 11), sharex=True, constrained_layout=True
    )
    colors = dict(zip(F1_VALUES, ("C0", "C1", "C2", "C3"), strict=True))
    if suptitle is not None:
        figure.suptitle(suptitle)
    for row_index, (geometry_name, _) in enumerate(GEOMETRIES):
        for column, (field, title) in enumerate(zip(fields, titles, strict=True)):
            axis = axes[row_index, column]
            for f1 in F1_VALUES:
                rows = [
                    row for row in sweep_rows
                    if row["geometry"] == geometry_name and row["f1"] == f1
                ]
                axis.plot(
                    [row["d_min_over_a"] for row in rows],
                    [100 * row[field] for row in rows],
                    color=colors[f1],
                    label=fr"$f_1={f1}$",
                )
            axis.set_title(f"{geometry_name}: {title}")
            axis.set_ylabel("RMS measure (%)")
            axis.grid(True, alpha=.3)
            panel = [
                row[field] for row in sweep_rows
                if row["geometry"] == geometry_name
            ]
            if all(value > 0 for value in panel):
                axis.set_yscale("log")
            if row_index == 2:
                axis.set_xlabel(r"$d_{\min}/a$")
            if row_index == 0 and column == 1:
                axis.legend()
    figure.savefig(FIGURES / output, dpi=200)
    plt.close(figure)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    regression_rows = _regression_rows()
    sweep_rows = []
    for geometry_name, geometry in GEOMETRIES:
        for f1 in F1_VALUES:
            for distance in DISTANCES:
                row, _ = _summary_row(
                    geometry_name, geometry, float(distance), f1
                )
                required = [
                    value for key, value in row.items()
                    if key not in {"geometry"} and "angle" not in key
                ]
                if not np.all(np.isfinite(required)):
                    raise RuntimeError("non-finite mandatory T06 metric")
                if row["reconstruction_max_abs"] > 5e-15:
                    raise RuntimeError("quartet reconstruction failed")
                if row["closed_form_max_abs"] > 5e-15:
                    raise RuntimeError("closed four-body identity failed")
                sweep_rows.append(row)
    _write_csv(DATA / "t06_quartet_regression.csv", regression_rows)
    _write_csv(DATA / "t06_quartet_sweep.csv", sweep_rows)
    _plot(
        sweep_rows,
        ("rms_a_vs_c", "rms_b_vs_c"),
        ("A (Silva--Bruus) vs C", "B (isolated-pair MS) vs C"),
        "t06_quartet_model_errors.png",
    )
    _plot(
        sweep_rows,
        ("relative_three_body_sum_amplitude", "relative_four_body_amplitude"),
        (
            r"$F_{\rm RMS}(\Phi_\Sigma^{(3)})/F_{\rm RMS}(F^C)$",
            r"$F_{\rm RMS}(\Phi^{(4)})/F_{\rm RMS}(F^C)$",
        ),
        "t06_quartet_body_decomposition.png",
        "RMS amplitude ratios; not additive force fractions",
    )
    print(
        f"T06: wrote {len(regression_rows)} regression rows and "
        f"{len(sweep_rows)} sweep rows"
    )


if __name__ == "__main__":
    main()
