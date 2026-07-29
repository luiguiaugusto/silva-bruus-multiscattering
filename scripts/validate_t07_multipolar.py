#!/usr/bin/env python3
"""Generate deterministic T07 multipolar convergence artifacts."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    corrected_nodal_pair_force_magnitude,
    decompose_multipolar_cluster,
    equilateral_trimer,
    irregular_quartet,
    linear_quartet,
    linear_trimer,
    rms_vector_magnitude,
    scalene_trimer,
    solve_multipolar_nodal_interaction_forces,
    square_quartet,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
TOLERANCE = 1e-3
BASE_ORDERS = (1, 3, 5, 7, 9)


def _pair(distance):
    return np.array([[-distance / 2, 0.0, 0.0], [distance / 2, 0.0, 0.0]])


def _relative_change(current, previous, reference=None):
    difference = rms_vector_magnitude(current - previous)
    denominator = rms_vector_magnitude(current if reference is None else reference)
    return float(difference / denominator) if denominator > 0 else 0.0


def _write_csv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _analytic_rows():
    rows = []
    for ka in (0.1, 0.05, 0.025):
        result = solve_multipolar_nodal_interaction_forces(
            _pair(2.1), ka, 1.0, 1.0, 0.0, 0.8, 5
        )
        force = float(result.forces_xy[0, 0])
        analytic = -corrected_nodal_pair_force_magnitude(ka, 1.0, 2.1, 1.0, 0.8)
        rows.append({
            "ka": ka, "distance_ratio": 2.1, "f1": 0.8, "lmax": 5,
            "full_mode_count": 36,
            "active_mode_count": len(result.solution.active_modes),
            "model_d_radial_force": format(force, ".17g"),
            "equation_30_radial_force": format(analytic, ".17g"),
            "relative_difference": format(abs(force - analytic) / abs(force), ".17g"),
            "physical_residual": format(result.solution.residual_relative, ".17g"),
            "balanced_condition": format(result.solution.condition_number, ".17g"),
            "raw_condition": format(result.solution.physical_condition_number, ".17g"),
        })
    return rows


def _dimer_rows():
    rows = []
    series = {}
    for distance in (2.0, 2.05, 2.1, 2.5, 3.0):
        for f1 in (0.1, 0.4, 0.8, 1.0):
            values = []
            previous = None
            for lmax in BASE_ORDERS:
                result = solve_multipolar_nodal_interaction_forces(
                    _pair(distance), 0.1, 1.0, 1.0, 0.0, f1, lmax
                )
                error = 0.0 if previous is None else _relative_change(result.forces_xy, previous)
                values.append([lmax, result, error, previous is not None])
                previous = result.forces_xy
            if values[-1][2] > TOLERANCE or values[-2][2] > TOLERANCE:
                lmax = 11
                result = solve_multipolar_nodal_interaction_forces(
                    _pair(distance), 0.1, 1.0, 1.0, 0.0, f1, lmax
                )
                values.append([lmax, result, _relative_change(result.forces_xy, previous), True])
            confirmed = 0
            for index in range(1, len(values) - 1):
                if values[index][2] <= TOLERANCE and values[index + 1][2] <= TOLERANCE:
                    confirmed = values[index][0]
                    break
            series[(distance, f1)] = values
            for lmax, result, error, applicable in values:
                rows.append({
                    "distance_ratio": distance, "f1": f1, "ka": 0.1,
                    "lmax": lmax, "full_mode_count": (lmax + 1) ** 2,
                    "active_mode_count": len(result.solution.active_modes),
                    "force_0_x": format(result.forces_xy[0, 0], ".17g"),
                    "force_0_y": format(result.forces_xy[0, 1], ".17g"),
                    "force_1_x": format(result.forces_xy[1, 0], ".17g"),
                    "force_1_y": format(result.forces_xy[1, 1], ".17g"),
                    "force_rms": format(rms_vector_magnitude(result.forces_xy), ".17g"),
                    "successive_force_error": format(error, ".17g"),
                    "successive_error_applicable": str(applicable).lower(),
                    "minimum_confirmed_lmax": confirmed,
                    "physical_residual": format(result.solution.residual_relative, ".17g"),
                    "balanced_condition": format(result.solution.condition_number, ".17g"),
                    "raw_condition": format(result.solution.physical_condition_number, ".17g"),
                })
    return rows, series


def _cluster_rows():
    geometries = (
        ("trimer_linear", linear_trimer(2.1)),
        ("trimer_equilateral", equilateral_trimer(2.1)),
        ("trimer_scalene", scalene_trimer(2.1)),
        ("quartet_linear", linear_quartet(2.1)),
        ("quartet_square", square_quartet(2.1)),
        ("quartet_irregular", irregular_quartet(2.1)),
    )
    rows = []
    series = {}
    for name, positions in geometries:
        values = []
        previous = None
        orders = list(BASE_ORDERS)
        for lmax in orders:
            expansion = decompose_multipolar_cluster(
                positions, 0.1, 1.0, 1.0, 0.0, 0.8, lmax
            )
            values.append(expansion)
        def connected_not_converged():
            current, prior = values[-1], values[-2]
            arrays = [(current.model_d_forces_xy, prior.model_d_forces_xy),
                      (current.irreducible_three_body_sum_xy, prior.irreducible_three_body_sum_xy)]
            if len(positions) == 4:
                arrays.append((current.irreducible_four_body_xy, prior.irreducible_four_body_xy))
            return any(_relative_change(a, b) > TOLERANCE for a, b in arrays)
        if connected_not_converged():
            orders.append(11)
            values.append(decompose_multipolar_cluster(
                positions, 0.1, 1.0, 1.0, 0.0, 0.8, 11
            ))
        if name == "quartet_linear" and connected_not_converged():
            orders.append(13)
            values.append(decompose_multipolar_cluster(
                positions, 0.1, 1.0, 1.0, 0.0, 0.8, 13
            ))
        series[name] = (orders, values)
        for index, (lmax, expansion) in enumerate(zip(orders, values)):
            applicable = index > 0
            prior = values[index - 1] if applicable else None
            force_error = _relative_change(expansion.model_d_forces_xy, prior.model_d_forces_xy) if applicable else 0.0
            three_error = _relative_change(expansion.irreducible_three_body_sum_xy, prior.irreducible_three_body_sum_xy) if applicable else 0.0
            three_delta = _relative_change(expansion.irreducible_three_body_sum_xy, prior.irreducible_three_body_sum_xy, expansion.model_d_forces_xy) if applicable else 0.0
            four_applicable = applicable and len(positions) == 4
            four_error = _relative_change(expansion.irreducible_four_body_xy, prior.irreducible_four_body_xy) if four_applicable else 0.0
            four_delta = _relative_change(expansion.irreducible_four_body_xy, prior.irreducible_four_body_xy, expansion.model_d_forces_xy) if four_applicable else 0.0
            subset_results = [result for _, result in expansion.subset_results]
            rows.append({
                "geometry": name, "particle_count": len(positions), "ka": 0.1,
                "f1": 0.8, "distance_ratio": 2.1, "lmax": lmax,
                "full_modes_per_particle": (lmax + 1) ** 2,
                "active_modes_per_particle": len(expansion.full_result.solution.active_modes),
                "total_force_rms": format(rms_vector_magnitude(expansion.model_d_forces_xy), ".17g"),
                "three_body_sum_rms": format(rms_vector_magnitude(expansion.irreducible_three_body_sum_xy), ".17g"),
                "four_body_rms": format(rms_vector_magnitude(expansion.irreducible_four_body_xy), ".17g"),
                "successive_total_error": format(force_error, ".17g"),
                "successive_three_body_error": format(three_error, ".17g"),
                "three_body_change_over_total": format(three_delta, ".17g"),
                "successive_four_body_error": format(four_error, ".17g"),
                "four_body_change_over_total": format(four_delta, ".17g"),
                "successive_error_applicable": str(applicable).lower(),
                "four_body_applicable": str(four_applicable).lower(),
                "max_physical_residual": format(max(x.solution.residual_relative for x in subset_results), ".17g"),
                "max_balanced_condition": format(max(x.solution.condition_number for x in subset_results), ".17g"),
                "max_raw_condition": format(max(x.solution.physical_condition_number for x in subset_results), ".17g"),
            })
    return rows, series


def _plot_dimers(series):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, 5))
    styles = {0.1: "-", 0.4: "--", 0.8: "-.", 1.0: ":"}
    for color, distance in zip(colors, (2.0, 2.05, 2.1, 2.5, 3.0)):
        for f1 in (0.1, 0.4, 0.8, 1.0):
            values = series[(distance, f1)]
            orders = [value[0] for value in values[1:]]
            errors = [value[2] for value in values[1:]]
            axes[0].semilogy(orders, errors, color=color, linestyle=styles[f1], linewidth=1.2)
        final = values[-1][1]
        axes[1].semilogy([value[0] for value in values], [value[1].solution.condition_number for value in values], color=color, linewidth=1.3, label=f"d/a={distance:g}")
    for f1, style in styles.items():
        axes[0].plot([], [], color="gray", linestyle=style, label=fr"$f_1={f1:g}$")
    axes[0].axhline(TOLERANCE, color="black", linewidth=1, label=r"$10^{-3}$ criterion")
    axes[0].set(xlabel=r"$L_{\max}$", ylabel=r"successive force error $\epsilon_L^F$", title="Dimer force convergence (line style: contrast)")
    axes[1].set(xlabel=r"$L_{\max}$", ylabel="balanced condition number", title="Balanced systems (shown for $f_1=1$)")
    for axis in axes: axis.grid(True, alpha=0.25); axis.legend(fontsize=8)
    fig.savefig(FIGURES / "t07_dimer_convergence.png", dpi=220, metadata={"Software": "acoustic_ms T07"})
    plt.close(fig)


def _plot_clusters(series):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    colors = plt.cm.tab10.colors
    for color, (name, (orders, values)) in zip(colors, series.items()):
        xs = orders[1:]
        total = [_relative_change(values[i].model_d_forces_xy, values[i-1].model_d_forces_xy) for i in range(1, len(values))]
        three = [_relative_change(values[i].irreducible_three_body_sum_xy, values[i-1].irreducible_three_body_sum_xy) for i in range(1, len(values))]
        axes[0].semilogy(xs, total, color=color, marker="o", linewidth=1.2, label=name)
        axes[1].semilogy(xs, three, color=color, marker="o", linewidth=1.2, label=name + r" $\Phi^{(3)}$")
        if values[0].model_d_forces_xy.shape[0] == 4:
            four = [_relative_change(values[i].irreducible_four_body_xy, values[i-1].irreducible_four_body_xy) for i in range(1, len(values))]
            axes[1].semilogy(xs, four, color=color, marker="s", linestyle="--", linewidth=1.1, label=name + r" $\Phi^{(4)}$")
    for axis in axes:
        axis.axhline(TOLERANCE, color="black", linewidth=1)
        axis.set_xlabel(r"$L_{\max}$"); axis.grid(True, alpha=0.25); axis.legend(fontsize=7, ncol=2)
    axes[0].set_ylabel("successive relative change"); axes[0].set_title("Total Model-D force")
    axes[1].set_ylabel("successive relative change"); axes[1].set_title("Connected terms; squares denote four-body")
    fig.savefig(FIGURES / "t07_cluster_convergence.png", dpi=220, metadata={"Software": "acoustic_ms T07"})
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    analytic = _analytic_rows(); dimers, dimer_series = _dimer_rows(); clusters, cluster_series = _cluster_rows()
    _write_csv(DATA / "t07_pair_analytic_validation.csv", list(analytic[0]), analytic)
    _write_csv(DATA / "t07_dimer_convergence.csv", list(dimers[0]), dimers)
    _write_csv(DATA / "t07_cluster_convergence.csv", list(clusters[0]), clusters)
    _plot_dimers(dimer_series); _plot_clusters(cluster_series)
    print(f"T07: {len(analytic)} analytic rows, {len(dimers)} dimer rows, {len(clusters)} cluster rows")


if __name__ == "__main__":
    main()
