#!/usr/bin/env python3
"""Leakage-safe calibration and holdout analysis for T08 raw CSVs."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata

from acoustic_ms.transferability import (
    conservative_threshold,
    fit_transferability_power_law,
    select_predictor_by_group_cv,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
PREDICTORS = ("eta", "lambda_max", "rho_l1")


def _read_cases():
    rows = list(csv.DictReader((DATA / "t08_cases.csv").read_text(encoding="utf-8").splitlines()))
    if len(rows) != 312 or len({row["case_id"] for row in rows}) != 312:
        raise RuntimeError("t08_cases.csv must contain 312 unique cases")
    return rows


def _true(row, field):
    return row[field] == "true"


def _write(path, rows):
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".17g") if isinstance(value, float) else value for key, value in row.items()})


def _eligible_a(rows, split):
    return [row for row in rows if row["split"] == split and _true(row, "total_converged") and _true(row, "metric_applicable") and float(row["epsilon_a"]) > 0]


def _eligible_b(rows, split):
    return [row for row in rows if row["split"] == split and int(row["particle_count"]) > 2 and _true(row, "joint_converged") and _true(row, "collective_residual_resolved") and float(row["epsilon_b"]) > 0]


def _fit_row(predictor, response, rows):
    fit = fit_transferability_power_law(
        [float(row[predictor]) for row in rows], [float(row[response]) for row in rows]
    )
    value = fit.power_law
    return {
        "record_type": "calibration_fit", "predictor": predictor, "response": response,
        "scope": "calibration_n_le_4", "point_count": value.point_count,
        "prefactor": value.prefactor, "exponent": value.exponent,
        "r_squared_log": value.r_squared_log, "rmse_log": value.rmse_log,
        "max_abs_log_residual": value.max_abs_log_residual, "spearman": fit.spearman,
    }


def _prediction_metrics(rows, predictor, fit_row, scope):
    observed = np.array([float(row["epsilon_a"]) for row in rows])
    predicted = float(fit_row["prefactor"]) * np.array([float(row[predictor]) for row in rows]) ** float(fit_row["exponent"])
    log_residual = np.log(observed / predicted)
    factors = np.exp(np.abs(log_residual))
    ranks_observed = rankdata(observed)
    ranks_predicted = rankdata(predicted)
    spearman = float(np.corrcoef(ranks_observed, ranks_predicted)[0, 1]) if len(rows) > 1 else 1.0
    return {
        "record_type": "holdout_performance", "predictor": predictor,
        "response": "epsilon_a", "scope": scope, "point_count": len(rows),
        "rmse_log": float(np.sqrt(np.mean(log_residual**2))),
        "median_error_factor": float(np.median(factors)),
        "p90_error_factor": float(np.percentile(factors, 90)),
        "maximum_error_factor": float(np.max(factors)),
        "fraction_within_factor_2": float(np.mean(factors <= 2)), "spearman": spearman,
    }


def _analysis(rows):
    calibration_a = _eligible_a(rows, "calibration")
    calibration_b = _eligible_b(rows, "calibration")
    fit_rows = []
    fits = {}
    for predictor in PREDICTORS:
        for response, eligible in (("epsilon_a", calibration_a), ("epsilon_b", calibration_b)):
            row = _fit_row(predictor, response, eligible)
            fit_rows.append(row); fits[(predictor, response)] = row
    selected, cv_scores = select_predictor_by_group_cv(calibration_a)
    for predictor in PREDICTORS:
        fit_rows.append({
            "record_type": "cross_validation", "predictor": predictor,
            "response": "epsilon_a", "scope": "leave_n_family_out",
            "point_count": len(calibration_a), "rmse_log": cv_scores[predictor],
            "selected": str(predictor == selected).lower(),
        })
    holdout_all = _eligible_a(rows, "holdout")
    for predictor in PREDICTORS:
        for particle_count, scope in ((None, "holdout_all"), (6, "holdout_n6"), (10, "holdout_n10")):
            subset = holdout_all if particle_count is None else [row for row in holdout_all if int(row["particle_count"]) == particle_count]
            fit_rows.append(_prediction_metrics(subset, predictor, fits[(predictor, "epsilon_a")], scope))

    thresholds = []
    for tolerance in (0.01, 0.05, 0.10):
        threshold, available, calibration_count = conservative_threshold(calibration_a, selected, tolerance)
        eligible_holdout = holdout_all
        safe = [row for row in eligible_holdout if available and float(row[selected]) <= threshold]
        false_safe = [row for row in safe if float(row["epsilon_a"]) > tolerance]
        thresholds.append({
            "predictor": selected, "tolerance": tolerance,
            "threshold": threshold, "threshold_available": str(available).lower(),
            "calibration_count_below": calibration_count,
            "holdout_eligible": len(eligible_holdout), "holdout_predicted_safe": len(safe),
            "holdout_coverage": len(safe) / len(eligible_holdout) if eligible_holdout else 0.0,
            "false_safe_count": len(false_safe),
            "worst_safe_error": max((float(row["epsilon_a"]) for row in safe), default=0.0),
        })
    selected_holdout = next(row for row in fit_rows if row["record_type"] == "holdout_performance" and row["predictor"] == selected and row["scope"] == "holdout_all")
    threshold_five = next(row for row in thresholds if float(row["tolerance"]) == 0.05)
    holdout_rows = [row for row in rows if row["split"] == "holdout"]
    convergence_coverage = sum(_true(row, "total_converged") for row in holdout_rows) / len(holdout_rows)
    criterion = (
        convergence_coverage >= 0.8
        and float(selected_holdout["rmse_log"]) <= np.log(2)
        and float(selected_holdout["fraction_within_factor_2"]) >= 0.8
        and int(threshold_five["holdout_predicted_safe"]) >= 1
        and int(threshold_five["false_safe_count"]) == 0
    )
    for row in thresholds:
        row["criterion_supported"] = str(criterion).lower()
        row["holdout_total_convergence_coverage"] = convergence_coverage
    return fit_rows, thresholds, selected, fits


def _plot_predictors(rows, fits):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    colors = {"linear": "tab:blue", "compact": "tab:orange", "irregular": "tab:green", "pair": "tab:purple"}
    for column, predictor in enumerate(PREDICTORS):
        for row_index, response in enumerate(("epsilon_a", "epsilon_b")):
            axis = axes[row_index, column]
            markers = {2: "P", 3: "o", 4: "^", 6: "s", 10: "D"}
            for row in rows:
                particle_count = int(row["particle_count"])
                if response == "epsilon_b" and particle_count == 2:
                    continue
                if float(row[response]) <= 0:
                    continue
                if response == "epsilon_a":
                    eligible = _true(row, "total_converged")
                else:
                    eligible = _true(row, "joint_converged") and _true(row, "collective_residual_resolved")
                kwargs = {"s": 20, "alpha": 0.7}
                if not eligible:
                    kwargs.update(marker="x", color=colors[row["family"]])
                elif row["split"] == "calibration":
                    kwargs.update(marker=markers[particle_count], color=colors[row["family"]])
                else:
                    kwargs.update(marker=markers[particle_count], facecolors="none", edgecolors=colors[row["family"]])
                axis.scatter(float(row[predictor]), float(row[response]), **kwargs)
            fit = fits[(predictor, response)]
            x_values = np.logspace(np.log10(min(float(row[predictor]) for row in rows)), np.log10(max(float(row[predictor]) for row in rows)), 200)
            axis.plot(x_values, float(fit["prefactor"]) * x_values ** float(fit["exponent"]), color="black", linewidth=1.4)
            axis.set_xscale("log"); axis.set_yscale("log"); axis.grid(True, alpha=0.25)
            axis.set_xlabel(predictor); axis.set_ylabel(response)
    for family, color in colors.items():
        axes[0, 0].scatter([], [], color=color, marker="o", label=family)
    for particle_count, marker in {2: "P", 3: "o", 4: "^", 6: "s", 10: "D"}.items():
        axes[0, 0].scatter([], [], color="black", marker=marker, label=f"N={particle_count}")
    axes[0, 0].scatter([], [], facecolors="none", edgecolors="black", marker="s", label="holdout")
    axes[0, 0].scatter([], [], color="black", marker="x", label="unconfirmed")
    axes[0, 0].legend(fontsize=7, ncol=3)
    axes[0, 1].set_title("T08 predictor comparison (calibration-only fits)", fontsize=11)
    fig.savefig(FIGURES / "t08_predictor_comparison.png", dpi=220, metadata={"Software": "acoustic_ms T08"})
    plt.close(fig)


def _plot_transferability(rows, selected, fit):
    eligible = [row for row in rows if _true(row, "total_converged") and float(row["epsilon_a"]) > 0]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    series = (
        ("calibration", None, "o", "tab:blue"),
        ("holdout N=6", 6, "s", "tab:orange"),
        ("holdout N=10", 10, "D", "tab:red"),
    )
    for label, particle_count, marker, color in series:
        subset = [
            row for row in eligible
            if (
                (particle_count is None and row["split"] == "calibration")
                or (
                    row["split"] == "holdout"
                    and int(row["particle_count"]) == particle_count
                )
            )
        ]
        x = np.array([float(row[selected]) for row in subset])
        y = np.array([float(row["epsilon_a"]) for row in subset])
        predicted = float(fit["prefactor"]) * x ** float(fit["exponent"])
        axes[0].scatter(x, y, marker=marker, color=color, alpha=0.7, label=label)
        axes[1].scatter(predicted, y, marker=marker, color=color, alpha=0.7, label=label)
    for particle_count, marker, color in ((6, "s", "tab:orange"), (10, "D", "tab:red")):
        unconfirmed = [
            row for row in rows
            if (
                not _true(row, "total_converged")
                and int(row["particle_count"]) == particle_count
                and float(row["epsilon_a"]) > 0
            )
        ]
        if not unconfirmed:
            continue
        x = np.array([float(row[selected]) for row in unconfirmed])
        y = np.array([float(row["epsilon_a"]) for row in unconfirmed])
        predicted = float(fit["prefactor"]) * x ** float(fit["exponent"])
        axes[0].scatter(
            x, y, marker=marker, facecolors="none", edgecolors=color,
            linewidths=1.3, label=f"N={particle_count} unconfirmed",
        )
        axes[1].scatter(
            predicted, y, marker=marker, facecolors="none", edgecolors=color,
            linewidths=1.3, label=f"N={particle_count} unconfirmed",
        )
    xline = np.logspace(np.log10(min(float(row[selected]) for row in eligible)), np.log10(max(float(row[selected]) for row in eligible)), 200)
    axes[0].plot(xline, float(fit["prefactor"]) * xline ** float(fit["exponent"]), color="black", label="calibration fit")
    for tolerance in (0.01, 0.05, 0.10): axes[0].axhline(tolerance, linewidth=0.9, linestyle="--", label=f"{100*tolerance:g}%")
    observed = np.array([float(row["epsilon_a"]) for row in eligible]); predicted = float(fit["prefactor"]) * np.array([float(row[selected]) for row in eligible]) ** float(fit["exponent"])
    bounds = [min(observed.min(), predicted.min()), max(observed.max(), predicted.max())]
    axes[1].plot(bounds, bounds, color="black", label="identity"); axes[1].plot(bounds, np.array(bounds)*2, color="gray", linestyle="--"); axes[1].plot(bounds, np.array(bounds)/2, color="gray", linestyle="--", label="factor 2")
    for axis in axes: axis.set_xscale("log"); axis.set_yscale("log"); axis.grid(True, alpha=0.25); axis.legend(fontsize=8)
    axes[0].set(xlabel=selected, ylabel=r"observed $\varepsilon_A$", title="Calibration fit; N=6,10 external holdout")
    axes[1].set(xlabel=r"predicted $\varepsilon_A$", ylabel=r"observed $\varepsilon_A$", title="Observed vs predicted; holdout excluded from fit")
    fig.savefig(FIGURES / "t08_transferability.png", dpi=220, metadata={"Software": "acoustic_ms T08"})
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = _read_cases()
    fits, thresholds, selected, fit_map = _analysis(rows)
    _write(DATA / "t08_predictor_fits.csv", fits)
    _write(DATA / "t08_validity_thresholds.csv", thresholds)
    _plot_predictors(rows, fit_map)
    _plot_transferability(rows, selected, fit_map[(selected, "epsilon_a")])
    print(f"T08 selected {selected}; criterion_supported={thresholds[0]['criterion_supported']}")


if __name__ == "__main__":
    main()
