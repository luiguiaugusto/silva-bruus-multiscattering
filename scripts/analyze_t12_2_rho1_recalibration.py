#!/usr/bin/env python3
"""Run the pre-registered T12.2 rho1 recalibration without new force solves."""

from __future__ import annotations

import csv
from pathlib import Path
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    classify_logo_safety,
    confirmatory_metrics,
    evaluate_recalibration_gate,
    fit_log_linear,
    grouped_bootstrap_calibration,
    logo_power_law_predictions,
    power_law_threshold,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
RESOLVED_INPUT = DATA / "t12_1_resolved_comparison.csv"
T12_CONVERGENCE = DATA / "t12_model_e_convergence.csv"
T12_1_CONVERGENCE = DATA / "t12_1_extended_convergence.csv"
PREDICTIONS_PATH = DATA / "t12_2_logo_predictions.csv"
FITS_PATH = DATA / "t12_2_logo_fits.csv"
METRICS_PATH = DATA / "t12_2_metrics.csv"
SAFETY_PATH = DATA / "t12_2_safety_audit.csv"
FINAL_PATH = DATA / "t12_2_final_calibration.csv"
GATE_PATH = DATA / "t12_2_gate.csv"
FIGURE_PATH = FIGURES / "t12_2_rho1_recalibration.png"

FROZEN_PREFACTOR = 2.6353684041458636
FROZEN_EXPONENT = 1.1088518115798773
TOLERANCES = (0.01, 0.05, 0.10)
EXPECTED_GROUPS = (
    "n2_pair",
    "n3_compact",
    "n3_irregular",
    "n3_linear",
    "n4_compact",
    "n4_irregular",
    "n4_linear",
)
INTERNAL_UNCONFIRMED = {
    "n2_pair_f1.0_d2.1",
    "n3_irregular_f1.0_d2.1",
}
BOOTSTRAP_SEED = 1202
BOOTSTRAP_VALID_SAMPLES = 10_000
EPSILON_FLOOR = 0.0


def _format(value: float) -> str:
    return format(float(value), ".17g")


def _bool(value: bool) -> str:
    return str(bool(value)).lower()


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot publish empty table {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _final_rows(path: Path) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read(path):
        grouped.setdefault(row["case_id"], []).append(row)
    return {case_id: rows[-1] for case_id, rows in grouped.items()}


def _load_confirmed_sentinels() -> list[dict[str, str]]:
    rows = _read(RESOLVED_INPUT)
    if len(rows) != 28 or len({row["case_id"] for row in rows}) != 28:
        raise RuntimeError("T12.2 requires exactly 28 unique T12 sentinels")
    if [int(row["sentinel_order"]) for row in rows] != list(range(1, 29)):
        raise RuntimeError("the canonical sentinel order changed")
    if {int(row["particle_count"]) for row in rows} != {2, 3, 4}:
        raise RuntimeError("T12.2 is restricted to N=2,3,4")
    groups = tuple(dict.fromkeys(row["stratum"] for row in rows))
    if groups != EXPECTED_GROUPS:
        raise RuntimeError("the seven frozen (N,family) groups changed")
    if any(row["interaction_confirmed"] != "true" for row in rows):
        raise RuntimeError("all 28 interaction-force references must be confirmed")
    if any(float(row["epsilon_a_e"]) <= EPSILON_FLOOR for row in rows):
        raise RuntimeError("all target errors must be strictly positive")

    old = _final_rows(T12_CONVERGENCE)
    extended = _final_rows(T12_1_CONVERGENCE)
    internal_unconfirmed: set[str] = set()
    for row in rows:
        final = extended.get(row["case_id"], old[row["case_id"]])
        for channel in ("total", "interaction", "external_scattered", "scattered_scattered"):
            row[f"{channel}_confirmed"] = final[f"{channel}_confirmed"]
        if row["total_confirmed"] != "true" or row["interaction_confirmed"] != "true":
            raise RuntimeError("a total/interaction reference is not confirmed")
        if row["scattered_scattered_confirmed"] != "true":
            internal_unconfirmed.add(row["case_id"])
    if internal_unconfirmed != INTERNAL_UNCONFIRMED:
        raise RuntimeError("the two internal unconfirmed statuses changed")
    return rows


def _metric_row(
    model: str,
    scope: str,
    observed: np.ndarray,
    predicted: np.ndarray,
) -> dict[str, object]:
    metrics = confirmatory_metrics(observed, predicted)
    return {
        "record_type": "global_oof" if scope == "all" else "group_oof",
        "model": model,
        "scope": scope,
        "point_count": metrics.point_count,
        "rmse_log": _format(metrics.rmse_log),
        "mae_log": _format(metrics.mae_log),
        "median_absolute_log_ratio": _format(metrics.median_absolute_log_ratio),
        "fraction_within_factor_2": _format(metrics.fraction_within_factor_2),
        "fraction_within_factor_1_5": _format(metrics.fraction_within_factor_1_5),
        "spearman": _format(metrics.spearman),
        "maximum_log_underestimation": _format(metrics.maximum_log_underestimation),
        "evidence_role": "confirmatory_oof" if scope == "all" else "diagnostic_oof",
    }


def _plot(
    sentinels: list[dict[str, str]],
    predictions: list[dict[str, object]],
    safety_rows: list[dict[str, object]],
    final_prefactor: float,
    final_exponent: float,
) -> None:
    colors = {
        "n2_pair": "#4c78a8",
        "n3_compact": "#54a24b",
        "n3_irregular": "#b279a2",
        "n3_linear": "#f58518",
        "n4_compact": "#72b7b2",
        "n4_irregular": "#e45756",
        "n4_linear": "#ff9da6",
    }
    markers = {2: "o", 3: "s", 4: "^"}
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.8), constrained_layout=True)
    for row in predictions:
        group = str(row["held_out_group"])
        particle_count = int(row["particle_count"])
        special = row["scattered_scattered_confirmed"] != "true"
        for field, label, face in (
            ("p0_prediction", "P0 frozen", "none"),
            ("candidate_oof_prediction", "recalibrated OOF", colors[group]),
        ):
            axes[0, 0].scatter(
                float(row[field]),
                float(row["observed_epsilon_a_e"]),
                marker="*" if special else markers[particle_count],
                facecolors=face,
                edgecolors=colors[group],
                s=75 if special else 43,
                linewidths=1.0,
                alpha=0.9,
                label=label if row is predictions[0] else None,
            )
    limits = [
        float(row[field]) for row in predictions
        for field in ("p0_prediction", "candidate_oof_prediction", "observed_epsilon_a_e")
    ]
    grid = np.geomspace(min(limits) * 0.7, max(limits) * 1.4, 300)
    axes[0, 0].plot(grid, grid, color="black", lw=1.1, label="identity")
    axes[0, 0].fill_between(grid, grid / 2.0, grid * 2.0, color="0.88", alpha=0.7, label="factor 2")
    axes[0, 0].set(xscale="log", yscale="log", xlabel="OOF prediction", ylabel=r"observed $\varepsilon_A^E$", title="Strict LOGO predictions")
    axes[0, 0].legend(fontsize=8)

    rho = np.asarray([float(row["rho_l1"]) for row in sentinels])
    observed = np.asarray([float(row["epsilon_a_e"]) for row in sentinels])
    rho_grid = np.geomspace(min(rho) * 0.75, max(rho) * 1.3, 400)
    for row in sentinels:
        special = row["case_id"] in INTERNAL_UNCONFIRMED
        axes[0, 1].scatter(
            float(row["rho_l1"]), float(row["epsilon_a_e"]),
            marker="*" if special else markers[int(row["particle_count"])],
            color=colors[row["stratum"]], s=75 if special else 43,
        )
    axes[0, 1].plot(rho_grid, FROZEN_PREFACTOR * rho_grid**FROZEN_EXPONENT, "k--", lw=1.2, label="P0 frozen")
    axes[0, 1].plot(rho_grid, final_prefactor * rho_grid**final_exponent, "k-", lw=1.5, label="candidate fit on 28")
    axes[0, 1].set(xscale="log", yscale="log", xlabel=r"$\rho_1$", ylabel=r"observed $\varepsilon_A^E$", title="Candidate calibration (not external validation)")
    axes[0, 1].legend(fontsize=8)

    group_index = {group: index for index, group in enumerate(EXPECTED_GROUPS)}
    for row in predictions:
        index = group_index[str(row["held_out_group"])]
        axes[1, 0].scatter(index - 0.10, float(row["p0_log_residual"]), marker="o", facecolors="none", edgecolors="#e45756", s=34)
        axes[1, 0].scatter(index + 0.10, float(row["candidate_log_residual"]), marker="s", color="#4c78a8", s=30)
    axes[1, 0].axhline(0.0, color="black", lw=0.9)
    axes[1, 0].set(xticks=range(len(EXPECTED_GROUPS)), xticklabels=EXPECTED_GROUPS, ylabel="log(predicted / observed)", title="OOF residuals by held-out group")
    axes[1, 0].tick_params(axis="x", rotation=32)
    axes[1, 0].scatter([], [], marker="o", facecolors="none", edgecolors="#e45756", label="P0")
    axes[1, 0].scatter([], [], marker="s", color="#4c78a8", label="recalibrated")
    axes[1, 0].legend(fontsize=8)

    tolerances = [float(row["tolerance"]) for row in safety_rows]
    width = 0.24
    positions = np.arange(len(tolerances))
    axes[1, 1].bar(positions - width, [int(row["predicted_safe_count"]) for row in safety_rows], width, color="#4c78a8", label="predicted safe")
    axes[1, 1].bar(positions, [int(row["true_safe_count"]) for row in safety_rows], width, color="#54a24b", label="true safe")
    axes[1, 1].bar(positions + width, [int(row["false_safe_count"]) for row in safety_rows], width, color="#e45756", label="false safe")
    for index, row in enumerate(safety_rows):
        axes[1, 1].text(index + width, int(row["false_safe_count"]) + 0.12, str(row["false_safe_count"]), ha="center", fontsize=8)
    axes[1, 1].set(xticks=positions, xticklabels=[f"{100*tolerance:.0f}%" for tolerance in tolerances], ylabel="OOF case count", title="Fold-derived safety-threshold audit")
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(True, alpha=0.2)
    figure.savefig(FIGURE_PATH, dpi=220)
    plt.close(figure)


def main() -> None:
    sentinels = _load_confirmed_sentinels()
    case_ids = [row["case_id"] for row in sentinels]
    groups = [row["stratum"] for row in sentinels]
    rho = np.asarray([float(row["rho_l1"]) for row in sentinels])
    observed = np.asarray([float(row["epsilon_a_e"]) for row in sentinels])
    fits, logo_predictions = logo_power_law_predictions(case_ids, rho, observed, groups)
    prediction_by_case = {item.case_id: item for item in logo_predictions}
    if set(prediction_by_case) != set(case_ids):
        raise RuntimeError("each sentinel must receive exactly one OOF prediction")
    classifications, safety_audits = classify_logo_safety(logo_predictions, fits, TOLERANCES)
    class_by_key = {(item.case_id, item.tolerance): item for item in classifications}

    p0 = FROZEN_PREFACTOR * rho**FROZEN_EXPONENT
    candidate = np.asarray([prediction_by_case[case_id].predicted for case_id in case_ids])
    p0_metrics = confirmatory_metrics(observed, p0)
    candidate_metrics = confirmatory_metrics(observed, candidate)
    criteria, decision = evaluate_recalibration_gate(
        candidate_metrics,
        p0_metrics,
        fits,
        safety_audits,
        predictions_finite_positive=bool(np.all(np.isfinite(candidate)) and np.all(candidate > 0.0)),
        integrity_passed=True,
    )

    prediction_rows: list[dict[str, object]] = []
    for row, p0_value, candidate_value in zip(sentinels, p0, candidate):
        result: dict[str, object] = {
            "sentinel_order": row["sentinel_order"],
            "case_id": row["case_id"],
            "particle_count": row["particle_count"],
            "family": row["family"],
            "held_out_group": row["stratum"],
            "f1": row["f1"],
            "distance_ratio": row["distance_ratio"],
            "rho_l1": row["rho_l1"],
            "observed_epsilon_a_e": row["epsilon_a_e"],
            "p0_prediction": _format(p0_value),
            "candidate_oof_prediction": _format(candidate_value),
            "p0_ratio_predicted_over_observed": _format(p0_value / float(row["epsilon_a_e"])),
            "candidate_ratio_predicted_over_observed": _format(candidate_value / float(row["epsilon_a_e"])),
            "p0_log_residual": _format(np.log(p0_value) - np.log(float(row["epsilon_a_e"]))),
            "candidate_log_residual": _format(np.log(candidate_value) - np.log(float(row["epsilon_a_e"]))),
            "total_confirmed": row["total_confirmed"],
            "interaction_confirmed": row["interaction_confirmed"],
            "external_scattered_confirmed": row["external_scattered_confirmed"],
            "scattered_scattered_confirmed": row["scattered_scattered_confirmed"],
            "epsilon_floor": _format(EPSILON_FLOOR),
        }
        for tolerance in TOLERANCES:
            item = class_by_key[(row["case_id"], tolerance)]
            suffix = f"{int(100 * tolerance)}pct"
            result[f"fold_threshold_{suffix}"] = _format(item.threshold)
            result[f"predicted_safe_{suffix}"] = _bool(item.predicted_safe)
            result[f"observed_safe_{suffix}"] = _bool(item.observed_safe)
            result[f"false_safe_{suffix}"] = _bool(item.false_safe)
            result[f"false_unsafe_{suffix}"] = _bool(item.false_unsafe)
        prediction_rows.append(result)

    fit_rows = []
    for fit in fits:
        fit_rows.append({
            "held_out_group": fit.held_out_group,
            "training_count": fit.training_count,
            "test_count": fit.test_count,
            "prefactor": _format(fit.prefactor),
            "exponent": _format(fit.exponent),
            "threshold_1pct": _format(power_law_threshold(0.01, fit.prefactor, fit.exponent)),
            "threshold_5pct": _format(power_law_threshold(0.05, fit.prefactor, fit.exponent)),
            "threshold_10pct": _format(power_law_threshold(0.10, fit.prefactor, fit.exponent)),
            "coefficients_positive": _bool(fit.prefactor > 0.0 and fit.exponent > 0.0),
        })

    metric_rows = [
        _metric_row("P0_frozen_rho1", "all", observed, p0),
        _metric_row("recalibrated_rho1", "all", observed, candidate),
    ]
    for group in EXPECTED_GROUPS:
        selected = np.asarray([value == group for value in groups])
        metric_rows.append(_metric_row("P0_frozen_rho1", group, observed[selected], p0[selected]))
        metric_rows.append(_metric_row("recalibrated_rho1", group, observed[selected], candidate[selected]))

    safety_rows = [{
        "tolerance": _format(item.tolerance),
        "predicted_safe_count": item.predicted_safe_count,
        "predicted_safe_group_count": item.predicted_safe_group_count,
        "true_safe_count": item.true_safe_count,
        "false_safe_count": item.false_safe_count,
        "false_unsafe_count": item.false_unsafe_count,
        "worst_false_safe_excess": _format(item.worst_false_safe_excess),
        "coverage_sufficient": _bool(item.coverage_sufficient),
        "zero_false_safe": _bool(item.false_safe_count == 0),
    } for item in safety_audits]

    final_fit = fit_log_linear(rho, observed)
    bootstrap = grouped_bootstrap_calibration(
        rho,
        observed,
        groups,
        tolerances=TOLERANCES,
        seed=BOOTSTRAP_SEED,
        valid_samples=BOOTSTRAP_VALID_SAMPLES,
    )
    final_thresholds = [
        power_law_threshold(value, final_fit.prefactor, final_fit.coefficient)
        for value in TOLERANCES
    ]
    final_rows = [{
        "calibration_role": "candidate_frozen_for_external_validation",
        "point_count": len(sentinels),
        "group_count": len(EXPECTED_GROUPS),
        "prefactor": _format(final_fit.prefactor),
        "exponent": _format(final_fit.coefficient),
        "threshold_1pct": _format(final_thresholds[0]),
        "threshold_5pct": _format(final_thresholds[1]),
        "threshold_10pct": _format(final_thresholds[2]),
        "fold_prefactor_minimum": _format(min(item.prefactor for item in fits)),
        "fold_prefactor_maximum": _format(max(item.prefactor for item in fits)),
        "fold_exponent_minimum": _format(min(item.exponent for item in fits)),
        "fold_exponent_maximum": _format(max(item.exponent for item in fits)),
        "bootstrap_seed": bootstrap.seed,
        "bootstrap_valid_samples": bootstrap.valid_samples,
        "bootstrap_attempts": bootstrap.attempts,
        "bootstrap_prefactor_low": _format(bootstrap.prefactor_interval[0]),
        "bootstrap_prefactor_high": _format(bootstrap.prefactor_interval[1]),
        "bootstrap_exponent_low": _format(bootstrap.exponent_interval[0]),
        "bootstrap_exponent_high": _format(bootstrap.exponent_interval[1]),
        "bootstrap_threshold_1pct_low": _format(bootstrap.threshold_intervals[0][0]),
        "bootstrap_threshold_1pct_high": _format(bootstrap.threshold_intervals[0][1]),
        "bootstrap_threshold_5pct_low": _format(bootstrap.threshold_intervals[1][0]),
        "bootstrap_threshold_5pct_high": _format(bootstrap.threshold_intervals[1][1]),
        "bootstrap_threshold_10pct_low": _format(bootstrap.threshold_intervals[2][0]),
        "bootstrap_threshold_10pct_high": _format(bootstrap.threshold_intervals[2][1]),
        "p0_prefactor": _format(FROZEN_PREFACTOR),
        "p0_exponent": _format(FROZEN_EXPONENT),
        "p0_threshold_1pct": _format(power_law_threshold(0.01, FROZEN_PREFACTOR, FROZEN_EXPONENT)),
        "p0_threshold_5pct": _format(power_law_threshold(0.05, FROZEN_PREFACTOR, FROZEN_EXPONENT)),
        "p0_threshold_10pct": _format(power_law_threshold(0.10, FROZEN_PREFACTOR, FROZEN_EXPONENT)),
        "epsilon_floor": _format(EPSILON_FLOOR),
    }]

    gate_rows = [{
        "record_type": "criterion",
        "criterion": item.name,
        "observed": _format(item.observed),
        "threshold": _format(item.threshold),
        "passed": _bool(item.passed),
        "justification": item.justification,
        "decision": "pending",
    } for item in criteria]
    gate_rows.append({
        "record_type": "decision",
        "criterion": "final_decision",
        "observed": "0",
        "threshold": "0",
        "passed": _bool(decision == "GO_T13_WITH_RECALIBRATED_RHO1"),
        "justification": "all ten pre-registered criteria must pass",
        "decision": decision,
    })

    _atomic_write(PREDICTIONS_PATH, prediction_rows)
    _atomic_write(FITS_PATH, fit_rows)
    _atomic_write(METRICS_PATH, metric_rows)
    _atomic_write(SAFETY_PATH, safety_rows)
    _atomic_write(FINAL_PATH, final_rows)
    _atomic_write(GATE_PATH, gate_rows)
    _plot(sentinels, prediction_rows, safety_rows, final_fit.prefactor, final_fit.coefficient)
    print(f"T12.2 decision: {decision}")


if __name__ == "__main__":
    main()
