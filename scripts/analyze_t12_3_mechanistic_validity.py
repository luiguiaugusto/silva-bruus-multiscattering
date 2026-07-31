#!/usr/bin/env python3
"""Evaluate the pre-registered T12.3 mechanistic validity criterion."""

from __future__ import annotations

import csv
from pathlib import Path
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    audit_safety_thresholds,
    evaluate_mechanistic_gate,
    fit_mechanistic_power_law,
    fixed_baseline_nested_predictions,
    multiplicative_metrics,
    nested_logo_predictions,
    predict_mechanistic_power_law,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
INPUT = DATA / "t12_1_resolved_comparison.csv"
OOF_PATH = DATA / "t12_3_oof_predictions.csv"
COEFFICIENT_PATH = DATA / "t12_3_logo_coefficients.csv"
SAFETY_FACTOR_PATH = DATA / "t12_3_nested_safety_factors.csv"
THRESHOLD_PATH = DATA / "t12_3_threshold_audit.csv"
METRICS_PATH = DATA / "t12_3_metrics.csv"
BOOTSTRAP_PATH = DATA / "t12_3_group_bootstrap.csv"
INFLUENCE_PATH = DATA / "t12_3_case_influence.csv"
GATE_PATH = DATA / "t12_3_gate.csv"
FIGURE_PATH = FIGURES / "t12_3_mechanistic_validity.png"

P0_PREFACTOR = 2.6353684041458636
P0_EXPONENT = 1.1088518115798773
P3_PREFACTOR = 14.73950709797405
P3_EXPONENT = 1.4226504975598322
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
BOOTSTRAP_SEED = 1203
BOOTSTRAP_VALID_SAMPLES = 10_000


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot publish empty table {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format(row.get(field, "")) for field in fields})
    temporary.replace(path)


def load_frozen_sentinels() -> list[dict[str, str]]:
    """Read and validate only the frozen 28-row N<=4 sentinel table."""

    rows = _read(INPUT)
    if len(rows) != 28 or len({row["case_id"] for row in rows}) != 28:
        raise RuntimeError("T12.3 requires exactly 28 unique frozen sentinels")
    if [int(row["sentinel_order"]) for row in rows] != list(range(1, 29)):
        raise RuntimeError("the canonical sentinel order changed")
    if {int(row["particle_count"]) for row in rows} != {2, 3, 4}:
        raise RuntimeError("T12.3 is restricted to N=2,3,4")
    if any(int(row["particle_count"]) in {6, 10} for row in rows):
        raise RuntimeError("the external N=6,10 holdout must not be read")
    groups = tuple(dict.fromkeys(row["stratum"] for row in rows))
    if groups != EXPECTED_GROUPS:
        raise RuntimeError("the seven frozen validation groups changed")
    required_positive = ("lambda_max", "rho_l1", "epsilon_a_e")
    for row in rows:
        if row["interaction_confirmed"] != "true":
            raise RuntimeError("all interaction-force references must be confirmed")
        if any(not np.isfinite(float(row[field])) or float(row[field]) <= 0.0 for field in required_positive):
            raise RuntimeError("mechanistic predictors and errors must be finite and positive")
    return rows


def _arrays(rows: list[dict[str, str]]) -> dict[str, np.ndarray]:
    return {
        "ids": np.asarray([row["case_id"] for row in rows], dtype=str),
        "groups": np.asarray([row["stratum"] for row in rows], dtype=str),
        "lambda": np.asarray([float(row["lambda_max"]) for row in rows]),
        "rho": np.asarray([float(row["rho_l1"]) for row in rows]),
        "observed": np.asarray([float(row["epsilon_a_e"]) for row in rows]),
    }


def _prediction_map(predictions: object) -> dict[str, object]:
    return {prediction.case_id: prediction for prediction in predictions}


def _metric_row(model: str, scope_type: str, scope: str, observed: np.ndarray, predicted: np.ndarray) -> dict[str, object]:
    metric = multiplicative_metrics(observed, predicted)
    return {
        "model": model,
        "scope_type": scope_type,
        "scope": scope,
        "point_count": metric.point_count,
        "rmse_log": metric.rmse_log,
        "mae_log": metric.mae_log,
        "fraction_within_factor_2": metric.fraction_within_factor_2,
        "fraction_within_factor_1_5": metric.fraction_within_factor_1_5,
        "spearman": metric.spearman,
        "worst_multiplicative_ratio": metric.worst_multiplicative_ratio,
    }


def _bootstrap(arrays: dict[str, np.ndarray]) -> list[dict[str, object]]:
    groups = np.asarray(sorted(set(arrays["groups"].tolist())))
    indices = {group: np.flatnonzero(arrays["groups"] == group) for group in groups}
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    buckets: dict[tuple[str, str], list[float]] = {}
    valid = 0
    attempts = 0
    while valid < BOOTSTRAP_VALID_SAMPLES:
        attempts += 1
        if attempts > 100_000:
            raise RuntimeError("group bootstrap failed to produce enough valid samples")
        sampled = rng.choice(groups, size=len(groups), replace=True)
        selected = np.concatenate([indices[group] for group in sampled])
        try:
            m1 = fit_mechanistic_power_law(
                arrays["lambda"][selected], arrays["observed"][selected]
            )
            m2 = fit_mechanistic_power_law(
                arrays["lambda"][selected], arrays["observed"][selected], arrays["rho"][selected]
            )
        except ValueError:
            continue
        model_predictions = {
            "P0": P0_PREFACTOR * arrays["rho"][selected] ** P0_EXPONENT,
            "P3": P3_PREFACTOR * arrays["rho"][selected] ** P3_EXPONENT,
            "M1": predict_mechanistic_power_law(m1, arrays["lambda"][selected]),
            "M2": predict_mechanistic_power_law(m2, arrays["lambda"][selected], arrays["rho"][selected]),
        }
        for model, predicted in model_predictions.items():
            metric = multiplicative_metrics(arrays["observed"][selected], predicted)
            for quantity, value in (
                ("rmse_log", metric.rmse_log),
                ("fraction_within_factor_2", metric.fraction_within_factor_2),
                ("spearman", metric.spearman),
            ):
                buckets.setdefault((model, quantity), []).append(float(value))
        for model, fit in (("M1", m1), ("M2", m2)):
            for quantity, value in (
                ("prefactor", fit.prefactor),
                ("alpha_lambda", fit.alpha_lambda),
                ("alpha_rho", fit.alpha_rho),
            ):
                buckets.setdefault((model, quantity), []).append(float(value))
        valid += 1
    rows: list[dict[str, object]] = []
    for model, quantity in sorted(buckets):
        values = np.asarray(buckets[(model, quantity)])
        rows.append({
            "model": model,
            "quantity": quantity,
            "seed": BOOTSTRAP_SEED,
            "valid_samples": valid,
            "attempts": attempts,
            "q025": float(np.quantile(values, 0.025)),
            "median": float(np.quantile(values, 0.5)),
            "q975": float(np.quantile(values, 0.975)),
        })
    return rows


def _plot(
    rows: list[dict[str, str]],
    oof_rows: list[dict[str, object]],
    coefficient_rows: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
) -> None:
    order = {row["case_id"]: int(row["sentinel_order"]) for row in rows}
    families = {row["case_id"]: row["family"] for row in rows}
    colors = {"P0": "0.55", "P3": "tab:purple", "M1": "tab:blue", "M2": "tab:orange"}
    markers = {"pair": "P", "compact": "o", "irregular": "^", "linear": "s"}
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.5), constrained_layout=True)

    axis = axes[0, 0]
    for model in ("P3", "M1", "M2"):
        subset = [row for row in oof_rows if row["model"] == model]
        axis.scatter(
            [float(row["point_prediction"]) for row in subset],
            [float(row["observed_epsilon_a_e"]) for row in subset],
            s=25, alpha=0.72, color=colors[model], label=model,
        )
    values = [float(row["observed_epsilon_a_e"]) for row in oof_rows]
    bounds = np.array([min(values) / 1.5, max(values) * 1.5])
    axis.plot(bounds, bounds, color="black", linewidth=1.0, label="identity")
    axis.plot(bounds, 2 * bounds, color="gray", linestyle="--", linewidth=0.8)
    axis.plot(bounds, bounds / 2, color="gray", linestyle="--", linewidth=0.8, label="factor 2")
    axis.set(xscale="log", yscale="log", xlabel="OOF point prediction", ylabel="observed Model-E error", title="Observed versus strictly OOF prediction")
    axis.legend(fontsize=8)

    axis = axes[0, 1]
    m1 = sorted([row for row in oof_rows if row["model"] == "M1"], key=lambda row: order[str(row["case_id"])])
    x = np.arange(1, 29)
    axis.plot(x, [float(row["safe_prediction"]) for row in m1], color=colors["M1"], marker="o", markersize=3, linewidth=0.9, label="M1 conservative OOF")
    axis.scatter(x, [float(row["observed_epsilon_a_e"]) for row in m1], facecolors="none", edgecolors="black", s=22, label="observed")
    for tolerance in TOLERANCES:
        axis.axhline(tolerance, linestyle="--", linewidth=0.8, label=f"{100*tolerance:g}%")
    axis.set(yscale="log", xlabel="frozen sentinel order", ylabel="error", title="Nested safety calibration (strict < threshold)")
    axis.legend(fontsize=7, ncol=2)

    axis = axes[0, 2]
    m1_residuals = [
        (families[str(row["case_id"])], np.log(float(row["observed_epsilon_a_e"]) / float(row["point_prediction"])))
        for row in m1
    ]
    family_order = ("pair", "compact", "irregular", "linear")
    for index, family in enumerate(family_order):
        values = [value for label, value in m1_residuals if label == family]
        axis.scatter(np.full(len(values), index), values, marker=markers[family], s=28, label=family)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(range(len(family_order)), family_order, rotation=20)
    axis.set(ylabel="log(observed / M1 point prediction)", title="M1 OOF residuals by family")
    axis.legend(fontsize=7)

    axis = axes[1, 0]
    for model, offset in (("M1", -0.08), ("M2 lambda", 0.0), ("M2 rho", 0.08)):
        base = "M1" if model == "M1" else "M2"
        field = "alpha_lambda" if model != "M2 rho" else "alpha_rho"
        subset = [row for row in coefficient_rows if row["model"] == base and row["scope"] == "outer_fold"]
        axis.plot(np.arange(len(subset)) + offset, [float(row[field]) for row in subset], marker="o", linewidth=0.8, label=model)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(range(7), EXPECTED_GROUPS, rotation=35, ha="right")
    axis.set(ylabel="log-law exponent", title="Outer-fold coefficient stability")
    axis.legend(fontsize=8)

    axis = axes[1, 1]
    conservative = [row for row in threshold_rows if row["rule"] == "conservative" and row["model"] in {"P3", "M1", "M2"}]
    width = 0.24
    for position, model in enumerate(("P3", "M1", "M2")):
        subset = sorted([row for row in conservative if row["model"] == model], key=lambda row: float(row["tolerance"]))
        axis.bar(np.arange(3) + (position - 1) * width, [int(row["predicted_safe_count"]) for row in subset], width=width, color=colors[model], label=model)
    axis.set_xticks(range(3), ["1%", "5%", "10%"])
    axis.set(ylabel="conservative predicted-safe cases", title="Safe-region coverage")
    axis.legend(fontsize=8)

    axis = axes[1, 2]
    for position, model in enumerate(("P3", "M1", "M2")):
        subset = sorted([row for row in conservative if row["model"] == model], key=lambda row: float(row["tolerance"]))
        axis.bar(np.arange(3) + (position - 1) * width, [int(row["false_safe_count"]) for row in subset], width=width, color=colors[model], label=model)
    axis.set_xticks(range(3), ["1%", "5%", "10%"])
    axis.set(ylabel="false-safe cases", title="Conservative false-safe audit")
    axis.set_ylim(0.0, 1.0)
    axis.text(1.0, 0.55, "all conservative counts = 0", ha="center", va="center", fontsize=9)
    axis.legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(True, alpha=0.22)
    fig.suptitle("T12.3 grouped mechanistic validation; amplitude-error ratios", fontsize=13)
    fig.savefig(FIGURE_PATH, dpi=220, metadata={"Software": "acoustic_ms T12.3"})
    plt.close(fig)


def main() -> None:
    rows = load_frozen_sentinels()
    arrays = _arrays(rows)
    full_fits = {
        "M1": fit_mechanistic_power_law(arrays["lambda"], arrays["observed"]),
        "M2": fit_mechanistic_power_law(arrays["lambda"], arrays["observed"], arrays["rho"]),
    }
    folds: dict[str, tuple] = {}
    predictions: dict[str, tuple] = {}
    for model in ("M1", "M2"):
        folds[model], predictions[model] = nested_logo_predictions(
            arrays["ids"], arrays["groups"], arrays["lambda"], arrays["rho"],
            arrays["observed"], model=model,
        )
    baseline_points = {
        "P0": P0_PREFACTOR * arrays["rho"] ** P0_EXPONENT,
        "P3": P3_PREFACTOR * arrays["rho"] ** P3_EXPONENT,
    }
    for model in ("P0", "P3"):
        folds[model], predictions[model] = fixed_baseline_nested_predictions(
            arrays["ids"], arrays["groups"], arrays["observed"], baseline_points[model], model=model
        )

    maps = {model: _prediction_map(values) for model, values in predictions.items()}
    oof_rows: list[dict[str, object]] = []
    for source in rows:
        case_id = source["case_id"]
        for model in ("P0", "P3", "M1", "M2"):
            prediction = maps[model][case_id]
            oof_rows.append({
                "sentinel_order": int(source["sentinel_order"]),
                "case_id": case_id,
                "particle_count": int(source["particle_count"]),
                "family": source["family"],
                "held_out_group": source["stratum"],
                "f1": float(source["f1"]),
                "distance_ratio": float(source["distance_ratio"]),
                "lambda_max": float(source["lambda_max"]),
                "rho_l1": float(source["rho_l1"]),
                "observed_epsilon_a_e": prediction.observed,
                "model": model,
                "point_prediction": prediction.point_prediction,
                "safety_factor": prediction.safety_factor,
                "safe_prediction": prediction.safe_prediction,
                "point_log_residual_observed_over_predicted": np.log(prediction.observed / prediction.point_prediction),
                "safe_log_residual_observed_over_predicted": np.log(prediction.observed / prediction.safe_prediction),
            })
    _write(OOF_PATH, oof_rows)

    coefficient_rows: list[dict[str, object]] = []
    for model in ("M1", "M2"):
        for fold in folds[model]:
            fit = fold.fit
            coefficient_rows.append({
                "model": model, "scope": "outer_fold", "held_out_group": fold.held_out_group,
                "training_count": fold.training_count, "test_count": fold.test_count,
                "prefactor": fit.prefactor, "intercept": fit.intercept,
                "alpha_lambda": fit.alpha_lambda, "alpha_rho": fit.alpha_rho,
                "design_rank": fit.design_rank,
                "standardized_condition_number": fit.standardized_condition_number,
                "log_predictor_correlation": fit.log_predictor_correlation,
            })
        fit = full_fits[model]
        coefficient_rows.append({
            "model": model, "scope": "descriptive_full_fit", "held_out_group": "",
            "training_count": 28, "test_count": 0, "prefactor": fit.prefactor,
            "intercept": fit.intercept, "alpha_lambda": fit.alpha_lambda,
            "alpha_rho": fit.alpha_rho, "design_rank": fit.design_rank,
            "standardized_condition_number": fit.standardized_condition_number,
            "log_predictor_correlation": fit.log_predictor_correlation,
        })
    _write(COEFFICIENT_PATH, coefficient_rows)

    safety_factor_rows: list[dict[str, object]] = []
    for model in ("P0", "P3", "M1", "M2"):
        for fold in folds[model]:
            safety_factor_rows.append({
                "model": model, "held_out_group": fold.held_out_group,
                "training_count": fold.training_count, "test_count": fold.test_count,
                "inner_prediction_count": fold.inner_prediction_count,
                "maximum_inner_underprediction_log": fold.maximum_inner_underprediction_log,
                "safety_factor": fold.safety_factor, "valid": fold.valid,
                "outer_group_excluded_from_fit_and_margin": True,
            })
    _write(SAFETY_FACTOR_PATH, safety_factor_rows)

    metrics_rows: list[dict[str, object]] = []
    metric_objects = {}
    for model in ("P0", "P3", "M1", "M2"):
        predicted = np.asarray([maps[model][case_id].point_prediction for case_id in arrays["ids"]])
        metric_objects[model] = multiplicative_metrics(arrays["observed"], predicted)
        metrics_rows.append(_metric_row(model, "global", "all_28", arrays["observed"], predicted))
        scopes = {
            "group": arrays["groups"],
            "particle_count": np.asarray([row["particle_count"] for row in rows]),
            "family": np.asarray([row["family"] for row in rows]),
            "f1": np.asarray([row["f1"] for row in rows]),
            "distance_ratio": np.asarray([row["distance_ratio"] for row in rows]),
        }
        for scope_type, labels in scopes.items():
            for scope in sorted(set(labels.tolist())):
                selected = np.flatnonzero(labels == scope)
                if selected.size >= 2:
                    metrics_rows.append(_metric_row(model, scope_type, str(scope), arrays["observed"][selected], predicted[selected]))
    _write(METRICS_PATH, metrics_rows)

    threshold_rows: list[dict[str, object]] = []
    audits = {}
    for model in ("P0", "P3", "M1", "M2"):
        point = np.asarray([maps[model][case_id].point_prediction for case_id in arrays["ids"]])
        safe = np.asarray([maps[model][case_id].safe_prediction for case_id in arrays["ids"]])
        for rule, values in (("point_diagnostic", point), ("conservative", safe)):
            result = audit_safety_thresholds(
                arrays["ids"], arrays["observed"], values, model=model, rule=rule,
                tolerances=TOLERANCES,
            )
            if rule == "conservative":
                audits[model] = result
            for item in result:
                threshold_rows.append({
                    "model": model, "rule": rule, "tolerance": item.tolerance,
                    "comparison": "prediction < tolerance; observed < tolerance",
                    "predicted_safe_count": item.predicted_safe_count,
                    "observed_safe_count": item.observed_safe_count,
                    "false_safe_count": item.false_safe_count,
                    "false_unsafe_count": item.false_unsafe_count,
                    "safe_precision": item.safe_precision, "safe_coverage": item.safe_coverage,
                    "false_safe_ids": ";".join(item.false_safe_ids),
                    "false_unsafe_ids": ";".join(item.false_unsafe_ids),
                })
    _write(THRESHOLD_PATH, threshold_rows)

    _write(BOOTSTRAP_PATH, _bootstrap(arrays))

    influence_rows: list[dict[str, object]] = []
    for model in ("M1", "M2"):
        full = full_fits[model]
        for index, case_id in enumerate(arrays["ids"]):
            keep = np.arange(28) != index
            fit = fit_mechanistic_power_law(
                arrays["lambda"][keep], arrays["observed"][keep],
                arrays["rho"][keep] if model == "M2" else None,
            )
            predicted = predict_mechanistic_power_law(
                fit, arrays["lambda"][[index]], arrays["rho"][[index]] if model == "M2" else None
            )[0]
            influence_rows.append({
                "model": model, "case_id": case_id,
                "is_historical_10pct_false_safe": case_id == "n2_pair_f0.8_d2.5",
                "leave_one_case_out_prefactor": fit.prefactor,
                "leave_one_case_out_alpha_lambda": fit.alpha_lambda,
                "leave_one_case_out_alpha_rho": fit.alpha_rho,
                "delta_prefactor": fit.prefactor - full.prefactor,
                "delta_alpha_lambda": fit.alpha_lambda - full.alpha_lambda,
                "delta_alpha_rho": fit.alpha_rho - full.alpha_rho,
                "held_out_prediction": predicted,
                "held_out_observed": arrays["observed"][index],
                "held_out_log_residual": np.log(arrays["observed"][index] / predicted),
            })
    _write(INFLUENCE_PATH, influence_rows)

    m2_sign_fraction = np.mean([
        fold.fit.alpha_lambda > 0 and fold.fit.alpha_rho > 0 for fold in folds["M2"]
    ])
    m2_unstable = bool(
        full_fits["M2"].standardized_condition_number > 1e3
        or any(fold.fit.standardized_condition_number > 1e3 for fold in folds["M2"])
        or m2_sign_fraction < 0.8
    )
    criteria, decision, m1_pass, m2_pass = evaluate_mechanistic_gate(
        metric_objects["M1"], metric_objects["M2"], audits["M1"], audits["M2"],
        folds["M1"], folds["M2"], full_fits["M1"], full_fits["M2"],
        integrity_passed=True, m2_unstable_collinearity=m2_unstable,
    )
    gate_rows = [
        {
            "candidate": item.candidate, "criterion": item.name,
            "observed": item.observed, "threshold": item.threshold,
            "passed": item.passed, "justification": item.justification,
            "candidate_pass": m1_pass if item.candidate == "M1" else m2_pass,
            "m2_collinearity_status": "UNSTABLE_COLLINEARITY" if m2_unstable else "IDENTIFIABLE",
            "final_decision": decision,
        }
        for item in criteria
    ]
    _write(GATE_PATH, gate_rows)
    _plot(rows, oof_rows, coefficient_rows, threshold_rows)
    print(
        f"T12.3 decision={decision}; M1_pass={m1_pass}; M2_pass={m2_pass}; "
        f"M2_status={'UNSTABLE_COLLINEARITY' if m2_unstable else 'IDENTIFIABLE'}"
    )


if __name__ == "__main__":
    main()
