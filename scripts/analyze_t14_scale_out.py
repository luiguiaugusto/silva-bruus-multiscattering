#!/usr/bin/env python3
"""Analyze the frozen T14 campaign without importing or calling Model E."""

from __future__ import annotations

from dataclasses import asdict
import csv
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    EXPECTED_SCALE_OUT_CASE_IDS,
    TOLERANCES,
    audit_external_threshold,
    evaluate_scale_out_gate,
    external_eligibility_mask,
    external_prediction_metrics,
    normalized_rms_error_xyz,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
MANIFEST = DATA / "t14_scale_manifest.csv"
FROZEN_PREDICTIONS = DATA / "t14_frozen_predictions.csv"
PROTOCOL = DATA / "t14_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_prior_artifact_hashes.csv"
RAW = DATA / "t14_model_e_convergence.csv"
FORCES = DATA / "t14_forces.csv"
SUMMARY = DATA / "t14_case_summary.csv"
PREDICTIONS = DATA / "t14_scale_predictions.csv"
METRICS = DATA / "t14_metrics.csv"
THRESHOLDS = DATA / "t14_threshold_audit.csv"
PAIRS = DATA / "t14_matched_scale_pairs.csv"
PERFORMANCE = DATA / "t14_performance.csv"
GATE = DATA / "t14_gate.csv"
FIGURE = FIGURES / "t14_scale_out_validation.png"
T13_SUMMARY = DATA / "t13_case_summary.csv"
PHASE_A_PATHS = (MANIFEST, FROZEN_PREDICTIONS, PROTOCOL, PRIOR_HASHES)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _truth(value: str) -> bool:
    if value not in {"true", "false"}:
        raise ValueError(f"invalid Boolean: {value}")
    return value == "true"


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot publish empty table {path.name}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format(row.get(field, "")) for field in fields})
    temporary.replace(path)


def _vectors(value: str) -> np.ndarray:
    result = np.asarray([
        [float(item) for item in row.split(":")] for row in value.split(";")
    ])
    if result.ndim != 2 or result.shape[1] != 3 or not np.all(np.isfinite(result)):
        raise ValueError("invalid serialized vectors")
    return result


def _phase_a_integrity() -> bool:
    for path in PHASE_A_PATHS:
        expected = subprocess.run(
            ["git", "show", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT,
            check=True, capture_output=True,
        ).stdout
        if sha256(path.read_bytes()).digest() != sha256(expected).digest():
            return False
    return True


def _prior_integrity() -> bool:
    for row in _read(PRIOR_HASHES):
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            return False
        if sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            return False
    return True


def _validated_raw() -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read(RAW):
        grouped.setdefault(row["case_id"], []).append(row)
    if tuple(grouped) != EXPECTED_SCALE_OUT_CASE_IDS:
        raise RuntimeError("T14 raw case identity differs from preregistration")
    for case_id, rows in grouped.items():
        rows.sort(key=lambda row: int(row["lmax"]))
        orders = [int(row["lmax"]) for row in rows]
        if orders != list(range(2, orders[-1] + 1)) or orders[-1] > 13:
            raise RuntimeError(f"invalid Lmax sequence for {case_id}")
        for row in rows:
            for field, value in row.items():
                if field.endswith("forces_xyz") or field == "coordinates_xyz":
                    continue
                if field in {"case_id", "family", "stop_reason", "production_solver"}:
                    continue
                if value.lower() in {"true", "false"}:
                    continue
                if not np.isfinite(float(value)):
                    raise RuntimeError(f"nonfinite {field} for {case_id}")
    return grouped


def _metric_row(
    model: str,
    scope_type: str,
    scope: str,
    observed: np.ndarray,
    predicted: np.ndarray,
) -> tuple[dict[str, object], object | None]:
    if len(observed) < 2:
        return {
            "model": model, "scope_type": scope_type, "scope": scope,
            "applicable": False, "reason": "fewer_than_two_eligible_cases",
            "point_count": len(observed), "rmse_log": 0.0, "mae_log": 0.0,
            "median_factor": 0.0, "p90_factor": 0.0, "maximum_factor": 0.0,
            "fraction_within_factor_2": 0.0,
            "fraction_within_factor_1_5": 0.0, "spearman": 0.0,
            "mean_log_bias": 0.0, "median_log_bias": 0.0,
        }, None
    metric = external_prediction_metrics(observed, predicted)
    return {
        "model": model, "scope_type": scope_type, "scope": scope,
        "applicable": True, "reason": "applicable", **asdict(metric),
    }, metric


def analyze() -> None:
    if not RAW.exists():
        raise RuntimeError("T14 raw campaign is absent; analysis never calls Model E")
    phase_a_integrity = _phase_a_integrity()
    prior_integrity = _prior_integrity()
    manifest = _read(MANIFEST)
    if len(manifest) != 24 or tuple(row["case_id"] for row in manifest) != EXPECTED_SCALE_OUT_CASE_IDS:
        raise RuntimeError("T14 manifest identity failed")
    grouped = _validated_raw()
    frozen = _read(FROZEN_PREDICTIONS)
    prediction_map = {(row["case_id"], row["model"]): row for row in frozen}
    if len(prediction_map) != 48:
        raise RuntimeError("T14 must contain 48 unique frozen predictions")

    summaries: list[dict[str, object]] = []
    force_rows: list[dict[str, object]] = []
    performance_rows: list[dict[str, object]] = []
    for source in manifest:
        case_id = source["case_id"]
        rows = grouped[case_id]
        final = rows[-1]
        positions = _vectors(final["coordinates_xyz"])
        model_a = _vectors(final["model_a_forces_xyz"])
        total = _vectors(final["total_forces_xyz"])
        interaction = _vectors(final["interaction_forces_xyz"])
        external_scattered = _vectors(final["external_scattered_forces_xyz"])
        scattered_scattered = _vectors(final["scattered_scattered_forces_xyz"])
        epsilon, applicable = normalized_rms_error_xyz(interaction, model_a)
        confirmed = _truth(final["interaction_confirmed"])
        diagnostics = _truth(final["diagnostics_pass"])
        eligible = bool(external_eligibility_mask(
            [confirmed], [diagnostics], [applicable]
        )[0])
        m1 = prediction_map[(case_id, "M1")]
        p3 = prediction_map[(case_id, "P3")]
        ineligibility = "eligible"
        if not confirmed:
            ineligibility = "interaction_unconfirmed_at_13"
        elif not diagnostics:
            ineligibility = "final_numerical_diagnostics_failed"
        elif not applicable:
            ineligibility = "interaction_force_numerically_null"
        summary = {
            "scale_order": int(source["scale_order"]),
            "case_id": case_id,
            "particle_count": int(source["particle_count"]),
            "family": source["family"],
            "target_level": int(source["target_level"]),
            "f1": float(source["f1"]),
            "distance_ratio": float(source["distance_ratio"]),
            "minimum_distance": float(source["minimum_distance"]),
            "cluster_diameter": float(source["cluster_diameter"]),
            "k_cluster_diameter": float(source["k_cluster_diameter"]),
            "lambda_max": float(source["lambda_max"]),
            "rho_l1": float(source["rho_l1"]),
            "final_lmax": int(final["final_lmax"]),
            "stop_reason": final["stop_reason"],
            "interaction_confirmed": confirmed,
            "total_confirmed": _truth(final["total_confirmed"]),
            "external_scattered_confirmed": _truth(final["external_scattered_confirmed"]),
            "scattered_scattered_confirmed": _truth(final["scattered_scattered_confirmed"]),
            "diagnostics_pass": diagnostics,
            "error_applicable": applicable,
            "eligible": eligible,
            "ineligibility_reason": ineligibility,
            "epsilon_a_e": epsilon,
            "m1_point_prediction": float(m1["point_prediction"]),
            "m1_conservative_prediction": float(m1["conservative_prediction"]),
            "p3_point_prediction": float(p3["point_prediction"]),
            "p3_conservative_prediction": float(p3["conservative_prediction"]),
            "m1_log_residual": np.log(float(m1["point_prediction"]) / epsilon) if eligible else 0.0,
            "final_system_dimension": int(final["system_dimension"]),
            "case_wall_seconds": float(final["case_accumulated_seconds"]),
            "peak_memory_kib": max(int(row["process_peak_memory_kib"]) for row in rows),
            "maximum_balanced_condition": max(float(row["balanced_condition_number"]) for row in rows),
            "maximum_backward_error": max(float(row["balanced_backward_error"]) for row in rows),
            "maximum_incident_closure_error": max(float(row["effective_incident_closure_error"]) for row in rows),
            "maximum_scattering_closure_error": max(float(row["scattering_closure_error"]) for row in rows),
            "maximum_force_decomposition_residual": max(float(row["force_decomposition_residual"]) for row in rows),
            "maximum_abs_fz": max(float(row["max_abs_fz"]) for row in rows),
        }
        summaries.append(summary)
        for particle in range(len(positions)):
            output: dict[str, object] = {
                "scale_order": int(source["scale_order"]), "case_id": case_id,
                "particle_index": particle,
            }
            for prefix, values in (
                ("position", positions), ("model_a", model_a),
                ("model_e_total", total), ("model_e_interaction", interaction),
                ("model_e_external_scattered", external_scattered),
                ("model_e_scattered_scattered", scattered_scattered),
            ):
                for axis, value in zip(("x", "y", "z"), values[particle]):
                    output[f"{prefix}_{axis}"] = value
            force_rows.append(output)
        performance_rows.append({
            "scale_order": int(source["scale_order"]), "case_id": case_id,
            "particle_count": int(source["particle_count"]), "family": source["family"],
            "target_level": int(source["target_level"]),
            "orders_computed": len(rows), "final_lmax": int(final["final_lmax"]),
            "final_system_dimension": int(final["system_dimension"]),
            "case_wall_seconds": float(final["case_accumulated_seconds"]),
            "mean_order_wall_seconds": float(np.mean([float(row["order_wall_seconds"]) for row in rows])),
            "maximum_order_wall_seconds": max(float(row["order_wall_seconds"]) for row in rows),
            "peak_memory_kib": max(int(row["process_peak_memory_kib"]) for row in rows),
            "stop_reason": final["stop_reason"],
        })

    revealed = []
    for row in frozen:
        summary = next(item for item in summaries if item["case_id"] == row["case_id"])
        revealed.append({
            **row, "eligible": summary["eligible"],
            "observed_epsilon_a_e": summary["epsilon_a_e"],
            "observed_safe_1pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.01),
            "observed_safe_5pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.05),
            "observed_safe_10pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.10),
            "response_source": "Model_E_interaction_force_only",
        })

    scopes: list[tuple[str, str, list[int]]] = [("global", "all", list(range(24)))]
    for n in (15, 28):
        scopes.append(("particle_count", f"N={n}", [i for i, row in enumerate(summaries) if row["particle_count"] == n]))
    for family in ("linear", "compact", "irregular"):
        scopes.append(("family", family, [i for i, row in enumerate(summaries) if row["family"] == family]))
    for level in range(1, 5):
        scopes.append(("target_level", str(level), [i for i, row in enumerate(summaries) if row["target_level"] == level]))
    metric_rows = []
    metric_map = {}
    audit_map: dict[str, list[object]] = {"M1": [], "P3": []}
    threshold_rows = []
    for model in ("M1", "P3"):
        for scope_type, scope, indices in scopes:
            selected = [i for i in indices if summaries[i]["eligible"]]
            observed = np.asarray([float(summaries[i]["epsilon_a_e"]) for i in selected])
            predicted = np.asarray([
                float(prediction_map[(str(summaries[i]["case_id"]), model)]["point_prediction"])
                for i in selected
            ])
            output, metric = _metric_row(model, scope_type, scope, observed, predicted)
            metric_rows.append(output)
            if metric is not None:
                metric_map[(model, scope_type, scope)] = metric
        for scope, n in (("all", None), ("N=15", 15), ("N=28", 28)):
            selected = [row for row in summaries if row["eligible"] and (n is None or row["particle_count"] == n)]
            identifiers = [str(row["case_id"]) for row in selected]
            observed = np.asarray([float(row["epsilon_a_e"]) for row in selected])
            conservative = np.asarray([
                float(prediction_map[(str(row["case_id"]), model)]["conservative_prediction"])
                for row in selected
            ])
            for tolerance in TOLERANCES:
                audit = audit_external_threshold(
                    identifiers, observed, conservative, model=model, scope=scope,
                    tolerance=tolerance,
                )
                audit_map[model].append(audit)
                output = asdict(audit)
                output["false_safe_ids"] = ";".join(audit.false_safe_ids)
                output["false_unsafe_ids"] = ";".join(audit.false_unsafe_ids)
                output["safe_margin_to_tolerance"] = tolerance - audit.worst_predicted_safe_error if audit.predicted_safe_count else 0.0
                threshold_rows.append(output)

    pair_rows = []
    for family in ("linear", "compact", "irregular"):
        for level in range(1, 5):
            first = next(row for row in summaries if row["particle_count"] == 15 and row["family"] == family and row["target_level"] == level)
            second = next(row for row in summaries if row["particle_count"] == 28 and row["family"] == family and row["target_level"] == level)
            applicable = bool(first["eligible"] and second["eligible"] and float(first["epsilon_a_e"]) > 0.0)
            ratio = float(second["epsilon_a_e"]) / float(first["epsilon_a_e"]) if applicable else 0.0
            pair_rows.append({
                "family": family, "target_level": level,
                "lambda_max": first["lambda_max"],
                "n15_case_id": first["case_id"], "n28_case_id": second["case_id"],
                "n15_epsilon_a_e": first["epsilon_a_e"], "n28_epsilon_a_e": second["epsilon_a_e"],
                "ratio_28_over_15": ratio,
                "log_ratio_28_over_15": np.log(ratio) if applicable else 0.0,
                "applicable": applicable,
                "sign": "increase" if applicable and ratio > 1 else "decrease" if applicable and ratio < 1 else "equal" if applicable else "inapplicable",
            })

    eligible_count = sum(bool(row["eligible"]) for row in summaries)
    eligible_by_n = {n: sum(bool(row["eligible"]) and row["particle_count"] == n for row in summaries) for n in (15, 28)}
    eligible_by_family = {family: sum(bool(row["eligible"]) and row["family"] == family for row in summaries) for family in ("linear", "compact", "irregular")}
    eligible_by_level = {level: sum(bool(row["eligible"]) and row["target_level"] == level for row in summaries) for level in range(1, 5)}
    predicted_safe_eligible = {}
    for tolerance, field in zip(TOLERANCES, ("safe_1pct", "safe_5pct", "safe_10pct")):
        predicted_safe_eligible[tolerance] = sum(
            row["model"] == "M1" and _truth(row[field])
            and bool(next(item for item in summaries if item["case_id"] == row["case_id"])["eligible"])
            for row in frozen
        )
    criteria, decision, next_gate = evaluate_scale_out_gate(
        eligible_count=eligible_count, eligible_by_n=eligible_by_n,
        eligible_by_family=eligible_by_family, eligible_by_level=eligible_by_level,
        predicted_safe_eligible=predicted_safe_eligible,
        manifest_intact=tuple(row["case_id"] for row in summaries) == EXPECTED_SCALE_OUT_CASE_IDS,
        phase_a_integrity=phase_a_integrity, prior_integrity=prior_integrity,
        maximum_lmax=max(int(row["final_lmax"]) for row in summaries),
        protocol_immutable=phase_a_integrity, resource_limit=False,
        m1_global=metric_map.get(("M1", "global", "all")),
        m1_by_n={n: metric_map[("M1", "particle_count", f"N={n}")] for n in (15, 28) if ("M1", "particle_count", f"N={n}") in metric_map},
        m1_audits=audit_map["M1"],
    )
    gate_rows = [{**asdict(item), "decision": "", "next_gate": ""} for item in criteria]
    gate_rows.append({
        "stage": "decision", "name": "final_decision", "observed": 1.0,
        "threshold": 1.0, "passed": True, "justification": "literal preregistered T14 gate result",
        "decision": decision, "next_gate": next_gate,
    })

    _write(FORCES, force_rows)
    _write(SUMMARY, summaries)
    _write(PREDICTIONS, revealed)
    _write(METRICS, metric_rows)
    _write(THRESHOLDS, threshold_rows)
    _write(PAIRS, pair_rows)
    _write(PERFORMANCE, performance_rows)
    _write(GATE, gate_rows)
    _plot(summaries, pair_rows, threshold_rows, performance_rows, decision)
    print(f"T14 analysis: eligible={eligible_count}/24 decision={decision}")


def _plot(summaries, pairs, thresholds, performance, decision) -> None:
    eligible = [row for row in summaries if row["eligible"]]
    colors = {15: "#1b9e77", 28: "#d95f02"}
    markers = {"linear": "o", "compact": "s", "irregular": "^"}
    figure, axes = plt.subplots(2, 3, figsize=(13.5, 8.3), constrained_layout=True)
    ax = axes[0, 0]
    for row in eligible:
        ax.scatter(row["m1_point_prediction"], row["epsilon_a_e"], color=colors[row["particle_count"]], marker=markers[row["family"]], s=34)
    line = np.logspace(-4, 0, 200)
    ax.plot(line, line, color="black", lw=1.3, label="identity")
    ax.plot(line, 2 * line, "--", color="0.45", lw=1)
    ax.plot(line, line / 2, "--", color="0.45", lw=1, label="factor 2")
    ax.set(xscale="log", yscale="log", xlabel="frozen M1 prediction", ylabel="observed $\\epsilon_A^E$", title="Frozen point prediction")
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[0, 1]
    for n in (15, 28):
        selected = [row for row in eligible if row["particle_count"] == n]
        ax.scatter([row["lambda_max"] for row in selected], [row["epsilon_a_e"] for row in selected], color=colors[n], label=f"N={n}", s=30)
    for tolerance in TOLERANCES:
        ax.axhline(tolerance, color="0.55", ls=":", lw=0.8)
    ax.set(xscale="log", yscale="log", xlabel="$\\Lambda_{\\max}$", ylabel="$\\epsilon_A^E$", title="Controlled coupling scale")
    ax.legend(fontsize=8)

    ax = axes[0, 2]
    t13 = _read(T13_SUMMARY)
    for n in (6, 10):
        selected = [row for row in t13 if int(row["particle_count"]) == n and row["eligible"] == "true"]
        ax.scatter([n] * len(selected), [float(row["m1_log_residual"]) for row in selected], color="0.65", s=18, alpha=0.7)
    for n in (15, 28):
        selected = [row for row in eligible if row["particle_count"] == n]
        ax.scatter([n] * len(selected), [row["m1_log_residual"] for row in selected], color=colors[n], s=28, label=f"N={n}")
    ax.axhline(0, color="black", lw=1)
    ax.set(xlabel="particle count N", ylabel="$\\log(\\widehat\\epsilon/\\epsilon)$", title="Frozen M1 residual versus scale")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    x = np.arange(len(pairs))
    values = [row["ratio_28_over_15"] for row in pairs]
    ax.bar(x, values, color=[{"linear": "#7570b3", "compact": "#66a61e", "irregular": "#e7298a"}[row["family"]] for row in pairs])
    ax.axhline(1, color="black", lw=1)
    ax.set(xticks=x, xticklabels=[f"{row['family'][0].upper()}{row['target_level']}" for row in pairs], ylabel="$R_{28/15}$", title="Matched size ratios")

    ax = axes[1, 1]
    for n in (15, 28):
        selected = [row for row in performance if row["particle_count"] == n]
        ax.scatter([row["final_system_dimension"] for row in selected], [row["case_wall_seconds"] for row in selected], color=colors[n], label=f"N={n}", s=30)
    ax.set(xlabel="final balanced dimension", ylabel="case solve time (s)", title="Computational scale")
    ax.legend(fontsize=8)

    ax = axes[1, 2]
    global_m1 = [row for row in thresholds if row["model"] == "M1" and row["scope"] == "all"]
    x = np.arange(3)
    ax.bar(x - 0.2, [row["predicted_safe_count"] for row in global_m1], 0.4, label="predicted safe")
    ax.bar(x + 0.2, [row["false_safe_count"] for row in global_m1], 0.4, label="false safe")
    ax.set(xticks=x, xticklabels=["1%", "5%", "10%"], ylabel="eligible cases", title=f"Conservative gate\n{decision}")
    ax.legend(fontsize=8)
    figure.suptitle("T14 scale-out audit of frozen $\\Lambda_{\\max}$ criterion", fontsize=12)
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE, dpi=220, metadata={"Software": "acoustic_ms T14"})
    plt.close(figure)


if __name__ == "__main__":
    analyze()
