#!/usr/bin/env python3
"""Analyze the frozen T14.1 campaign without importing or calling Model E."""

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
    EXPECTED_LARGE_N_CASE_IDS,
    TOLERANCES,
    audit_external_threshold,
    classify_large_n_trend,
    evaluate_large_n_gate,
    external_eligibility_mask,
    external_prediction_metrics,
    normalized_rms_error_xyz,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
MANIFEST = DATA / "t14_1_large_n_manifest.csv"
LOCAL = DATA / "t14_1_local_coupling.csv"
FROZEN = DATA / "t14_1_frozen_predictions.csv"
PROTOCOL = DATA / "t14_1_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_1_prior_artifact_hashes.csv"
RAW = DATA / "t14_1_model_e_convergence.csv"
FORCES = DATA / "t14_1_forces.csv"
SUMMARY = DATA / "t14_1_case_summary.csv"
PREDICTIONS = DATA / "t14_1_large_n_predictions.csv"
METRICS = DATA / "t14_1_metrics.csv"
THRESHOLDS = DATA / "t14_1_threshold_audit.csv"
PAIRS = DATA / "t14_1_matched_large_n_pairs.csv"
COMBINED = DATA / "t14_1_combined_scale_sequence.csv"
PERFORMANCE = DATA / "t14_1_performance.csv"
GATE = DATA / "t14_1_gate.csv"
FIGURE = FIGURES / "t14_1_large_n_validation.png"
T14_SUMMARY = DATA / "t14_case_summary.csv"
PHASE_A = (MANIFEST, LOCAL, FROZEN, PROTOCOL, PRIOR_HASHES)


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
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format(row.get(field, "")) for field in fields})
    temporary.replace(path)


def _vectors(value: str) -> np.ndarray:
    result = np.asarray([[float(item) for item in row.split(":")] for row in value.split(";")])
    if result.ndim != 2 or result.shape[1] != 3 or not np.all(np.isfinite(result)):
        raise ValueError("invalid serialized vectors")
    return result


def _phase_a_integrity() -> bool:
    for path in PHASE_A:
        expected = subprocess.run(["git", "show", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT, check=True, capture_output=True).stdout
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
    if tuple(grouped) != EXPECTED_LARGE_N_CASE_IDS:
        raise RuntimeError("T14.1 raw identity differs from preregistration")
    for case_id, rows in grouped.items():
        rows.sort(key=lambda row: int(row["lmax"]))
        orders = [int(row["lmax"]) for row in rows]
        if orders != list(range(2, orders[-1] + 1)) or orders[-1] > 13:
            raise RuntimeError(f"invalid Lmax sequence for {case_id}")
        for row in rows:
            for field, value in row.items():
                if field.endswith("forces_xyz") or field == "coordinates_xyz":
                    continue
                if field in {"case_id", "family", "stop_reason", "production_solver", "coordinate_sha256", "local_coupling_sha256"}:
                    continue
                if value.lower() in {"true", "false"}:
                    continue
                if not np.isfinite(float(value)):
                    raise RuntimeError(f"nonfinite {field} for {case_id}")
    return grouped


def _metric_row(model: str, scope_type: str, scope: str, observed: np.ndarray, predicted: np.ndarray) -> tuple[dict[str, object], object | None]:
    if len(observed) < 2:
        return {
            "model": model, "scope_type": scope_type, "scope": scope,
            "applicable": False, "reason": "fewer_than_two_eligible_cases",
            "point_count": len(observed), "rmse_log": 0.0, "mae_log": 0.0,
            "median_factor": 0.0, "p90_factor": 0.0, "maximum_factor": 0.0,
            "fraction_within_factor_2": 0.0, "fraction_within_factor_1_5": 0.0,
            "spearman": 0.0, "mean_log_bias": 0.0, "median_log_bias": 0.0,
        }, None
    metric = external_prediction_metrics(observed, predicted)
    return {"model": model, "scope_type": scope_type, "scope": scope, "applicable": True, "reason": "applicable", **asdict(metric)}, metric


def analyze() -> None:
    if not RAW.exists():
        raise RuntimeError("T14.1 raw campaign is absent; analysis never calls Model E")
    phase_a_integrity = _phase_a_integrity()
    prior_integrity = _prior_integrity()
    manifest = _read(MANIFEST)
    if len(manifest) != 24 or tuple(row["case_id"] for row in manifest) != EXPECTED_LARGE_N_CASE_IDS:
        raise RuntimeError("T14.1 manifest identity failed")
    grouped = _validated_raw()
    frozen = _read(FROZEN)
    prediction_map = {(row["case_id"], row["model"]): row for row in frozen}
    if len(prediction_map) != 48:
        raise RuntimeError("T14.1 must contain 48 unique frozen predictions")

    summaries: list[dict[str, object]] = []
    force_rows: list[dict[str, object]] = []
    performance_rows: list[dict[str, object]] = []
    resource_limit = False
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
        campaign_complete = _truth(final["campaign_complete"])
        integrity = (
            final["coordinate_sha256"] == source["coordinate_sha256"]
            and final["local_coupling_sha256"] == source["local_coupling_sha256"]
            and np.isclose(float(final["lambda_max"]), float(source["lambda_max"]), rtol=5e-13, atol=5e-15)
            and np.isclose(float(final["rho_l1"]), float(source["rho_l1"]), rtol=3e-12, atol=3e-13)
        )
        eligible = bool(external_eligibility_mask([confirmed and campaign_complete and integrity], [diagnostics], [applicable])[0])
        resource_failed = final["stop_reason"] == "resource_precheck_failed"
        resource_limit |= resource_failed
        reason = "eligible"
        if resource_failed:
            reason = "resource_precheck_failed"
        elif not campaign_complete:
            reason = "campaign_incomplete"
        elif not integrity:
            reason = "phase_a_identity_mismatch"
        elif not confirmed:
            reason = "interaction_unconfirmed_at_13"
        elif not diagnostics:
            reason = "final_numerical_diagnostics_failed"
        elif not applicable:
            reason = "interaction_force_numerically_null"
        m1 = prediction_map[(case_id, "M1")]
        p3 = prediction_map[(case_id, "P3")]
        summary = {
            "scale_order": int(source["scale_order"]), "case_id": case_id,
            "particle_count": int(source["particle_count"]), "family": source["family"],
            "target_level": int(source["target_level"]), "f1": float(source["f1"]),
            "distance_ratio": float(source["distance_ratio"]),
            "minimum_distance": float(source["minimum_distance"]),
            "cluster_diameter": float(source["cluster_diameter"]),
            "k_cluster_diameter": float(source["k_cluster_diameter"]),
            "lambda_min": float(source["lambda_min"]), "lambda_mean": float(source["lambda_mean"]),
            "lambda_median": float(source["lambda_median"]), "lambda_std": float(source["lambda_std"]),
            "lambda_p10": float(source["lambda_p10"]), "lambda_p90": float(source["lambda_p90"]),
            "lambda_max": float(source["lambda_max"]),
            "lambda_mean_over_max": float(source["lambda_mean_over_max"]),
            "fraction_ge_0_9_max": float(source["fraction_ge_0_9_max"]),
            "first_argmax": int(source["first_argmax"]), "rho_l1": float(source["rho_l1"]),
            "final_lmax": int(final["final_lmax"]), "stop_reason": final["stop_reason"],
            "interaction_confirmed": confirmed, "total_confirmed": _truth(final["total_confirmed"]),
            "external_scattered_confirmed": _truth(final["external_scattered_confirmed"]),
            "scattered_scattered_confirmed": _truth(final["scattered_scattered_confirmed"]),
            "diagnostics_pass": diagnostics, "phase_a_identity_pass": integrity,
            "error_applicable": applicable, "eligible": eligible,
            "ineligibility_reason": reason, "epsilon_a_e": epsilon,
            "m1_point_prediction": float(m1["point_prediction"]),
            "m1_conservative_prediction": float(m1["conservative_prediction"]),
            "p3_point_prediction": float(p3["point_prediction"]),
            "p3_conservative_prediction": float(p3["conservative_prediction"]),
            "m1_log_residual": np.log(float(m1["point_prediction"]) / epsilon) if eligible else 0.0,
            "final_system_dimension": int(final["system_dimension"]),
            "case_wall_seconds": float(final["case_accumulated_seconds"]),
            "peak_memory_kib": max(int(row["process_peak_memory_kib"]) for row in rows),
            "maximum_estimated_memory_bytes": max(int(row["estimated_memory_bytes"]) for row in rows),
            "maximum_balanced_condition": max(float(row["balanced_condition_number"]) for row in rows),
            "maximum_backward_error": max(float(row["balanced_backward_error"]) for row in rows),
            "maximum_incident_closure_error": max(float(row["effective_incident_closure_error"]) for row in rows),
            "maximum_scattering_closure_error": max(float(row["scattering_closure_error"]) for row in rows),
            "maximum_force_decomposition_residual": max(float(row["force_decomposition_residual"]) for row in rows),
            "maximum_abs_fz": max(float(row["max_abs_fz"]) for row in rows),
        }
        summaries.append(summary)
        for particle in range(len(positions)):
            output: dict[str, object] = {"scale_order": int(source["scale_order"]), "case_id": case_id, "particle_index": particle}
            for prefix, values in (
                ("position", positions), ("model_a", model_a), ("model_e_total", total),
                ("model_e_interaction", interaction), ("model_e_external_scattered", external_scattered),
                ("model_e_scattered_scattered", scattered_scattered),
            ):
                for axis, value in zip(("x", "y", "z"), values[particle]):
                    output[f"{prefix}_{axis}"] = value
            force_rows.append(output)
        performance_rows.append({
            "scale_order": int(source["scale_order"]), "case_id": case_id,
            "particle_count": int(source["particle_count"]), "family": source["family"],
            "target_level": int(source["target_level"]), "orders_computed": len(rows),
            "final_lmax": int(final["final_lmax"]), "final_system_dimension": int(final["system_dimension"]),
            "case_wall_seconds": float(final["case_accumulated_seconds"]),
            "mean_order_wall_seconds": float(np.mean([float(row["order_wall_seconds"]) for row in rows])),
            "maximum_order_wall_seconds": max(float(row["order_wall_seconds"]) for row in rows),
            "peak_memory_kib": max(int(row["process_peak_memory_kib"]) for row in rows),
            "maximum_estimated_memory_bytes": max(int(row["estimated_memory_bytes"]) for row in rows),
            "blas_threads": int(final["blas_threads"]), "workers": int(final["workers"]),
            "resource_precheck_failed": resource_failed, "stop_reason": final["stop_reason"],
        })

    revealed = []
    for row in frozen:
        summary = next(item for item in summaries if item["case_id"] == row["case_id"])
        revealed.append({
            **row, "eligible": summary["eligible"], "observed_epsilon_a_e": summary["epsilon_a_e"],
            "observed_safe_1pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.01),
            "observed_safe_5pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.05),
            "observed_safe_10pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.10),
            "response_source": "Model_E_interaction_force_only",
        })

    scopes: list[tuple[str, str, list[int]]] = [("global", "all", list(range(24)))]
    for n in (45, 105):
        scopes.append(("particle_count", f"N={n}", [i for i, row in enumerate(summaries) if row["particle_count"] == n]))
    for family in ("linear", "compact", "irregular"):
        scopes.append(("family", family, [i for i, row in enumerate(summaries) if row["family"] == family]))
    for level in range(1, 5):
        scopes.append(("target_level", str(level), [i for i, row in enumerate(summaries) if row["target_level"] == level]))
    metric_rows: list[dict[str, object]] = []
    metric_map: dict[tuple[str, str, str], object] = {}
    audit_map: dict[str, list[object]] = {"M1": [], "P3": []}
    threshold_rows: list[dict[str, object]] = []
    for model in ("M1", "P3"):
        for scope_type, scope, indices in scopes:
            selected = [i for i in indices if summaries[i]["eligible"]]
            observed = np.asarray([float(summaries[i]["epsilon_a_e"]) for i in selected])
            predicted = np.asarray([float(prediction_map[(str(summaries[i]["case_id"]), model)]["point_prediction"]) for i in selected])
            output, metric = _metric_row(model, scope_type, scope, observed, predicted)
            metric_rows.append(output)
            if metric is not None:
                metric_map[(model, scope_type, scope)] = metric
        for scope, n in (("all", None), ("N=45", 45), ("N=105", 105)):
            selected = [row for row in summaries if row["eligible"] and (n is None or row["particle_count"] == n)]
            identifiers = [str(row["case_id"]) for row in selected]
            observed = np.asarray([float(row["epsilon_a_e"]) for row in selected])
            conservative = np.asarray([float(prediction_map[(str(row["case_id"]), model)]["conservative_prediction"]) for row in selected])
            for tolerance in TOLERANCES:
                audit = audit_external_threshold(identifiers, observed, conservative, model=model, scope=scope, tolerance=tolerance)
                audit_map[model].append(audit)
                output = asdict(audit)
                output["false_safe_ids"] = ";".join(audit.false_safe_ids)
                output["false_unsafe_ids"] = ";".join(audit.false_unsafe_ids)
                output["safe_margin_to_tolerance"] = tolerance - audit.worst_predicted_safe_error if audit.predicted_safe_count else 0.0
                threshold_rows.append(output)

    t14 = _read(T14_SUMMARY)
    pair_rows: list[dict[str, object]] = []
    combined_rows: list[dict[str, object]] = []
    ratios_105_45 = []
    for family in ("linear", "compact", "irregular"):
        for level in range(1, 5):
            n45 = next(row for row in summaries if row["particle_count"] == 45 and row["family"] == family and row["target_level"] == level)
            n105 = next(row for row in summaries if row["particle_count"] == 105 and row["family"] == family and row["target_level"] == level)
            n28 = next(row for row in t14 if int(row["particle_count"]) == 28 and row["family"] == family and int(row["target_level"]) == level)
            applicable = bool(n45["eligible"] and n105["eligible"] and float(n45["epsilon_a_e"]) > 0.0)
            ratio = float(n105["epsilon_a_e"]) / float(n45["epsilon_a_e"]) if applicable else 0.0
            if applicable:
                ratios_105_45.append(ratio)
            pair_rows.append({
                "family": family, "target_level": level, "lambda_max": n45["lambda_max"],
                "n28_case_id": n28["case_id"], "n45_case_id": n45["case_id"], "n105_case_id": n105["case_id"],
                "n28_epsilon_a_e": float(n28["epsilon_a_e"]), "n45_epsilon_a_e": n45["epsilon_a_e"],
                "n105_epsilon_a_e": n105["epsilon_a_e"],
                "ratio_45_over_28": float(n45["epsilon_a_e"]) / float(n28["epsilon_a_e"]) if n45["eligible"] and n28["eligible"] == "true" else 0.0,
                "ratio_105_over_45": ratio,
                "ratio_105_over_28": float(n105["epsilon_a_e"]) / float(n28["epsilon_a_e"]) if n105["eligible"] and n28["eligible"] == "true" else 0.0,
                "log_ratio_105_over_45": np.log(ratio) if applicable else 0.0,
                "applicable": applicable,
                "sign": "increase" if applicable and ratio > 1 else "decrease" if applicable and ratio < 1 else "equal" if applicable else "inapplicable",
            })
            sequence = []
            for n in (15, 28):
                old = next(row for row in t14 if int(row["particle_count"]) == n and row["family"] == family and int(row["target_level"]) == level)
                sequence.append((n, old["case_id"], float(old["epsilon_a_e"]), old["eligible"] == "true", "T14"))
            sequence.extend(((45, n45["case_id"], float(n45["epsilon_a_e"]), bool(n45["eligible"]), "T14.1"), (105, n105["case_id"], float(n105["epsilon_a_e"]), bool(n105["eligible"]), "T14.1")))
            eligible_sequence = [(n, error) for n, _, error, flag, _ in sequence if flag and error > 0.0]
            if len(eligible_sequence) >= 2:
                exponent, intercept = np.polyfit(np.log([item[0] for item in eligible_sequence]), np.log([item[1] for item in eligible_sequence]), 1)
            else:
                exponent, intercept = 0.0, 0.0
            for n, case_id, error, flag, source_name in sequence:
                combined_rows.append({
                    "family": family, "target_level": level, "particle_count": n,
                    "case_id": case_id, "epsilon_a_e": error, "eligible": flag,
                    "source_task": source_name, "descriptive_log_n_exponent": exponent,
                    "descriptive_log_n_prefactor": np.exp(intercept) if eligible_sequence else 0.0,
                })

    eligible_count = sum(bool(row["eligible"]) for row in summaries)
    eligible_by_n = {n: sum(bool(row["eligible"]) and row["particle_count"] == n for row in summaries) for n in (45, 105)}
    eligible_by_family = {family: sum(bool(row["eligible"]) and row["family"] == family for row in summaries) for family in ("linear", "compact", "irregular")}
    eligible_by_level = {level: sum(bool(row["eligible"]) and row["target_level"] == level for row in summaries) for level in range(1, 5)}
    predicted_safe_eligible = {}
    for tolerance, field in zip(TOLERANCES, ("safe_1pct", "safe_5pct", "safe_10pct")):
        predicted_safe_eligible[tolerance] = sum(
            row["model"] == "M1" and _truth(row[field])
            and bool(next(item for item in summaries if item["case_id"] == row["case_id"])["eligible"])
            for row in frozen
        )
    criteria, decision, next_gate = evaluate_large_n_gate(
        eligible_count=eligible_count, eligible_by_n=eligible_by_n,
        eligible_by_family=eligible_by_family, eligible_by_level=eligible_by_level,
        predicted_safe_eligible=predicted_safe_eligible,
        manifest_intact=tuple(row["case_id"] for row in summaries) == EXPECTED_LARGE_N_CASE_IDS,
        phase_a_integrity=phase_a_integrity, prior_integrity=prior_integrity,
        maximum_lmax=max(int(row["final_lmax"]) for row in summaries),
        protocol_immutable=phase_a_integrity, resource_limit=resource_limit,
        m1_global=metric_map.get(("M1", "global", "all")),
        m1_by_n={n: metric_map[("M1", "particle_count", f"N={n}")] for n in (45, 105) if ("M1", "particle_count", f"N={n}") in metric_map},
        m1_audits=audit_map["M1"],
    )
    trend = classify_large_n_trend(ratios_105_45)
    gate_rows = [{**asdict(item), "decision": "", "next_gate": "", "large_n_trend": ""} for item in criteria]
    gate_rows.append({
        "stage": "diagnostic", "name": "matched_large_n_trend", "observed": len(ratios_105_45),
        "threshold": 10, "passed": trend != "INCONCLUSIVE_LARGE_N_TREND",
        "justification": "frozen diagnostic N=105/N=45 trend classification",
        "decision": "", "next_gate": "", "large_n_trend": trend,
    })
    gate_rows.append({
        "stage": "decision", "name": "final_decision", "observed": 1.0,
        "threshold": 1.0, "passed": True, "justification": "literal preregistered T14.1 gate result",
        "decision": decision, "next_gate": next_gate, "large_n_trend": trend,
    })

    _write(FORCES, force_rows)
    _write(SUMMARY, summaries)
    _write(PREDICTIONS, revealed)
    _write(METRICS, metric_rows)
    _write(THRESHOLDS, threshold_rows)
    _write(PAIRS, pair_rows)
    _write(COMBINED, combined_rows)
    _write(PERFORMANCE, performance_rows)
    _write(GATE, gate_rows)
    _plot(summaries, threshold_rows, pair_rows, combined_rows, performance_rows, decision, trend)
    print(f"T14.1 analysis: eligible={eligible_count}/24 decision={decision} trend={trend}")


def _plot(summaries, thresholds, pairs, combined, performance, decision, trend) -> None:
    eligible = [row for row in summaries if row["eligible"]]
    colors = {45: "#1b9e77", 105: "#d95f02"}
    markers = {"linear": "o", "compact": "s", "irregular": "^"}
    figure, axes = plt.subplots(2, 3, figsize=(13.7, 8.4), constrained_layout=True)
    ax = axes[0, 0]
    for row in eligible:
        ax.scatter(row["m1_point_prediction"], row["epsilon_a_e"], color=colors[row["particle_count"]], marker=markers[row["family"]], s=34)
    line = np.logspace(-4, 0, 200)
    ax.plot(line, line, color="black", lw=1.3, label="identity")
    ax.plot(line, 2 * line, "--", color="0.45", lw=1)
    ax.plot(line, line / 2, "--", color="0.45", lw=1, label="factor 2")
    ax.set(xscale="log", yscale="log", xlabel="frozen M1 prediction", ylabel="observed $\\epsilon_A^E$", title="Frozen M1 transfer")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for n in (45, 105):
        selected = [row for row in eligible if row["particle_count"] == n]
        ax.scatter([row["lambda_max"] for row in selected], [row["epsilon_a_e"] for row in selected], color=colors[n], label=f"N={n}", s=30)
    for tolerance in TOLERANCES:
        ax.axhline(tolerance, color="0.55", ls=":", lw=0.8)
    ax.set(xscale="log", yscale="log", xlabel="$\\Lambda_{\\max}$", ylabel="$\\epsilon_A^E$", title="Large-N controlled coupling")
    ax.legend(fontsize=8)
    ax = axes[0, 2]
    for row in eligible:
        ax.scatter(row["lambda_mean_over_max"], row["epsilon_a_e"], color=colors[row["particle_count"]], marker=markers[row["family"]], s=30)
    ax.set(yscale="log", xlabel="$\\overline\\Lambda/\\Lambda_{\\max}$", ylabel="$\\epsilon_A^E$", title="Local-coupling structure")
    ax = axes[1, 0]
    x = np.arange(len(pairs))
    ax.bar(x, [row["ratio_105_over_45"] for row in pairs], color=[{"linear": "#7570b3", "compact": "#66a61e", "irregular": "#e7298a"}[row["family"]] for row in pairs])
    ax.axhline(1.0, color="black", lw=1)
    ax.set(xticks=x, xticklabels=[f"{row['family'][0].upper()}{row['target_level']}" for row in pairs], ylabel="$R_{105/45}$", title=f"Matched trend: {trend.replace('_', ' ').lower()}")
    ax = axes[1, 1]
    for family, color in (("linear", "#7570b3"), ("compact", "#66a61e"), ("irregular", "#e7298a")):
        selected = [row for row in combined if row["family"] == family and row["target_level"] == 2 and row["eligible"]]
        ax.plot([row["particle_count"] for row in selected], [row["epsilon_a_e"] for row in selected], "o-", color=color, label=family)
    ax.set(xscale="log", yscale="log", xlabel="particle count N", ylabel="$\\epsilon_A^E$", title="N=15 to 105, level 2")
    ax.legend(fontsize=8)
    ax = axes[1, 2]
    for n in (45, 105):
        selected = [row for row in performance if row["particle_count"] == n]
        ax.scatter([row["final_system_dimension"] for row in selected], [row["case_wall_seconds"] for row in selected], color=colors[n], label=f"N={n}", s=30)
    label = decision.split("_", 1)[0]
    ax.set(xlabel="final balanced dimension", ylabel="case wall time (s)", title=f"Resource profile — {label}")
    ax.legend(fontsize=8)
    figure.suptitle("T14.1 large-N audit of frozen $\\Lambda_{\\max}$ criterion", fontsize=12)
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE, dpi=220, metadata={"Software": "acoustic_ms T14.1"})
    plt.close(figure)


if __name__ == "__main__":
    analyze()
