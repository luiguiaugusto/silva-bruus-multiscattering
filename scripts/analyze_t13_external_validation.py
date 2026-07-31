#!/usr/bin/env python3
"""Analyze the frozen T13 holdout without calling the Model-E solver."""

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
    EXPECTED_CASE_IDS,
    TOLERANCES,
    audit_external_threshold,
    compare_model_e_forces,
    evaluate_external_validation_gate,
    external_prediction_metrics,
    mechanism_diagnostics,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
MANIFEST = DATA / "t13_holdout_manifest.csv"
FROZEN_PREDICTIONS = DATA / "t13_frozen_predictions.csv"
FROZEN_PROTOCOL = DATA / "t13_frozen_protocol.csv"
RAW = DATA / "t13_model_e_convergence.csv"
T08_FORCES = DATA / "t08_forces.csv"
T12_DEVELOPMENT = DATA / "t12_1_resolved_comparison.csv"
FORCES = DATA / "t13_forces.csv"
SUMMARY = DATA / "t13_case_summary.csv"
PREDICTIONS = DATA / "t13_external_predictions.csv"
METRICS = DATA / "t13_metrics.csv"
THRESHOLDS = DATA / "t13_threshold_audit.csv"
GATE = DATA / "t13_gate.csv"
FIGURE = FIGURES / "t13_external_validation.png"
PHASE_A_PATHS = (MANIFEST, FROZEN_PREDICTIONS, FROZEN_PROTOCOL)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


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
        raise ValueError("invalid serialized vector field")
    return result


def _truth(value: str) -> bool:
    if value not in {"true", "false"}:
        raise ValueError(f"invalid Boolean field: {value}")
    return value == "true"


def _phase_a_integrity() -> bool:
    for path in PHASE_A_PATHS:
        relative = str(path.relative_to(ROOT))
        expected = subprocess.run(
            ["git", "show", f"HEAD:{relative}"], cwd=ROOT, check=True,
            capture_output=True,
        ).stdout
        if sha256(expected).digest() != sha256(path.read_bytes()).digest():
            return False
    return subprocess.run(
        ["git", "diff", "--exit-code", "HEAD", "--", "results"],
        cwd=ROOT, capture_output=True,
    ).returncode == 0


def _frozen_force_map() -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read(T08_FORCES):
        grouped.setdefault(row["case_id"], []).append(row)
    result = {}
    for case_id, rows in grouped.items():
        rows.sort(key=lambda row: int(row["particle_index"]))
        positions = np.asarray([[float(row[key]) for key in ("x", "y", "z")] for row in rows])
        model_a = np.asarray([[float(row["a_x"]), float(row["a_y"])] for row in rows])
        model_d = np.asarray([[float(row["d_x"]), float(row["d_y"])] for row in rows])
        result[case_id] = positions, model_a, model_d
    return result


def _validated_raw() -> tuple[list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    rows = _read(RAW)
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
    if tuple(grouped) != EXPECTED_CASE_IDS or len(grouped) != 24:
        raise RuntimeError("raw T13 cases differ from the frozen manifest")
    for case_id, case_rows in grouped.items():
        case_rows.sort(key=lambda row: int(row["lmax"]))
        orders = [int(row["lmax"]) for row in case_rows]
        if orders != list(range(2, orders[-1] + 1)) or orders[-1] > 21:
            raise RuntimeError(f"invalid order sequence for {case_id}")
        if any(row["case_id"] != case_id for row in case_rows):
            raise RuntimeError("raw case grouping failed")
        for row in case_rows:
            for field, value in row.items():
                if field.endswith("forces_xyz") or field == "coordinates_xyz" or not value:
                    continue
                if field in {"case_id", "family", "stratum", "production_solver"}:
                    continue
                if value.lower() in {"true", "false"}:
                    continue
                if not np.isfinite(float(value)):
                    raise RuntimeError(f"nonfinite raw field {field} for {case_id}")
    return rows, grouped


def _scopes(summaries: list[dict[str, object]]) -> list[tuple[str, str, list[int]]]:
    scopes: list[tuple[str, str, list[int]]] = [("global", "all", list(range(len(summaries))))]
    for particle_count in (6, 10):
        scopes.append(("particle_count", f"N={particle_count}", [
            index for index, row in enumerate(summaries) if row["particle_count"] == particle_count
        ]))
    for family in ("linear", "compact", "irregular"):
        scopes.append(("family", family, [
            index for index, row in enumerate(summaries) if row["family"] == family
        ]))
    for level in (1, 2, 3, 4):
        scopes.append(("target_level", str(level), [
            index for index, row in enumerate(summaries) if row["target_level"] == level
        ]))
    return scopes


def analyze() -> None:
    if not RAW.exists():
        raise RuntimeError("T13 raw campaign is absent; analysis cannot call Model E")
    integrity = _phase_a_integrity()
    manifest = _read(MANIFEST)
    if len(manifest) != 24 or tuple(row["case_id"] for row in manifest) != EXPECTED_CASE_IDS:
        raise RuntimeError("frozen manifest identity failed")
    _, grouped = _validated_raw()
    frozen_forces = _frozen_force_map()
    frozen_predictions = _read(FROZEN_PREDICTIONS)
    prediction_map = {(row["case_id"], row["model"]): row for row in frozen_predictions}
    if len(prediction_map) != 48:
        raise RuntimeError("frozen prediction table must contain 48 unique rows")
    development_lambda = np.asarray([
        float(row["lambda_max"]) for row in _read(T12_DEVELOPMENT)
    ])
    development_min = float(np.min(development_lambda))
    development_max = float(np.max(development_lambda))

    force_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    final_arrays: dict[str, dict[str, np.ndarray]] = {}
    for manifest_row in manifest:
        case_id = manifest_row["case_id"]
        case_rows = grouped[case_id]
        final = case_rows[-1]
        positions, model_a, model_d = frozen_forces[case_id]
        if len(positions) != int(manifest_row["particle_count"]):
            raise RuntimeError(f"particle count mismatch for {case_id}")
        interaction = _vectors(final["interaction_forces_xyz"])
        external_scattered = _vectors(final["external_scattered_forces_xyz"])
        scattered_scattered = _vectors(final["scattered_scattered_forces_xyz"])
        total = _vectors(final["total_forces_xyz"])
        comparison = compare_model_e_forces(
            model_a, model_d, interaction, external_scattered, scattered_scattered
        )
        a_xyz = np.column_stack((model_a, np.zeros(len(model_a))))
        d_xyz = np.column_stack((model_d, np.zeros(len(model_d))))
        mechanism = mechanism_diagnostics(
            a_xyz, d_xyz, interaction, external_scattered, scattered_scattered
        )
        diagnostics = all(_truth(row["diagnostics_pass"]) for row in case_rows)
        confirmed = _truth(final["interaction_confirmed"])
        applicable = comparison.epsilon_a_e_applicable
        eligible = bool(confirmed and diagnostics and applicable)
        extended = int(final["final_lmax"]) > 13
        if int(final["final_lmax"]) < 13 and not all(
            _truth(final[f"{channel}_confirmed"])
            for channel in ("total", "interaction", "external_scattered", "scattered_scattered")
        ):
            raise RuntimeError(f"premature early stop for {case_id}")
        final_arrays[case_id] = {
            "a": a_xyz, "d": d_xyz, "total": total, "interaction": interaction,
            "external_scattered": external_scattered,
            "scattered_scattered": scattered_scattered,
        }
        for particle in range(len(positions)):
            row: dict[str, object] = {
                "holdout_order": int(manifest_row["holdout_order"]),
                "case_id": case_id,
                "particle_index": particle,
            }
            for prefix, values in (
                ("position", positions), ("model_a", a_xyz), ("model_d", d_xyz),
                ("model_e_total", total), ("model_e_interaction", interaction),
                ("model_e_external_scattered", external_scattered),
                ("model_e_scattered_scattered", scattered_scattered),
            ):
                for component, value in zip(("x", "y", "z"), values[particle]):
                    row[f"{prefix}_{component}"] = value
            force_rows.append(row)
        m1 = prediction_map[(case_id, "M1")]
        p3 = prediction_map[(case_id, "P3")]
        summaries.append({
            "holdout_order": int(manifest_row["holdout_order"]),
            "case_id": case_id,
            "particle_count": int(manifest_row["particle_count"]),
            "family": manifest_row["family"],
            "stratum": manifest_row["stratum"],
            "target_level": int(manifest_row["target_level"]),
            "f1": float(manifest_row["f1"]),
            "distance_ratio": float(manifest_row["distance_ratio"]),
            "lambda_max": float(manifest_row["lambda_max"]),
            "rho_l1": float(manifest_row["rho_l1"]),
            "model_d_lmax": int(manifest_row["reference_lmax_d"]),
            "model_e_final_lmax": int(final["final_lmax"]),
            "model_e_extended_beyond_13": extended,
            "interaction_confirmed": confirmed,
            "diagnostics_pass": diagnostics,
            "error_applicable": applicable,
            "eligible": eligible,
            "ineligibility_reason": "eligible" if eligible else (
                "interaction_unconfirmed" if not confirmed else
                "numerical_diagnostics_failed" if not diagnostics else "interaction_force_numerically_null"
            ),
            "epsilon_a_e": comparison.epsilon_a_e,
            "epsilon_d_e": comparison.epsilon_d_e,
            "epsilon_a_external_scattered": comparison.epsilon_a_external_scattered,
            "epsilon_d_external_scattered": comparison.epsilon_d_external_scattered,
            "m1_point_prediction": float(m1["point_prediction"]),
            "m1_conservative_prediction": float(m1["conservative_prediction"]),
            "p3_point_prediction": float(p3["point_prediction"]),
            "p3_conservative_prediction": float(p3["conservative_prediction"]),
            "m1_log_residual": np.log(float(m1["point_prediction"]) / comparison.epsilon_a_e) if eligible else 0.0,
            "lambda_development_min": development_min,
            "lambda_development_max": development_max,
            "lambda_extrapolation": not development_min <= float(manifest_row["lambda_max"]) <= development_max,
            "rms_model_e_interaction": comparison.rms_model_e_interaction,
            "x_d_minus_a": comparison.x_d_minus_a,
            "x_mie_external": comparison.x_mie_external,
            "x_scattered_scattered": comparison.x_scattered_scattered,
            "cancellation_ratio": comparison.cancellation_ratio,
            "cancellation_ratio_applicable": comparison.cancellation_ratio_applicable,
            "mechanism_mu_dm": mechanism.mu_dm.value,
            "mechanism_mu_dm_applicable": mechanism.mu_dm.applicable,
            "mechanism_mu_ds": mechanism.mu_ds.value,
            "mechanism_mu_ds_applicable": mechanism.mu_ds.applicable,
            "mechanism_mu_ms": mechanism.mu_ms.value,
            "mechanism_mu_ms_applicable": mechanism.mu_ms.applicable,
            "mechanism_projection_sum": mechanism.projection_sum.value,
            "mechanism_closure_rms": mechanism.closure_rms,
            "force_identity_relative_error": comparison.decomposition_relative_error,
            "force_identity_max_abs_error": comparison.decomposition_max_abs_error,
            "maximum_balanced_condition": max(float(row["balanced_condition_number"]) for row in case_rows),
            "maximum_backward_error": max(float(row["balanced_backward_error"]) for row in case_rows),
            "maximum_effective_incident_closure_error": max(float(row["effective_incident_closure_error"]) for row in case_rows),
            "maximum_scattering_closure_error": max(float(row["scattering_closure_error"]) for row in case_rows),
            "maximum_force_decomposition_residual": max(float(row["force_decomposition_residual"]) for row in case_rows),
            "maximum_abs_fz": max(float(row["max_abs_fz"]) for row in case_rows),
        })

    revealed_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    threshold_rows: list[dict[str, object]] = []
    scopes = _scopes(summaries)
    for frozen in frozen_predictions:
        summary = next(row for row in summaries if row["case_id"] == frozen["case_id"])
        revealed_rows.append({
            **frozen,
            "eligible": summary["eligible"],
            "observed_epsilon_a_e": summary["epsilon_a_e"],
            "observed_safe_1pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.01),
            "observed_safe_5pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.05),
            "observed_safe_10pct": bool(summary["eligible"] and float(summary["epsilon_a_e"]) < 0.10),
            "response_source": "Model_E_interaction_force_only",
        })
    metrics_by_key = {}
    audits_by_model: dict[str, list[object]] = {"M1": [], "P3": []}
    for model in ("M1", "P3"):
        for scope_type, scope, indices in scopes:
            eligible_indices = [index for index in indices if summaries[index]["eligible"]]
            if len(eligible_indices) < 2:
                raise RuntimeError(f"too few eligible cases for {model} {scope_type} {scope}")
            observed = np.asarray([float(summaries[index]["epsilon_a_e"]) for index in eligible_indices])
            predicted = np.asarray([
                float(prediction_map[(str(summaries[index]["case_id"]), model)]["point_prediction"])
                for index in eligible_indices
            ])
            metric = external_prediction_metrics(observed, predicted)
            row = {"model": model, "scope_type": scope_type, "scope": scope, **asdict(metric)}
            metric_rows.append(row)
            metrics_by_key[(model, scope_type, scope)] = metric
        for scope, particle_count in (("all", None), ("N=6", 6), ("N=10", 10)):
            selected = [row for row in summaries if row["eligible"] and (
                particle_count is None or row["particle_count"] == particle_count
            )]
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
                audits_by_model[model].append(audit)
                output = asdict(audit)
                output["false_safe_ids"] = ";".join(audit.false_safe_ids)
                output["false_unsafe_ids"] = ";".join(audit.false_unsafe_ids)
                threshold_rows.append(output)

    eligible_count = sum(bool(row["eligible"]) for row in summaries)
    eligible_n6 = sum(bool(row["eligible"]) and row["particle_count"] == 6 for row in summaries)
    eligible_n10 = sum(bool(row["eligible"]) and row["particle_count"] == 10 for row in summaries)
    predicted_safe_n = {}
    for tolerance in (0.05, 0.10):
        values = []
        for scope in ("N=6", "N=10"):
            audit = next(item for item in audits_by_model["M1"] if item.scope == scope and item.tolerance == tolerance)
            values.append(audit.predicted_safe_count)
        predicted_safe_n[tolerance] = tuple(values)
    criteria, decision, t14 = evaluate_external_validation_gate(
        eligible_count=eligible_count,
        eligible_n6=eligible_n6,
        eligible_n10=eligible_n10,
        diagnostics_all_passed=all(bool(row["diagnostics_pass"]) for row in summaries if row["eligible"]),
        manifest_intact=tuple(row["case_id"] for row in summaries) == EXPECTED_CASE_IDS,
        integrity_passed=integrity,
        m1_global=metrics_by_key.get(("M1", "global", "all")),
        m1_n6=metrics_by_key.get(("M1", "particle_count", "N=6")),
        m1_n10=metrics_by_key.get(("M1", "particle_count", "N=10")),
        m1_audits=audits_by_model["M1"],
        predicted_safe_n_by_tolerance=predicted_safe_n,
    )
    gate_rows = [{**asdict(item), "decision": "", "t14_recommendation": ""} for item in criteria]
    gate_rows.append({
        "stage": "decision", "name": "final_decision", "observed": 1.0,
        "threshold": 1.0, "passed": True, "justification": "literal preregistered gate result",
        "decision": decision, "t14_recommendation": t14,
    })

    _write(FORCES, force_rows)
    _write(SUMMARY, summaries)
    _write(PREDICTIONS, revealed_rows)
    _write(METRICS, metric_rows)
    _write(THRESHOLDS, threshold_rows)
    _write(GATE, gate_rows)
    _plot(summaries, threshold_rows, decision)
    print(f"T13 analysis: eligible={eligible_count}/24 decision={decision}")


def _plot(
    summaries: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
    decision: str,
) -> None:
    eligible = [row for row in summaries if row["eligible"]]
    colors = {6: "#2166ac", 10: "#b2182b"}
    markers = {"linear": "o", "compact": "s", "irregular": "^"}
    figure, axes = plt.subplots(2, 3, figsize=(13.2, 8.2), constrained_layout=True)
    ax = axes[0, 0]
    for row in eligible:
        ax.scatter(row["m1_point_prediction"], row["epsilon_a_e"], color=colors[int(row["particle_count"])], marker=markers[str(row["family"])], s=31)
    limits = np.logspace(-4, 0, 200)
    ax.plot(limits, limits, color="black", lw=1.4, label="identity")
    ax.plot(limits, 2 * limits, "--", color="0.45", lw=1)
    ax.plot(limits, limits / 2, "--", color="0.45", lw=1, label="factor 2")
    ax.set(xscale="log", yscale="log", xlabel="frozen M1 point prediction", ylabel="observed $\\epsilon_A^E$", title="Point prediction")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for row in eligible:
        ax.scatter(row["m1_conservative_prediction"], row["epsilon_a_e"], color=colors[int(row["particle_count"])], marker=markers[str(row["family"])], s=31)
    for tolerance in TOLERANCES:
        ax.axvline(tolerance, color="0.55", lw=0.8, ls=":")
        ax.axhline(tolerance, color="0.55", lw=0.8, ls=":")
    ax.set(xscale="log", yscale="log", xlabel="frozen conservative M1", ylabel="observed $\\epsilon_A^E$", title="Strict safety regions")

    ax = axes[0, 2]
    for particle_count in (6, 10):
        for family in ("linear", "compact", "irregular"):
            selected = [row for row in eligible if row["particle_count"] == particle_count and row["family"] == family]
            ax.scatter([row["target_level"] for row in selected], [row["m1_log_residual"] for row in selected], color=colors[particle_count], marker=markers[family], label=f"N={particle_count} {family}")
    ax.axhline(0.0, color="black", lw=1)
    ax.set(xlabel="target level", ylabel="$\\log(\\widehat\\epsilon/\\epsilon)$", title="M1 residuals")
    ax.legend(fontsize=6, ncol=2)

    ax = axes[1, 0]
    global_m1 = [row for row in threshold_rows if row["model"] == "M1" and row["scope"] == "all"]
    x = np.arange(3)
    ax.bar(x - 0.18, [int(row["predicted_safe_count"]) for row in global_m1], 0.36, label="predicted safe")
    ax.bar(x + 0.18, [int(row["false_safe_count"]) for row in global_m1], 0.36, label="false safe")
    ax.set(xticks=x, xticklabels=["1%", "5%", "10%"], ylabel="eligible cases", title="M1 safety audit")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    for model, xfield, color in (("M1", "m1_point_prediction", "#4d4d4d"), ("P3", "p3_point_prediction", "#7b3294")):
        factors = [np.exp(abs(np.log(float(row[xfield]) / float(row["epsilon_a_e"])))) for row in eligible]
        ax.scatter([float(row[xfield]) for row in eligible], factors, s=24, alpha=0.75, color=color, label=model)
    ax.axhline(2.0, color="black", ls="--", lw=1)
    ax.set(xscale="log", yscale="log", xlabel="frozen point prediction", ylabel="multiplicative factor", title="Transparent M1/P3 comparison")
    ax.legend(fontsize=8)

    ax = axes[1, 2]
    for row in summaries:
        color = colors[int(row["particle_count"])] if row["eligible"] else "0.6"
        marker = markers[str(row["family"])] if row["eligible"] else "x"
        ax.scatter(row["target_level"], row["model_e_final_lmax"], color=color, marker=marker, s=38)
    ax.axhline(13, color="black", ls=":", lw=1, label="standard cap")
    ax.set(xlabel="target level", ylabel="final $L_{\\max}$", title="Convergence eligibility")
    ax.legend(fontsize=8)
    figure.suptitle(f"T13 external validation of frozen $\\Lambda_{{\\max}}$ M1\n{decision}", fontsize=12)
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE, dpi=220, metadata={"Software": "acoustic_ms T13"})
    plt.close(figure)


def main() -> None:
    analyze()


if __name__ == "__main__":
    main()
