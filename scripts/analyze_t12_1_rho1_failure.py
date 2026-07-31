#!/usr/bin/env python3
"""Extend ten T12 sentinels and diagnose the frozen rho1 prediction failure."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    convergence_tail_diagnostics,
    fit_log_linear,
    leave_group_out_folds,
    mechanism_diagnostics,
    out_of_fold_metrics,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    spearman_correlation,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_t12_model_e_sentinels as t12  # noqa: E402
from analyze_t11_model_e import _successive  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
T12_RAW = DATA / "t12_model_e_convergence.csv"
T12_COMPARISON = DATA / "t12_model_comparison.csv"
EXTENDED_PATH = DATA / "t12_1_extended_convergence.csv"
SUMMARY_PATH = DATA / "t12_1_convergence_summary.csv"
RESOLVED_PATH = DATA / "t12_1_resolved_comparison.csv"
MECHANISM_PATH = DATA / "t12_1_mechanism_diagnostics.csv"
PREDICTOR_PATH = DATA / "t12_1_predictor_diagnostics.csv"
OOF_PATH = DATA / "t12_1_out_of_fold_predictions.csv"
FIGURE_PATH = FIGURES / "t12_1_rho1_failure_diagnostics.png"

TOLERANCE = 1.0e-5
MAX_LMAX = 21
INTERACTION_UNCONFIRMED = (
    "n2_pair_f1.0_d2.1",
    "n3_compact_f0.8_d2.1",
    "n3_irregular_f1.0_d2.1",
    "n3_linear_f1.0_d2.1",
    "n4_irregular_f0.8_d2.1",
    "n4_linear_f0.8_d2.1",
)
SS_UNCONFIRMED = (
    "n3_compact_f0.1_d2.1",
    "n4_compact_f0.1_d2.1",
    "n4_irregular_f0.1_d2.1",
    "n4_linear_f0.1_d2.1",
)
EXTENDED_CASES = INTERACTION_UNCONFIRMED + SS_UNCONFIRMED
SPECIAL_CASES = {
    "n2_pair_f1.0_d6.0": "largest_frozen_prediction_factor",
    "n2_pair_f0.8_d2.5": "ten_percent_false_safe",
    "n2_pair_f1.0_d2.1": "extended_large_error_check",
}
PREDICTORS = {
    "P1_eta": "eta",
    "P2_lambda_max": "lambda_max",
    "P3_rho_l1": "rho_l1",
    "P4_epsilon_a_d": "epsilon_a_t08",
}


def _format(value: float) -> str:
    return format(float(value), ".17g")


def _bool(value: bool) -> str:
    return str(bool(value)).lower()


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _group(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
    return grouped


def _confirmation(rows: list[dict[str, object]], short: str) -> int:
    for index in range(2, len(rows)):
        if (
            rows[index - 1][f"{short}_change_applicable"] == "true"
            and rows[index][f"{short}_change_applicable"] == "true"
            and float(rows[index - 1][f"{short}_successive_change"]) <= TOLERANCE
            and float(rows[index][f"{short}_successive_change"]) <= TOLERANCE
        ):
            return int(rows[index]["lmax"])
    return 0


def _all_confirmed(rows: list[dict[str, object]]) -> bool:
    return all(_confirmation(rows, short) > 0 for short, _ in t12.CHANNELS)


def _new_row(sentinel: object, positions: np.ndarray, result: object, previous: dict[str, object]) -> dict[str, object]:
    solution = result.solution
    if solution.production_solver != "balanced_sqrt":
        raise RuntimeError("Model E did not use balanced_sqrt")
    diagnostics = (
        solution.balanced_condition_number < 10.0
        and solution.balanced_backward_error < 1.0e-12
        and solution.effective_incident_closure_error < 1.0e-12
        and solution.scattering_closure_error < 1.0e-12
        and result.decomposition_residual < 1.0e-12
    )
    if not diagnostics:
        raise RuntimeError(f"numerical diagnostics failed for {sentinel.case_id} L={result.lmax}")
    row: dict[str, object] = {
        "case_id": sentinel.case_id,
        "particle_count": sentinel.particle_count,
        "family": sentinel.family,
        "rho_band": sentinel.rho_band,
        "f1": _format(sentinel.f1),
        "distance_ratio": _format(sentinel.distance_ratio),
        "rho_l1": _format(sentinel.rho_l1),
        "coordinates_xyz": t12._serialize(positions),
        "lmax": result.lmax,
        "final_lmax": result.lmax,
        "system_dimension": solution.balanced_system_matrix.shape[0],
        "active_modes_per_particle": len(solution.active_modes),
        "balanced_condition_number": _format(solution.balanced_condition_number),
        "balanced_backward_error": _format(solution.balanced_backward_error),
        "effective_incident_closure_error": _format(solution.effective_incident_closure_error),
        "scattering_closure_error": _format(solution.scattering_closure_error),
        "force_decomposition_residual": _format(result.decomposition_residual),
        "max_abs_fz": _format(max(
            float(np.max(np.abs(getattr(result, attribute)[:, 2])))
            for _, attribute in t12.CHANNELS
        )),
        "total_forces_xyz": t12._serialize(result.total_forces_xyz),
        "external_forces_xyz": t12._serialize(result.external_forces_xyz),
        "interaction_forces_xyz": t12._serialize(result.interaction_forces_xyz),
        "external_scattered_forces_xyz": t12._serialize(result.external_scattered_forces_xyz),
        "scattered_scattered_forces_xyz": t12._serialize(result.scattered_scattered_forces_xyz),
        "production_solver": solution.production_solver,
        "finite": "true",
        "campaign_complete": "true",
        "source": "t12_1",
    }
    for short, attribute in t12.CHANNELS:
        current = getattr(result, attribute)
        old = t12._deserialize(str(previous[f"{short}_forces_xyz"]))
        change, applicable, absolute = _successive(current, old)
        row[f"{short}_rms"] = _format(rms_vector_magnitude_xyz(current))
        row[f"{short}_successive_change"] = _format(change)
        row[f"{short}_absolute_change"] = _format(absolute)
        row[f"{short}_change_applicable"] = _bool(applicable)
        row[f"{short}_minimum_confirmed_lmax"] = 0
        row[f"{short}_confirmed"] = "false"
    return row


def _finalize_case(rows: list[dict[str, object]]) -> None:
    final_lmax = int(rows[-1]["lmax"])
    confirmations = {
        short: _confirmation(rows, short) for short, _ in t12.CHANNELS
    }
    for row in rows:
        row["final_lmax"] = final_lmax
        for short in confirmations:
            row[f"{short}_minimum_confirmed_lmax"] = confirmations[short]
            row[f"{short}_confirmed"] = _bool(confirmations[short] > 0)


def _validate_twelve_problem_cases(comparison: list[dict[str, str]]) -> None:
    by_id = {row["case_id"]: row for row in comparison}
    interaction = tuple(
        case_id for case_id in (item.case_id for item in t12.SENTINELS)
        if by_id[case_id]["interaction_confirmed"] != "true"
    )
    ss_only = tuple(
        case_id for case_id in (item.case_id for item in t12.SENTINELS)
        if by_id[case_id]["interaction_confirmed"] == "true"
        and by_id[case_id]["scattered_scattered_confirmed"] != "true"
    )
    if interaction != INTERACTION_UNCONFIRMED or ss_only != SS_UNCONFIRMED:
        raise RuntimeError("the ten preregistered extension cases disagree with T12")


def run_extension() -> None:
    original = _group(_read(T12_RAW))
    comparison = _read(T12_COMPARISON)
    _validate_twelve_problem_cases(comparison)
    sentinel_by_id = {item.case_id: item for item in t12.SENTINELS}
    positions_by_id = {
        case_id: t12._deserialize(rows[0]["coordinates_xyz"])
        for case_id, rows in original.items()
    }
    all_rows: list[dict[str, object]] = []
    for number, case_id in enumerate(EXTENDED_CASES, start=1):
        print(f"T12.1 extension {number:02d}/10: {case_id}", flush=True)
        inherited: list[dict[str, object]] = []
        for source_row in original[case_id]:
            copied: dict[str, object] = dict(source_row)
            copied["source"] = "t12"
            inherited.append(copied)
        rows = inherited
        sentinel = sentinel_by_id[case_id]
        positions = positions_by_id[case_id]
        for lmax in range(14, MAX_LMAX + 1):
            result = solve_model_e_nodal(
                positions,
                t12.KA,
                t12.RADIUS,
                t12.ENERGY_DENSITY,
                t12.F0,
                sentinel.f1,
                lmax,
            )
            rows.append(_new_row(sentinel, positions, result, rows[-1]))
            if _all_confirmed(rows):
                break
        _finalize_case(rows)
        all_rows.extend(rows)
    _atomic_write(EXTENDED_PATH, all_rows)


def _tail_summary(extended: dict[str, list[dict[str, str]]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_id in EXTENDED_CASES:
        case_rows = extended[case_id]
        for short, _ in t12.CHANNELS:
            diagnostic = convergence_tail_diagnostics(
                [int(row["lmax"]) for row in case_rows],
                [float(row[f"{short}_successive_change"]) for row in case_rows],
                [row[f"{short}_change_applicable"] == "true" for row in case_rows],
                tolerance=TOLERANCE,
                maximum_order=MAX_LMAX,
            )
            rows.append({
                "case_id": case_id,
                "channel": short,
                "final_lmax": case_rows[-1]["lmax"],
                "confirmation_order": diagnostic.confirmation_order or 0,
                "last_change": _format(diagnostic.last_change),
                "classification": diagnostic.classification,
                "q_median_last_four": _format(diagnostic.q_median),
                "q_minimum_last_four": _format(diagnostic.q_minimum),
                "q_maximum_last_four": _format(diagnostic.q_maximum),
                "q_count": diagnostic.q_count,
                "oscillatory_last_five": _bool(diagnostic.oscillatory),
            })
    return rows


def _frozen_vectors() -> tuple[dict[str, dict[str, str]], dict[str, tuple[np.ndarray, np.ndarray]]]:
    cases, forces = t12._load_t08()
    values: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for case_id, rows in forces.items():
        a = np.asarray([[float(row["a_x"]), float(row["a_y"]), 0.0] for row in rows])
        d = np.asarray([[float(row["d_x"]), float(row["d_y"]), 0.0] for row in rows])
        values[case_id] = (a, d)
    return cases, values


def _final_states() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    old_grouped = _group(_read(T12_RAW))
    extended = _group(_read(EXTENDED_PATH))
    original_comparison = {row["case_id"]: row for row in _read(T12_COMPARISON)}
    cases, frozen = _frozen_vectors()
    resolved: list[dict[str, object]] = []
    mechanisms: list[dict[str, object]] = []
    for order, sentinel in enumerate(t12.SENTINELS, start=1):
        rows = extended.get(sentinel.case_id, old_grouped[sentinel.case_id])
        final = rows[-1]
        interaction_confirmed = final["interaction_confirmed"] == "true"
        all_confirmed = all(final[f"{short}_confirmed"] == "true" for short, _ in t12.CHANNELS)
        a, d = frozen[sentinel.case_id]
        interaction = t12._deserialize(final["interaction_forces_xyz"])
        external = t12._deserialize(final["external_scattered_forces_xyz"])
        ss = t12._deserialize(final["scattered_scattered_forces_xyz"])
        comparison = t12.compare_model_e_forces(a[:, :2], d[:, :2], interaction, external, ss)
        source = original_comparison[sentinel.case_id]
        epsilon = comparison.epsilon_a_e if interaction_confirmed else 0.0
        predicted = t12.FROZEN_PREFACTOR * sentinel.rho_l1**t12.FROZEN_EXPONENT
        numerical_pass = (
            float(final["balanced_condition_number"]) < 10.0
            and float(final["balanced_backward_error"]) < 1.0e-12
            and float(final["effective_incident_closure_error"]) < 1.0e-12
            and float(final["scattering_closure_error"]) < 1.0e-12
            and float(final["force_decomposition_residual"]) < 1.0e-12
            and float(final["max_abs_fz"]) < 1.0e-12
        )
        resolved.append({
            "sentinel_order": order,
            "case_id": sentinel.case_id,
            "particle_count": sentinel.particle_count,
            "family": sentinel.family,
            "stratum": f"n{sentinel.particle_count}_{sentinel.family}",
            "rho_band": sentinel.rho_band,
            "f1": _format(sentinel.f1),
            "distance_ratio": _format(sentinel.distance_ratio),
            "eta": cases[sentinel.case_id]["eta"],
            "lambda_max": cases[sentinel.case_id]["lambda_max"],
            "rho_l1": _format(sentinel.rho_l1),
            "epsilon_a_t08": source["epsilon_a_t08"],
            "final_lmax_e": final["lmax"],
            "state_source": "t12_1" if sentinel.case_id in EXTENDED_CASES else "t12",
            "interaction_confirmed": _bool(interaction_confirmed),
            "all_channels_confirmed": _bool(all_confirmed),
            "epsilon_a_e": _format(epsilon),
            "epsilon_a_e_applicable": _bool(interaction_confirmed),
            "epsilon_a_e_reason": "applicable" if interaction_confirmed else "interaction_unconfirmed_at_21",
            "predicted_epsilon_a": _format(predicted),
            "prediction_factor": _format(
                max(epsilon / predicted, predicted / epsilon)
                if interaction_confirmed and epsilon > 0.0 else 0.0
            ),
            "frozen_log_residual": _format(
                np.log(epsilon) - np.log(predicted)
                if interaction_confirmed and epsilon > 0.0 else 0.0
            ),
            "diagnostics_pass": _bool(numerical_pass),
            "special_case": SPECIAL_CASES.get(sentinel.case_id, "none"),
            "historical_t12_gate": "NO-GO_T13",
            "t12_1_diagnostic_gate": "pending",
        })
        mechanism = mechanism_diagnostics(a, d, interaction, external, ss)
        mechanism_applicable = interaction_confirmed and all_confirmed
        record: dict[str, object] = {
            "sentinel_order": order,
            "case_id": sentinel.case_id,
            "particle_count": sentinel.particle_count,
            "family": sentinel.family,
            "rho_l1": _format(sentinel.rho_l1),
            "applicable": _bool(mechanism_applicable),
            "reason": "applicable" if mechanism_applicable else "all_relevant_channels_not_confirmed",
            "rms_c_d": _format(mechanism.rms_d if mechanism_applicable else 0.0),
            "rms_c_m": _format(mechanism.rms_m if mechanism_applicable else 0.0),
            "rms_c_s": _format(mechanism.rms_s if mechanism_applicable else 0.0),
            "rms_c": _format(mechanism.rms_c if mechanism_applicable else 0.0),
            "closure_rms": _format(mechanism.closure_rms if mechanism_applicable else 0.0),
        }
        for name in (
            "mu_dm", "mu_ds", "mu_ms", "mu_dc", "mu_mc", "mu_sc",
            "p_d", "p_m", "p_s", "projection_sum", "r_s_over_d",
            "r_m_over_d", "cancellation_ratio",
        ):
            value = getattr(mechanism, name)
            record[name] = _format(value.value if mechanism_applicable and value.applicable else 0.0)
            record[f"{name}_applicable"] = _bool(mechanism_applicable and value.applicable)
            record[f"{name}_reason"] = value.reason if mechanism_applicable else record["reason"]
        mechanisms.append(record)
    return resolved, mechanisms


def _candidate_predictions(
    resolved: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    eligible = [
        row for row in resolved
        if row["epsilon_a_e_applicable"] == "true" and float(row["epsilon_a_e"]) > 0.0
    ]
    observed = np.asarray([float(row["epsilon_a_e"]) for row in eligible])
    groups = [str(row["stratum"]) for row in eligible]
    predictions: dict[str, np.ndarray] = {
        "P0_frozen_rho_l1": np.asarray(
            [float(row["predicted_epsilon_a"]) for row in eligible]
        )
    }
    fold_parameters: list[dict[str, object]] = []
    for candidate, field in PREDICTORS.items():
        x = np.asarray([float(row[field]) for row in eligible])
        estimate = np.empty_like(observed)
        for group, train, test in leave_group_out_folds(groups):
            fit = fit_log_linear(x[train], observed[train])
            estimate[test] = fit.prefactor * x[test] ** fit.coefficient
            fold_parameters.append({
                "record_type": "fold_fit",
                "candidate": candidate,
                "scope": group,
                "point_count": fit.point_count,
                "prefactor": _format(fit.prefactor),
                "exponent": _format(fit.coefficient),
                "rmse_log": "0",
                "median_factor": "0",
                "p90_factor": "0",
                "maximum_factor": "0",
                "fraction_within_factor_2": "0",
                "spearman": "0",
                "reference_derived": _bool(candidate == "P4_epsilon_a_d"),
                "value": "0",
                "recommendation": "not_applicable",
            })
        predictions[candidate] = estimate

    oof_rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for candidate, estimate in predictions.items():
        metrics = out_of_fold_metrics(observed, estimate)
        diagnostics.append({
            "record_type": "global_oof",
            "candidate": candidate,
            "scope": "all",
            "point_count": metrics.point_count,
            "prefactor": _format(t12.FROZEN_PREFACTOR) if candidate == "P0_frozen_rho_l1" else "0",
            "exponent": _format(t12.FROZEN_EXPONENT) if candidate == "P0_frozen_rho_l1" else "0",
            "rmse_log": _format(metrics.rmse_log),
            "median_factor": _format(metrics.median_factor),
            "p90_factor": _format(metrics.p90_factor),
            "maximum_factor": _format(metrics.maximum_factor),
            "fraction_within_factor_2": _format(metrics.fraction_within_factor_2),
            "spearman": _format(metrics.spearman),
            "reference_derived": _bool(candidate == "P4_epsilon_a_d"),
            "value": "0",
            "recommendation": "not_applicable",
        })
        for row, prediction in zip(eligible, estimate):
            oof_rows.append({
                "sentinel_order": row["sentinel_order"],
                "case_id": row["case_id"],
                "stratum": row["stratum"],
                "candidate": candidate,
                "observed_epsilon_a_e": row["epsilon_a_e"],
                "predicted_epsilon_a_e": _format(prediction),
                "log_residual": _format(np.log(float(row["epsilon_a_e"])) - np.log(prediction)),
                "multiplicative_factor": _format(max(
                    float(row["epsilon_a_e"]) / prediction,
                    prediction / float(row["epsilon_a_e"]),
                )),
                "reference_derived": _bool(candidate == "P4_epsilon_a_d"),
            })
        for group in dict.fromkeys(groups):
            selection = np.asarray([item == group for item in groups])
            metrics_group = out_of_fold_metrics(observed[selection], estimate[selection])
            diagnostics.append({
                "record_type": "stratum_oof",
                "candidate": candidate,
                "scope": group,
                "point_count": metrics_group.point_count,
                "prefactor": "0",
                "exponent": "0",
                "rmse_log": _format(metrics_group.rmse_log),
                "median_factor": _format(metrics_group.median_factor),
                "p90_factor": _format(metrics_group.p90_factor),
                "maximum_factor": _format(metrics_group.maximum_factor),
                "fraction_within_factor_2": _format(metrics_group.fraction_within_factor_2),
                "spearman": _format(metrics_group.spearman),
                "reference_derived": _bool(candidate == "P4_epsilon_a_d"),
                "value": "0",
                "recommendation": "not_applicable",
            })
    diagnostics.extend(fold_parameters)
    for candidate, field in PREDICTORS.items():
        fit = fit_log_linear(
            [float(row[field]) for row in eligible],
            observed,
        )
        diagnostics.append({
            "record_type": "descriptive_global_fit",
            "candidate": candidate,
            "scope": "all",
            "point_count": fit.point_count,
            "prefactor": _format(fit.prefactor),
            "exponent": _format(fit.coefficient),
            "rmse_log": "0",
            "median_factor": "0",
            "p90_factor": "0",
            "maximum_factor": "0",
            "fraction_within_factor_2": "0",
            "spearman": "0",
            "reference_derived": _bool(candidate == "P4_epsilon_a_d"),
            "value": "0",
            "recommendation": "not_applicable",
        })
    return diagnostics, oof_rows


def _correlations(
    resolved: list[dict[str, object]],
    mechanisms: list[dict[str, object]],
) -> list[dict[str, object]]:
    mechanism_by_id = {row["case_id"]: row for row in mechanisms}
    selected = [
        row for row in resolved
        if row["epsilon_a_e_applicable"] == "true"
        and mechanism_by_id[str(row["case_id"])]["applicable"] == "true"
    ]
    residual = np.asarray([float(row["frozen_log_residual"]) for row in selected])
    fields = (
        ("particle_count", False),
        ("f1", False),
        ("distance_ratio", False),
        ("rho_l1", False),
        ("r_s_over_d", True),
        ("r_m_over_d", True),
        ("mu_ds", True),
        ("p_s", True),
        ("cancellation_ratio", True),
    )
    rows = []
    for field, derived in fields:
        values = np.asarray([
            float(mechanism_by_id[str(row["case_id"])][field])
            if derived else float(row[field])
            for row in selected
        ])
        rows.append({
            "record_type": "residual_spearman",
            "candidate": field,
            "scope": "all_mechanism_applicable",
            "point_count": len(selected),
            "prefactor": "0",
            "exponent": "0",
            "rmse_log": "0",
            "median_factor": "0",
            "p90_factor": "0",
            "maximum_factor": "0",
            "fraction_within_factor_2": "0",
            "spearman": _format(spearman_correlation(residual, values)),
            "reference_derived": _bool(derived),
            "value": "0",
            "recommendation": "not_applicable",
        })
    return rows


def _recommendation(
    resolved: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
) -> str:
    if any(row["interaction_confirmed"] != "true" for row in resolved):
        return "NEED_MORE_CONVERGENCE"
    if any(row["diagnostics_pass"] != "true" for row in resolved):
        return "NEED_MORE_CONVERGENCE"
    global_rows = {
        row["candidate"]: row for row in diagnostics if row["record_type"] == "global_oof"
    }
    p3 = global_rows["P3_rho_l1"]
    competitors = [global_rows[name] for name in ("P1_eta", "P2_lambda_max", "P3_rho_l1")]
    p3_rmse = float(p3["rmse_log"])
    compatible = (
        p3_rmse <= np.log(2.0)
        and float(p3["fraction_within_factor_2"]) >= 0.8
        and p3_rmse <= min(float(row["rmse_log"]) for row in competitors) + 0.05
    )
    if compatible:
        return "READY_T12_2_RHO1_RECALIBRATION_STUDY"
    if len(resolved) < 20:
        return "INCONCLUSIVE_SMALL_SENTINEL_SET"
    return "NEED_NEW_PHYSICS_INFORMED_DESCRIPTOR"


def _plot(
    extended: dict[str, list[dict[str, str]]],
    resolved: list[dict[str, object]],
    mechanisms: list[dict[str, object]],
    oof: list[dict[str, object]],
) -> None:
    colors = {"pair": "#4c78a8", "compact": "#54a24b", "irregular": "#b279a2", "linear": "#f58518"}
    markers = {2: "o", 3: "s", 4: "^"}
    resolved_by_id = {row["case_id"]: row for row in resolved}
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5), constrained_layout=True)
    for case_id in INTERACTION_UNCONFIRMED:
        rows = extended[case_id]
        axes[0, 0].plot(
            [int(row["lmax"]) for row in rows[1:]],
            [float(row["interaction_successive_change"]) for row in rows[1:]],
            lw=1.2,
            marker="o",
            ms=2.8,
            label=case_id,
        )
    axes[0, 0].axhline(TOLERANCE, color="black", ls="--", lw=1.1, label=r"$10^{-5}$")
    axes[0, 0].set(yscale="log", xlabel=r"$L_{\max}$", ylabel="successive interaction-force change", title="Extended convergence of six unresolved interactions")
    axes[0, 0].legend(fontsize=7, ncol=2)

    for candidate, color in (("P0_frozen_rho_l1", "#e45756"), ("P3_rho_l1", "#4c78a8")):
        selected = [row for row in oof if row["candidate"] == candidate]
        axes[0, 1].scatter(
            [float(row["predicted_epsilon_a_e"]) for row in selected],
            [float(row["observed_epsilon_a_e"]) for row in selected],
            s=28,
            color=color,
            alpha=0.8,
            label=candidate.split("_", 1)[0],
        )
    positive = [
        float(row[field]) for row in oof
        for field in ("predicted_epsilon_a_e", "observed_epsilon_a_e")
    ]
    grid = np.geomspace(min(positive) * 0.7, max(positive) * 1.4, 200)
    axes[0, 1].plot(grid, grid, "k-", lw=1.0, label="identity")
    axes[0, 1].fill_between(grid, grid / 2.0, grid * 2.0, color="0.85", label="factor 2")
    axes[0, 1].set(xscale="log", yscale="log", xlabel="out-of-fold prediction", ylabel=r"observed $\varepsilon_A^E$", title="Frozen P0 versus cross-validated P3")
    axes[0, 1].legend()

    for row in mechanisms:
        if row["r_s_over_d_applicable"] != "true":
            continue
        case = resolved_by_id[row["case_id"]]
        axes[1, 0].scatter(
            float(row["rho_l1"]),
            float(row["r_s_over_d"]),
            marker=markers[int(row["particle_count"])],
            color=colors[row["family"]],
            s=42,
        )
    axes[1, 0].set(xscale="log", yscale="log", xlabel=r"$\rho_1$", ylabel=r"$F_{\rm RMS}(C_S)/F_{\rm RMS}(C_D)$", title="Scattered–scattered amplitude diagnostic")

    for row in mechanisms:
        if row["mu_ds_applicable"] != "true":
            continue
        case = resolved_by_id[row["case_id"]]
        axes[1, 1].scatter(
            float(row["mu_ds"]),
            float(case["frozen_log_residual"]),
            marker=markers[int(row["particle_count"])],
            color=colors[row["family"]],
            s=42,
        )
    axes[1, 1].axhline(0.0, color="0.5", ls=":", lw=0.9)
    axes[1, 1].set(xlabel=r"$\mu_{DS}$", ylabel="frozen log residual", title="Residual versus D–SS alignment (preregistered)")
    for axis in axes.flat:
        axis.grid(True, alpha=0.22)
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def analyze() -> None:
    extended_rows = _read(EXTENDED_PATH)
    extended = _group(extended_rows)
    if tuple(extended) != EXTENDED_CASES:
        raise RuntimeError("extended convergence does not contain exactly the ten cases")
    for case_id, rows in extended.items():
        orders = [int(row["lmax"]) for row in rows]
        if orders[:12] != list(range(2, 14)) or orders != list(range(2, orders[-1] + 1)):
            raise RuntimeError(f"incomplete deterministic order sequence for {case_id}")
        if any(row["source"] != ("t12" if int(row["lmax"]) <= 13 else "t12_1") for row in rows):
            raise RuntimeError(f"invalid provenance for {case_id}")
    summary = _tail_summary(extended)
    resolved, mechanisms = _final_states()
    if len(resolved) != 28 or tuple(row["case_id"] for row in resolved) != tuple(item.case_id for item in t12.SENTINELS):
        raise RuntimeError("resolved comparison lost the frozen 28-case order")
    diagnostics, oof = _candidate_predictions(resolved)
    diagnostics.extend(_correlations(resolved, mechanisms))
    recommendation = _recommendation(resolved, diagnostics)
    for row in resolved:
        row["t12_1_diagnostic_gate"] = recommendation
    diagnostics.append({
        "record_type": "recommendation",
        "candidate": "t12_1",
        "scope": "all",
        "point_count": len(resolved),
        "prefactor": "0",
        "exponent": "0",
        "rmse_log": "0",
        "median_factor": "0",
        "p90_factor": "0",
        "maximum_factor": "0",
        "fraction_within_factor_2": "0",
        "spearman": "0",
        "reference_derived": "false",
        "value": "0",
        "recommendation": recommendation,
    })
    _atomic_write(SUMMARY_PATH, summary)
    _atomic_write(RESOLVED_PATH, resolved)
    _atomic_write(MECHANISM_PATH, mechanisms)
    _atomic_write(PREDICTOR_PATH, diagnostics)
    _atomic_write(OOF_PATH, oof)
    _plot(extended, resolved, mechanisms, oof)
    print(f"T12.1 recommendation: {recommendation}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="validate the raw extension and rebuild only derived artifacts",
    )
    arguments = parser.parse_args()
    if not arguments.analyze_only:
        run_extension()
    analyze()


if __name__ == "__main__":
    main()
