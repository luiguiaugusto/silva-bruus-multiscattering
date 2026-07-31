#!/usr/bin/env python3
"""Run and analyze the preregistered T12 Model-E sentinel campaign."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    cluster_family,
    compare_model_e_forces,
    nodal_pair_force_on_probe,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    solve_multipolar_nodal_interaction_forces,
)
from analyze_t11_model_e import _successive


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
T08_CASES = DATA / "t08_cases.csv"
T08_FORCES = DATA / "t08_forces.csv"
T08_FITS = DATA / "t08_predictor_fits.csv"
T08_THRESHOLDS = DATA / "t08_validity_thresholds.csv"
MANIFEST_PATH = DATA / "t12_sentinel_manifest.csv"
CONVERGENCE_PATH = DATA / "t12_model_e_convergence.csv"
COMPARISON_PATH = DATA / "t12_model_comparison.csv"
AUDIT_PATH = DATA / "t12_threshold_audit.csv"
FIGURE_PATH = FIGURES / "t12_model_e_sentinel_audit.png"

RADIUS = 1.0
ENERGY_DENSITY = 1.0
KA = 0.1
F0 = 0.0
CONVERGENCE_TOLERANCE = 1.0e-5
MAX_LMAX = 13
FROZEN_PREFACTOR = 2.6353684041458636
FROZEN_EXPONENT = 1.1088518115798773
FROZEN_THRESHOLDS = (
    (0.01, 0.0053990295322641655),
    (0.05, 0.02000077753569526),
    (0.10, 0.03914887870730305),
)
CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)


@dataclass(frozen=True)
class Sentinel:
    case_id: str
    particle_count: int
    family: str
    rho_band: int
    f1: float
    distance_ratio: float
    rho_l1: float
    reference_lmax: int


SENTINELS = (
    Sentinel("n2_pair_f1.0_d6.0", 2, "pair", 1, 1.0, 6.0, 0.002699514766132083, 5),
    Sentinel("n2_pair_f1.0_d4.0", 2, "pair", 2, 1.0, 4.0, 0.008414320011147664, 7),
    Sentinel("n2_pair_f0.8_d2.5", 2, "pair", 3, 0.8, 2.5, 0.026387876003953037, 9),
    Sentinel("n2_pair_f1.0_d2.1", 2, "pair", 4, 1.0, 2.1, 0.055167482766554234, 13),
    Sentinel("n3_compact_f0.4_d6.0", 3, "compact", 1, 0.4, 6.0, 0.0021596118129056665, 5),
    Sentinel("n3_compact_f0.1_d2.1", 3, "compact", 2, 0.1, 2.1, 0.011033496553310861, 9),
    Sentinel("n3_compact_f0.4_d2.5", 3, "compact", 3, 0.4, 2.5, 0.026387876003953033, 9),
    Sentinel("n3_compact_f0.8_d2.1", 3, "compact", 4, 0.8, 2.1, 0.08826797242648693, 13),
    Sentinel("n3_irregular_f0.1_d3.0", 3, "irregular", 1, 0.1, 3.0, 0.0028389581530360573, 7),
    Sentinel("n3_irregular_f0.8_d4.0", 3, "irregular", 2, 0.8, 4.0, 0.009938913727557815, 7),
    Sentinel("n3_irregular_f1.0_d3.0", 3, "irregular", 3, 1.0, 3.0, 0.02838958153036055, 7),
    Sentinel("n3_irregular_f1.0_d2.1", 3, "irregular", 4, 1.0, 2.1, 0.08067672477956526, 13),
    Sentinel("n3_linear_f0.1_d3.0", 3, "linear", 1, 0.1, 3.0, 0.002872347791196304, 7),
    Sentinel("n3_linear_f0.8_d4.0", 3, "linear", 2, 0.8, 4.0, 0.010030161281933956, 7),
    Sentinel("n3_linear_f1.0_d3.0", 3, "linear", 3, 1.0, 3.0, 0.02872347791196306, 7),
    Sentinel("n3_linear_f1.0_d2.1", 3, "linear", 4, 1.0, 2.1, 0.08176358624705106, 13),
    Sentinel("n4_compact_f0.4_d6.0", 4, "compact", 1, 0.4, 6.0, 0.002587643284047554, 5),
    Sentinel("n4_compact_f0.1_d2.1", 4, "compact", 2, 0.1, 2.1, 0.013024703240875866, 9),
    Sentinel("n4_compact_f0.4_d2.5", 4, "compact", 3, 0.4, 2.5, 0.031187722478002788, 9),
    Sentinel("n4_compact_f1.0_d2.5", 4, "compact", 4, 1.0, 2.5, 0.07796930619500694, 9),
    Sentinel("n4_irregular_f0.1_d3.0", 4, "irregular", 1, 0.1, 3.0, 0.003408578045201502, 7),
    Sentinel("n4_irregular_f0.1_d2.1", 4, "irregular", 2, 0.1, 2.1, 0.009652629379316477, 9),
    Sentinel("n4_irregular_f0.8_d3.0", 4, "irregular", 3, 0.8, 3.0, 0.027268624361612016, 7),
    Sentinel("n4_irregular_f0.8_d2.1", 4, "irregular", 4, 0.8, 2.1, 0.07722103503453182, 11),
    Sentinel("n4_linear_f0.1_d3.0", 4, "linear", 1, 0.1, 3.0, 0.0034000806868721262, 7),
    Sentinel("n4_linear_f0.1_d2.1", 4, "linear", 2, 0.1, 2.1, 0.009659874464446596, 9),
    Sentinel("n4_linear_f0.8_d3.0", 4, "linear", 3, 0.8, 3.0, 0.027200645494976996, 7),
    Sentinel("n4_linear_f0.8_d2.1", 4, "linear", 4, 0.8, 2.1, 0.07727899571557277, 13),
)


def _format(value: float) -> str:
    return format(float(value), ".17g")


def _boolean(value: bool) -> str:
    return str(bool(value)).lower()


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _atomic_write(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _serialize(vectors: np.ndarray) -> str:
    values = np.asarray(vectors, dtype=float)
    return ";".join(":".join(_format(value) for value in row) for row in values)


def _deserialize(value: str, columns: int = 3) -> np.ndarray:
    rows = [[float(item) for item in row.split(":")] for row in value.split(";")]
    result = np.asarray(rows, dtype=float)
    if result.ndim != 2 or result.shape[1] != columns or not np.all(np.isfinite(result)):
        raise ValueError("serialized vector field is invalid")
    return result


def _model_a(positions: np.ndarray, f1: float) -> np.ndarray:
    forces = np.zeros((len(positions), 2), dtype=float)
    for target in range(len(positions)):
        for source in range(len(positions)):
            if source != target:
                forces[target] += nodal_pair_force_on_probe(
                    positions[target, :2], positions[source, :2], KA, RADIUS,
                    ENERGY_DENSITY, f1,
                )
    return forces


def _load_t08() -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    cases = {row["case_id"]: row for row in _read(T08_CASES)}
    forces: dict[str, list[dict[str, str]]] = {}
    for row in _read(T08_FORCES):
        forces.setdefault(row["case_id"], []).append(row)
    for rows in forces.values():
        rows.sort(key=lambda row: int(row["particle_index"]))
    return cases, forces


def _frozen_fit_and_thresholds() -> None:
    selected = [
        row for row in _read(T08_FITS)
        if row["record_type"] == "calibration_fit"
        and row["predictor"] == "rho_l1"
        and row["response"] == "epsilon_a"
    ]
    if len(selected) != 1:
        raise RuntimeError("the frozen rho_l1 calibration row is missing or duplicated")
    row = selected[0]
    if float(row["prefactor"]) != FROZEN_PREFACTOR or float(row["exponent"]) != FROZEN_EXPONENT:
        raise RuntimeError("the frozen rho_l1 calibration changed")
    threshold_rows = _read(T08_THRESHOLDS)
    if len(threshold_rows) != 3:
        raise RuntimeError("the frozen threshold table must contain three rows")
    for row, (tolerance, threshold) in zip(threshold_rows, FROZEN_THRESHOLDS):
        if not np.isclose(float(row["tolerance"]), tolerance, rtol=0.0, atol=1e-16):
            raise RuntimeError("a frozen error tolerance changed")
        if float(row["threshold"]) != threshold:
            raise RuntimeError("a frozen rho_l1 threshold changed")


def validate_and_write_manifest() -> tuple[
    dict[str, dict[str, str]], dict[str, list[dict[str, str]]], dict[str, np.ndarray]
]:
    """Validate the preregistration against T08 and publish its manifest."""

    _frozen_fit_and_thresholds()
    cases, forces = _load_t08()
    positions_by_case: dict[str, np.ndarray] = {}
    manifest_rows: list[dict[str, object]] = []
    for order, sentinel in enumerate(SENTINELS, start=1):
        if sentinel.case_id not in cases or sentinel.case_id not in forces:
            raise RuntimeError(f"missing frozen T08 case {sentinel.case_id}")
        case = cases[sentinel.case_id]
        force_rows = forces[sentinel.case_id]
        positions = np.asarray(
            [[float(row[name]) for name in ("x", "y", "z")] for row in force_rows]
        )
        reconstructed = cluster_family(
            sentinel.particle_count, sentinel.family, sentinel.distance_ratio
        )
        checks = (
            case["split"] == "calibration",
            int(case["particle_count"]) == sentinel.particle_count <= 4,
            case["family"] == sentinel.family,
            float(case["ka"]) == KA,
            float(case["f0"]) == F0,
            float(case["f1"]) == sentinel.f1,
            float(case["distance_ratio"]) == sentinel.distance_ratio,
            float(case["rho_l1"]) == sentinel.rho_l1,
            int(case["reference_lmax"]) == sentinel.reference_lmax,
            len(force_rows) == sentinel.particle_count,
            np.all(np.isfinite(positions)),
            np.all(positions[:, 2] == 0.0),
            len(np.unique(positions, axis=0)) == sentinel.particle_count,
            np.allclose(positions, reconstructed, rtol=0.0, atol=5e-15),
        )
        if not all(checks):
            raise RuntimeError(f"T08 metadata or positions disagree for {sentinel.case_id}")
        positions_by_case[sentinel.case_id] = positions
        manifest_rows.append({
            "sentinel_order": order,
            "case_id": sentinel.case_id,
            "split": "calibration",
            "particle_count": sentinel.particle_count,
            "family": sentinel.family,
            "rho_band": sentinel.rho_band,
            "f0": _format(F0),
            "f1": _format(sentinel.f1),
            "ka": _format(KA),
            "distance_ratio": _format(sentinel.distance_ratio),
            "rho_l1": _format(sentinel.rho_l1),
            "reference_lmax": sentinel.reference_lmax,
            "source": "pre_registered_t08_calibration",
            "validated": "true",
        })
    _atomic_write(MANIFEST_PATH, list(manifest_rows[0]), manifest_rows)
    return cases, forces, positions_by_case


def audit_models_a_d(
    cases: dict[str, dict[str, str]],
    forces: dict[str, list[dict[str, str]]],
    positions_by_case: dict[str, np.ndarray],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Audit frozen A/D vectors against their public production APIs."""

    audited: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for sentinel in SENTINELS:
        rows = forces[sentinel.case_id]
        frozen_a = np.asarray([[float(row["a_x"]), float(row["a_y"])] for row in rows])
        frozen_d = np.asarray([[float(row["d_x"]), float(row["d_y"])] for row in rows])
        positions = positions_by_case[sentinel.case_id]
        recalculated_a = _model_a(positions, sentinel.f1)
        recalculated_d = solve_multipolar_nodal_interaction_forces(
            positions, KA, RADIUS, ENERGY_DENSITY, F0, sentinel.f1,
            int(cases[sentinel.case_id]["reference_lmax"]),
        ).forces_xy
        if not np.allclose(recalculated_a, frozen_a, rtol=5e-12, atol=5e-14):
            raise RuntimeError(f"Model A audit failed for {sentinel.case_id}")
        if not np.allclose(recalculated_d, frozen_d, rtol=5e-12, atol=5e-14):
            raise RuntimeError(f"Model D audit failed for {sentinel.case_id}")
        audited[sentinel.case_id] = (frozen_a, frozen_d)
    return audited


def _minimum_confirmation(results: list[object], attribute: str) -> int:
    for index in range(2, len(results)):
        first, first_applicable, _ = _successive(
            getattr(results[index - 1], attribute), getattr(results[index - 2], attribute)
        )
        second, second_applicable, _ = _successive(
            getattr(results[index], attribute), getattr(results[index - 1], attribute)
        )
        if (
            first_applicable and second_applicable
            and first <= CONVERGENCE_TOLERANCE
            and second <= CONVERGENCE_TOLERANCE
        ):
            return results[index].lmax
    return 0


def _all_channels_confirmed(results: list[object]) -> bool:
    return all(_minimum_confirmation(results, attribute) > 0 for _, attribute in CHANNELS)


def _campaign_rows(sentinel: Sentinel, positions: np.ndarray) -> list[dict[str, object]]:
    results = []
    for order in range(2, MAX_LMAX + 1):
        result = solve_model_e_nodal(
            positions, KA, RADIUS, ENERGY_DENSITY, F0, sentinel.f1, order
        )
        if result.solution.production_solver != "balanced_sqrt":
            raise RuntimeError("Model E did not use the balanced production solver")
        results.append(result)
        if order >= 5 and _all_channels_confirmed(results):
            break

    confirmation = {
        short: _minimum_confirmation(results, attribute)
        for short, attribute in CHANNELS
    }
    final_lmax = results[-1].lmax
    coordinates = _serialize(positions)
    rows: list[dict[str, object]] = []
    for index, result in enumerate(results):
        solution = result.solution
        previous = results[index - 1] if index else None
        row: dict[str, object] = {
            "case_id": sentinel.case_id,
            "particle_count": sentinel.particle_count,
            "family": sentinel.family,
            "rho_band": sentinel.rho_band,
            "f1": _format(sentinel.f1),
            "distance_ratio": _format(sentinel.distance_ratio),
            "rho_l1": _format(sentinel.rho_l1),
            "coordinates_xyz": coordinates,
            "lmax": result.lmax,
            "final_lmax": final_lmax,
            "system_dimension": solution.balanced_system_matrix.shape[0],
            "active_modes_per_particle": len(solution.active_modes),
            "balanced_condition_number": _format(solution.balanced_condition_number),
            "balanced_backward_error": _format(solution.balanced_backward_error),
            "effective_incident_closure_error": _format(solution.effective_incident_closure_error),
            "scattering_closure_error": _format(solution.scattering_closure_error),
            "force_decomposition_residual": _format(result.decomposition_residual),
            "max_abs_fz": _format(max(
                float(np.max(np.abs(getattr(result, attribute)[:, 2])))
                for _, attribute in CHANNELS
            )),
            "total_forces_xyz": _serialize(result.total_forces_xyz),
            "external_forces_xyz": _serialize(result.external_forces_xyz),
            "interaction_forces_xyz": _serialize(result.interaction_forces_xyz),
            "external_scattered_forces_xyz": _serialize(result.external_scattered_forces_xyz),
            "scattered_scattered_forces_xyz": _serialize(result.scattered_scattered_forces_xyz),
            "production_solver": solution.production_solver,
            "finite": "true",
            "campaign_complete": "true",
        }
        for short, attribute in CHANNELS:
            current = getattr(result, attribute)
            if previous is None:
                change, applicable, absolute = 0.0, False, 0.0
            else:
                change, applicable, absolute = _successive(current, getattr(previous, attribute))
            row[f"{short}_rms"] = _format(rms_vector_magnitude_xyz(current))
            row[f"{short}_successive_change"] = _format(change)
            row[f"{short}_absolute_change"] = _format(absolute)
            row[f"{short}_change_applicable"] = _boolean(applicable)
            row[f"{short}_minimum_confirmed_lmax"] = confirmation[short]
            row[f"{short}_confirmed"] = _boolean(confirmation[short] > 0)
        numeric_values = [
            float(value) for key, value in row.items()
            if key.endswith(("_rms", "_change", "_error", "_residual", "_number"))
        ]
        if not np.all(np.isfinite(numeric_values)):
            raise RuntimeError(f"non-finite convergence diagnostic for {sentinel.case_id}")
        rows.append(row)
    return rows


def run_campaign() -> None:
    cases, forces, positions_by_case = validate_and_write_manifest()
    audit_models_a_d(cases, forces, positions_by_case)
    rows: list[dict[str, object]] = []
    for index, sentinel in enumerate(SENTINELS, start=1):
        print(f"Model E sentinel {index:02d}/28: {sentinel.case_id}", flush=True)
        rows.extend(_campaign_rows(sentinel, positions_by_case[sentinel.case_id]))
    _atomic_write(CONVERGENCE_PATH, list(rows[0]), rows)


def _validate_complete_convergence(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
    if tuple(grouped) != tuple(item.case_id for item in SENTINELS):
        raise RuntimeError("convergence CSV does not contain the complete preregistration")
    for sentinel in SENTINELS:
        case_rows = grouped[sentinel.case_id]
        orders = [int(row["lmax"]) for row in case_rows]
        if orders != list(range(2, int(case_rows[-1]["final_lmax"]) + 1)):
            raise RuntimeError(f"incomplete order sequence for {sentinel.case_id}")
        if orders[-1] < 5 or orders[-1] > MAX_LMAX:
            raise RuntimeError(f"invalid final order for {sentinel.case_id}")
        if any(row["campaign_complete"] != "true" for row in case_rows):
            raise RuntimeError("convergence campaign is not marked complete")
    return grouped


def _rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def _performance(rows: list[dict[str, object]]) -> dict[str, float]:
    eligible = [row for row in rows if row["prediction_metric_applicable"] == "true"]
    observed = np.asarray([float(row["epsilon_a_e"]) for row in eligible])
    predicted = np.asarray([float(row["predicted_epsilon_a"]) for row in eligible])
    log_residual = np.log(observed) - np.log(predicted)
    factors = np.maximum(observed / predicted, predicted / observed)
    spearman = float(np.corrcoef(_rank(observed), _rank(predicted))[0, 1])
    return {
        "point_count": float(len(eligible)),
        "rmse_log": float(np.sqrt(np.mean(log_residual**2))),
        "median_factor": float(np.median(factors)),
        "p90_factor": float(np.percentile(factors, 90)),
        "maximum_factor": float(np.max(factors)),
        "fraction_within_factor_2": float(np.mean(factors <= 2.0)),
        "spearman": spearman,
    }


def _audit_row(
    comparison_rows: list[dict[str, object]], tolerance: float, scope: str,
    selector,
) -> dict[str, object]:
    selected = [
        row for row in comparison_rows
        if selector(row) and row["threshold_metric_applicable"] == "true"
    ]
    predicted = [row for row in selected if row[f"predicted_safe_{int(tolerance * 100)}pct"] == "true"]
    observed = [row for row in selected if float(row["epsilon_a_e"]) <= tolerance]
    false_safe = [row for row in predicted if float(row["epsilon_a_e"]) > tolerance]
    false_unsafe = [
        row for row in selected
        if row[f"predicted_safe_{int(tolerance * 100)}pct"] == "false"
        and float(row["epsilon_a_e"]) <= tolerance
    ]
    return {
        "record_type": "threshold",
        "scope": scope,
        "tolerance": _format(tolerance),
        "eligible_count": len(selected),
        "predicted_safe_count": len(predicted),
        "observed_safe_count": len(observed),
        "false_safe_count": len(false_safe),
        "false_unsafe_count": len(false_unsafe),
        "worst_predicted_safe_error": _format(max((float(row["epsilon_a_e"]) for row in predicted), default=0.0)),
        "max_rho_observed_safe": _format(max((float(row["rho_l1"]) for row in observed), default=0.0)),
        "point_count": 0,
        "rmse_log": "0",
        "median_factor": "0",
        "p90_factor": "0",
        "maximum_factor": "0",
        "fraction_within_factor_2": "0",
        "spearman": "0",
        "t12_gate_supported": "false",
        "gate_recommendation": "not_applicable",
        "gate_failed_conditions": "none",
    }


def _write_audit(comparison_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], bool]:
    rows: list[dict[str, object]] = []
    for tolerance, _ in FROZEN_THRESHOLDS:
        rows.append(_audit_row(comparison_rows, tolerance, "all", lambda row: True))
        for particle_count in (2, 3, 4):
            rows.append(_audit_row(
                comparison_rows, tolerance, f"n{particle_count}",
                lambda row, n=particle_count: int(row["particle_count"]) == n,
            ))
        for family in ("pair", "compact", "irregular", "linear"):
            rows.append(_audit_row(
                comparison_rows, tolerance, f"family_{family}",
                lambda row, name=family: row["family"] == name,
            ))

    performance = _performance(comparison_rows)
    performance_row = {
        "record_type": "frozen_prediction_performance",
        "scope": "all_applicable",
        "tolerance": "0",
        "eligible_count": 0,
        "predicted_safe_count": 0,
        "observed_safe_count": 0,
        "false_safe_count": 0,
        "false_unsafe_count": 0,
        "worst_predicted_safe_error": "0",
        "max_rho_observed_safe": "0",
        **{key: _format(value) if key != "point_count" else int(value) for key, value in performance.items()},
        "t12_gate_supported": "false",
        "gate_recommendation": "not_applicable",
        "gate_failed_conditions": "none",
    }
    rows.append(performance_row)

    interaction_coverage = np.mean([
        row["interaction_confirmed"] == "true" for row in comparison_rows
    ])
    diagnostics_pass = all(row["diagnostics_pass"] == "true" for row in comparison_rows)
    five_percent = next(
        row for row in rows
        if row["record_type"] == "threshold"
        and row["scope"] == "all"
        and float(row["tolerance"]) == 0.05
    )
    conditions = {
        "interaction_coverage": interaction_coverage >= 0.8,
        "diagnostics": diagnostics_pass,
        "rmse_log": performance["rmse_log"] <= np.log(2.0),
        "factor_two": performance["fraction_within_factor_2"] >= 0.8,
        "five_percent_coverage": int(five_percent["predicted_safe_count"]) >= 1,
        "five_percent_no_false_safe": int(five_percent["false_safe_count"]) == 0,
    }
    gate_supported = all(conditions.values())
    failed = [name for name, passed in conditions.items() if not passed]
    gate_row = {
        "record_type": "gate",
        "scope": "t12",
        "tolerance": "0",
        "eligible_count": 0,
        "predicted_safe_count": 0,
        "observed_safe_count": 0,
        "false_safe_count": 0,
        "false_unsafe_count": 0,
        "worst_predicted_safe_error": "0",
        "max_rho_observed_safe": "0",
        "point_count": int(performance["point_count"]),
        "rmse_log": _format(performance["rmse_log"]),
        "median_factor": _format(performance["median_factor"]),
        "p90_factor": _format(performance["p90_factor"]),
        "maximum_factor": _format(performance["maximum_factor"]),
        "fraction_within_factor_2": _format(performance["fraction_within_factor_2"]),
        "spearman": _format(performance["spearman"]),
        "t12_gate_supported": _boolean(gate_supported),
        "gate_recommendation": "GO_T13" if gate_supported else "NO-GO_T13",
        "gate_failed_conditions": ";".join(failed) if failed else "none",
    }
    rows.append(gate_row)
    _atomic_write(AUDIT_PATH, list(rows[0]), rows)
    return rows, gate_supported


def _plot(comparison_rows: list[dict[str, object]]) -> None:
    colors = {"pair": "#4c78a8", "linear": "#f58518", "compact": "#54a24b", "irregular": "#b279a2"}
    markers = {2: "o", 3: "s", 4: "^"}
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), constrained_layout=True)
    rho_grid = np.geomspace(
        min(float(row["rho_l1"]) for row in comparison_rows) * 0.8,
        max(float(row["rho_l1"]) for row in comparison_rows) * 1.2,
        300,
    )

    for row in comparison_rows:
        rho = float(row["rho_l1"])
        observed = float(row["epsilon_a_e"])
        applicable = row["threshold_metric_applicable"] == "true" and observed > 0.0
        if observed > 0.0:
            axes[0, 0].scatter(
                rho, observed,
                marker=markers[int(row["particle_count"])],
                facecolors=colors[row["family"]] if applicable else "none",
                edgecolors=colors[row["family"]], s=48, linewidths=1.2,
            )
        if applicable:
            predicted = float(row["predicted_epsilon_a"])
            axes[0, 1].scatter(
                predicted, observed, marker=markers[int(row["particle_count"])],
                color=colors[row["family"]], s=48,
            )
            if row["mechanism_decomposition_applicable"] == "true":
                axes[1, 0].scatter(rho, float(row["x_d_minus_a"]), marker="o", color="#4c78a8", s=25)
                axes[1, 0].scatter(rho, float(row["x_mie_external"]), marker="s", color="#e45756", s=25)
                axes[1, 0].scatter(rho, float(row["x_scattered_scattered"]), marker="^", color="#54a24b", s=25)
            epsilon_t08 = float(row["epsilon_a_t08"])
            axes[1, 1].scatter(
                epsilon_t08, observed,
                marker=markers[int(row["particle_count"])],
                color="#d62728" if row["false_safe_5pct"] == "true" else colors[row["family"]],
                s=55,
            )

    axes[0, 0].plot(rho_grid, FROZEN_PREFACTOR * rho_grid**FROZEN_EXPONENT, "k-", lw=1.5, label="T08 frozen fit")
    for tolerance, threshold in FROZEN_THRESHOLDS:
        axes[0, 0].axvline(threshold, color="0.65", ls=":", lw=0.9)
        axes[0, 0].axhline(tolerance, color="0.65", ls="--", lw=0.9)
    axes[0, 0].set(xscale="log", yscale="log", xlabel=r"$\rho_1$", ylabel=r"$\varepsilon_A^E$", title="Model E sentinel audit")
    axes[0, 0].legend(loc="best")

    positive = [
        value for row in comparison_rows
        for value in (float(row["predicted_epsilon_a"]), float(row["epsilon_a_e"]))
        if value > 0.0 and row["prediction_metric_applicable"] == "true"
    ]
    lower, upper = min(positive) * 0.75, max(positive) * 1.3
    diagonal = np.geomspace(lower, upper, 200)
    axes[0, 1].plot(diagonal, diagonal, "k-", lw=1.3, label="identity")
    axes[0, 1].fill_between(diagonal, diagonal / 2.0, diagonal * 2.0, color="0.85", alpha=0.7, label="factor 2")
    axes[0, 1].set(xscale="log", yscale="log", xlim=(lower, upper), ylim=(lower, upper), xlabel="frozen prediction", ylabel="observed Model E error", title="Observed versus frozen prediction")
    axes[0, 1].legend(loc="best")

    axes[1, 0].set(xscale="log", yscale="log", xlabel=r"$\rho_1$", ylabel="normalized RMS amplitude", title="A–D–E mechanism amplitudes")
    for marker, color, label in (("o", "#4c78a8", r"$X_{D-A}$"), ("s", "#e45756", r"$X_{\mathrm{Mie/ext-sc}}$"), ("^", "#54a24b", r"$X_{\mathrm{ss}}$")):
        axes[1, 0].scatter([], [], marker=marker, color=color, label=label)
    axes[1, 0].legend(loc="best")

    for tolerance, _ in FROZEN_THRESHOLDS:
        axes[1, 1].axhline(tolerance, color="0.7", ls="--", lw=0.8)
        axes[1, 1].axvline(tolerance, color="0.7", ls=":", lw=0.8)
    comparison_positive = [
        value for row in comparison_rows
        for value in (float(row["epsilon_a_t08"]), float(row["epsilon_a_e"]))
        if value > 0.0 and row["threshold_metric_applicable"] == "true"
    ]
    lower_compare, upper_compare = min(comparison_positive) * 0.7, max(comparison_positive) * 1.4
    compare_line = np.geomspace(lower_compare, upper_compare, 200)
    axes[1, 1].plot(compare_line, compare_line, "k-", lw=1.3, label="identity")
    axes[1, 1].set(xscale="log", yscale="log", xlim=(lower_compare, upper_compare), ylim=(lower_compare, upper_compare), xlabel=r"T08 $\varepsilon_A^D$", ylabel=r"T12 $\varepsilon_A^E$", title="Frozen D reference versus full E reference")
    axes[1, 1].legend(loc="best")
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.2)
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def analyze_only() -> bool:
    cases, forces, positions_by_case = validate_and_write_manifest()
    audited = audit_models_a_d(cases, forces, positions_by_case)
    grouped = _validate_complete_convergence(_read(CONVERGENCE_PATH))
    comparison_rows: list[dict[str, object]] = []
    threshold_map = {int(tolerance * 100): threshold for tolerance, threshold in FROZEN_THRESHOLDS}
    for sentinel in SENTINELS:
        case = cases[sentinel.case_id]
        final = grouped[sentinel.case_id][-1]
        a, d = audited[sentinel.case_id]
        total = _deserialize(final["total_forces_xyz"])
        external = _deserialize(final["external_forces_xyz"])
        interaction = _deserialize(final["interaction_forces_xyz"])
        external_scattered = _deserialize(final["external_scattered_forces_xyz"])
        scattered_scattered = _deserialize(final["scattered_scattered_forces_xyz"])
        comparison = compare_model_e_forces(a, d, interaction, external_scattered, scattered_scattered)
        channel_confirmed = {
            short: final[f"{short}_confirmed"] == "true" for short, _ in CHANNELS
        }
        threshold_metric_applicable = (
            channel_confirmed["interaction"] and comparison.epsilon_a_e_applicable
        )
        mechanism_applicable = (
            channel_confirmed["interaction"]
            and channel_confirmed["external_scattered"]
            and channel_confirmed["scattered_scattered"]
        )
        prediction_applicable = threshold_metric_applicable and comparison.epsilon_a_e > 0.0
        predicted = FROZEN_PREFACTOR * sentinel.rho_l1**FROZEN_EXPONENT
        factor = (
            max(comparison.epsilon_a_e / predicted, predicted / comparison.epsilon_a_e)
            if prediction_applicable else 0.0
        )
        epsilon_t08 = float(case["epsilon_a"])
        ratio_t08 = (
            comparison.epsilon_a_e / epsilon_t08
            if threshold_metric_applicable and epsilon_t08 > 0.0 else 0.0
        )
        row: dict[str, object] = {
            "case_id": sentinel.case_id,
            "particle_count": sentinel.particle_count,
            "family": sentinel.family,
            "rho_band": sentinel.rho_band,
            "f1": _format(sentinel.f1),
            "distance_ratio": _format(sentinel.distance_ratio),
            "rho_l1": _format(sentinel.rho_l1),
            "reference_lmax_d": sentinel.reference_lmax,
            "final_lmax_e": int(final["final_lmax"]),
            **{f"{short}_confirmed": _boolean(value) for short, value in channel_confirmed.items()},
            "convergence_status": "confirmed" if all(channel_confirmed.values()) else "unconfirmed",
            "threshold_metric_applicable": _boolean(threshold_metric_applicable),
            "threshold_metric_reason": "applicable" if threshold_metric_applicable else ("unconfirmed" if not channel_confirmed["interaction"] else "reference_unresolved"),
            "mechanism_decomposition_applicable": _boolean(mechanism_applicable),
            "mechanism_decomposition_reason": "applicable" if mechanism_applicable else "unconfirmed",
            "prediction_metric_applicable": _boolean(prediction_applicable),
            "prediction_metric_reason": "applicable" if prediction_applicable else ("zero_observed_error" if threshold_metric_applicable else "unconfirmed_or_unresolved"),
            "rms_model_a": _format(rms_vector_magnitude_xyz(np.column_stack((a, np.zeros(len(a)))))),
            "rms_model_d": _format(rms_vector_magnitude_xyz(np.column_stack((d, np.zeros(len(d)))))),
            "rms_model_e_total": _format(rms_vector_magnitude_xyz(total)),
            "rms_model_e_external": _format(rms_vector_magnitude_xyz(external)),
            "rms_model_e_interaction": _format(comparison.rms_model_e_interaction),
            "rms_model_e_external_scattered": _format(comparison.rms_model_e_external_scattered),
            "rms_model_e_scattered_scattered": _format(comparison.rms_model_e_scattered_scattered),
            "epsilon_a_t08": _format(epsilon_t08),
            "epsilon_a_e": _format(comparison.epsilon_a_e),
            "epsilon_d_e": _format(comparison.epsilon_d_e),
            "epsilon_a_e_external_scattered": _format(comparison.epsilon_a_external_scattered),
            "epsilon_d_e_external_scattered": _format(comparison.epsilon_d_external_scattered),
            "symmetric_rms_a_e": _format(comparison.symmetric_a_e),
            "x_d_minus_a": _format(comparison.x_d_minus_a),
            "x_mie_external": _format(comparison.x_mie_external),
            "x_scattered_scattered": _format(comparison.x_scattered_scattered),
            "cancellation_ratio": _format(comparison.cancellation_ratio),
            "cancellation_ratio_applicable": _boolean(comparison.cancellation_ratio_applicable),
            "decomposition_max_abs_error": _format(comparison.decomposition_max_abs_error),
            "decomposition_relative_error": _format(comparison.decomposition_relative_error),
            "predicted_epsilon_a": _format(predicted),
            "prediction_factor": _format(factor),
            "epsilon_a_e_over_epsilon_a_t08": _format(ratio_t08),
            "ratio_t08_applicable": _boolean(threshold_metric_applicable and epsilon_t08 > 0.0),
            "final_truncation_uncertainty": final["interaction_successive_change"],
            "max_abs_interaction_fz": _format(comparison.max_abs_interaction_fz),
        }
        for percent, threshold in threshold_map.items():
            tolerance = percent / 100.0
            predicted_safe = sentinel.rho_l1 <= threshold
            observed_safe = threshold_metric_applicable and comparison.epsilon_a_e <= tolerance
            row[f"predicted_safe_{percent}pct"] = _boolean(predicted_safe)
            row[f"observed_safe_{percent}pct"] = _boolean(observed_safe)
            row[f"false_safe_{percent}pct"] = _boolean(
                threshold_metric_applicable and predicted_safe and not observed_safe
            )
        rms_scale = max(
            comparison.rms_model_e_interaction,
            comparison.rms_model_e_external_scattered,
            comparison.rms_model_e_scattered_scattered,
        )
        fz_pass = comparison.max_abs_interaction_fz <= 1e-11 * rms_scale if rms_scale > 0.0 else comparison.max_abs_interaction_fz == 0.0
        diagnostics_pass = (
            float(final["balanced_condition_number"]) < 10.0
            and float(final["balanced_backward_error"]) < 1e-12
            and float(final["effective_incident_closure_error"]) < 1e-12
            and float(final["scattering_closure_error"]) < 1e-12
            and float(final["force_decomposition_residual"]) < 1e-12
            and comparison.decomposition_relative_error < 1e-12
            and fz_pass
        )
        row["diagnostics_pass"] = _boolean(diagnostics_pass)
        row["finite"] = _boolean(all(
            np.isfinite(float(value)) for key, value in row.items()
            if key not in {
                "case_id", "family", "convergence_status", "threshold_metric_reason",
                "mechanism_decomposition_reason", "prediction_metric_reason",
            } and value not in {"true", "false", "confirmed", "unconfirmed"}
        ))
        comparison_rows.append(row)

    _atomic_write(COMPARISON_PATH, list(comparison_rows[0]), comparison_rows)
    _, gate_supported = _write_audit(comparison_rows)
    _plot(comparison_rows)
    return gate_supported


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analyze-only", action="store_true",
        help="reuse a validated complete convergence CSV without Model-E solves",
    )
    arguments = parser.parse_args()
    if not arguments.analyze_only:
        run_campaign()
    gate_supported = analyze_only()
    print(f"T12 gate supported: {str(gate_supported).lower()}")


if __name__ == "__main__":
    main()
