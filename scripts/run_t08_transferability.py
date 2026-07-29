#!/usr/bin/env python3
"""Run the single expensive T08 force sweep or audit a stratified sample."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np

from acoustic_ms.cluster_families import TransferabilityConfiguration, enumerate_transferability_configurations
from acoustic_ms.metrics import rms_vector_magnitude
from acoustic_ms.model_d import solve_multipolar_nodal_interaction_forces
from acoustic_ms.scaling import coupling_eta, maximum_geometric_coupling
from acoustic_ms.silva_bruus import nodal_pair_force_on_probe
from acoustic_ms.transferability import (
    matched_multipolar_pairwise_baseline,
    normalized_rms_difference,
    spectral_radius_l1,
    two_step_converged,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
K = 0.1
RADIUS = 1.0
ENERGY = 1.0
F0 = 0.0
TOLERANCE = 1e-3


def _model_a(positions, f1):
    forces = np.zeros((len(positions), 2))
    for first, second in combinations(range(len(positions)), 2):
        forces[first] += nodal_pair_force_on_probe(
            positions[first, :2], positions[second, :2], K, RADIUS, ENERGY, f1
        )
        forces[second] += nodal_pair_force_on_probe(
            positions[second, :2], positions[first, :2], K, RADIUS, ENERGY, f1
        )
    return forces


def _successive(current, previous, reference):
    if previous is None:
        return 0.0, False
    return normalized_rms_difference(current, previous, reference)


def _last_two(values):
    applicable = [value for value, flag in values if flag]
    penultimate = applicable[-2] if len(applicable) >= 2 else 0.0
    last = applicable[-1] if applicable else 0.0
    return penultimate, last


def _evaluate(configuration: TransferabilityConfiguration):
    positions = configuration.positions_xyz
    model_a = _model_a(positions, configuration.f1)
    history = []
    previous_d = previous_b = previous_r = None
    difference_d, difference_b, difference_r = [], [], []
    d_l1 = None
    for lmax in (1, 3, 5, 7, 9, 11):
        model_d = solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY, F0, configuration.f1, lmax
        )
        matched = matched_multipolar_pairwise_baseline(
            positions, K, RADIUS, ENERGY, F0, configuration.f1, lmax
        )
        residual = model_d.forces_xy - matched.forces_xy
        delta_d = _successive(model_d.forces_xy, previous_d, model_d.forces_xy)
        delta_b = _successive(matched.forces_xy, previous_b, model_d.forces_xy)
        delta_r = _successive(residual, previous_r, model_d.forces_xy)
        difference_d.append(delta_d); difference_b.append(delta_b); difference_r.append(delta_r)
        history.append((lmax, model_d, matched, residual, delta_d, delta_b, delta_r))
        if lmax == 1:
            d_l1 = model_d
        previous_d = model_d.forces_xy.copy(); previous_b = matched.forces_xy.copy(); previous_r = residual.copy()
        if lmax >= 5 and all(two_step_converged([x[0] for x in values if x[1]]) for values in (difference_d, difference_b, difference_r)):
            break
    if configuration.particle_count <= 4 and not all(
        two_step_converged([x[0] for x in values if x[1]])
        for values in (difference_d, difference_b, difference_r)
    ):
        lmax = 13
        model_d = solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY, F0, configuration.f1, lmax
        )
        matched = matched_multipolar_pairwise_baseline(
            positions, K, RADIUS, ENERGY, F0, configuration.f1, lmax
        )
        residual = model_d.forces_xy - matched.forces_xy
        delta_d = _successive(model_d.forces_xy, previous_d, model_d.forces_xy)
        delta_b = _successive(matched.forces_xy, previous_b, model_d.forces_xy)
        delta_r = _successive(residual, previous_r, model_d.forces_xy)
        difference_d.append(delta_d); difference_b.append(delta_b); difference_r.append(delta_r)
        history.append((lmax, model_d, matched, residual, delta_d, delta_b, delta_r))

    total_converged = two_step_converged([x[0] for x in difference_d if x[1]])
    pair_converged = two_step_converged([x[0] for x in difference_b if x[1]])
    residual_converged = two_step_converged([x[0] for x in difference_r if x[1]])
    joint_converged = total_converged and pair_converged and residual_converged
    lmax, model_d, matched, residual, *_ = history[-1]
    denominator = rms_vector_magnitude(model_d.forces_xy)
    metric_applicable = denominator > 128 * np.finfo(float).eps * max(
        denominator, rms_vector_magnitude(model_a), rms_vector_magnitude(matched.forces_xy)
    )
    def metric(first, second):
        return rms_vector_magnitude(first - second) / denominator if metric_applicable else 0.0
    epsilon_a = metric(model_a, model_d.forces_xy)
    epsilon_b = metric(matched.forces_xy, model_d.forces_xy)
    y_two = metric(matched.forces_xy, model_a)
    y_collective = metric(model_d.forces_xy, matched.forces_xy)
    y_multipolar = metric(model_d.forces_xy, d_l1.forces_xy)
    _, delta_d_last = _last_two(difference_d)
    delta_r_previous, delta_r_last = _last_two(difference_r)
    uncertainty_r = max(delta_r_previous, delta_r_last)
    residual_resolved = bool(joint_converged and epsilon_b > 5 * uncertainty_r)
    identity_error = float(np.max(np.abs(
        (model_d.forces_xy - model_a)
        - ((matched.forces_xy - model_a) + (model_d.forces_xy - matched.forces_xy))
    )))
    case = {
        "case_id": configuration.case_id, "split": configuration.split,
        "particle_count": configuration.particle_count, "family": configuration.family,
        "radius": RADIUS, "k": K, "ka": K * RADIUS, "energy_density": ENERGY,
        "f0": F0, "f1": configuration.f1, "distance_ratio": configuration.distance_ratio,
        "pair_count": configuration.particle_count * (configuration.particle_count - 1) // 2,
        "reference_lmax": lmax, "maximum_allowed_lmax": 13 if configuration.particle_count <= 4 else 11,
        "total_converged": str(total_converged).lower(),
        "matched_pairwise_converged": str(pair_converged).lower(),
        "collective_residual_converged": str(residual_converged).lower(),
        "joint_converged": str(joint_converged).lower(),
        "collective_residual_resolved": str(residual_resolved).lower(),
        "metric_applicable": str(metric_applicable).lower(),
        "eta": coupling_eta(RADIUS, configuration.distance_ratio, configuration.f1),
        "lambda_max": maximum_geometric_coupling(positions, RADIUS, configuration.f1),
        "rho_l1": spectral_radius_l1(d_l1),
        "rms_a": rms_vector_magnitude(model_a),
        "rms_matched_pairwise": rms_vector_magnitude(matched.forces_xy),
        "rms_d": denominator, "epsilon_a": epsilon_a, "epsilon_b": epsilon_b,
        "y_two_body": y_two, "y_collective": y_collective, "y_multipolar": y_multipolar,
        "delta_d_penultimate": _last_two(difference_d)[0], "delta_d_last": delta_d_last,
        "delta_b_penultimate": _last_two(difference_b)[0], "delta_b_last": _last_two(difference_b)[1],
        "delta_r_penultimate": delta_r_previous, "delta_r_last": delta_r_last,
        "collective_uncertainty": uncertainty_r, "identity_max_abs_error": identity_error,
        "maximum_physical_residual": max(max(item[1].solution.residual_relative, item[2].maximum_residual) for item in history),
        "maximum_balanced_condition": max(max(item[1].solution.condition_number, item[2].maximum_balanced_condition) for item in history),
        "maximum_raw_condition": max(max(item[1].solution.physical_condition_number, item[2].maximum_raw_condition) for item in history),
    }
    force_rows = []
    for index, position in enumerate(positions):
        force_rows.append({
            "case_id": configuration.case_id, "particle_index": index,
            "x": position[0], "y": position[1], "z": position[2],
            "a_x": model_a[index, 0], "a_y": model_a[index, 1],
            "b_l_x": matched.forces_xy[index, 0], "b_l_y": matched.forces_xy[index, 1],
            "d_x": model_d.forces_xy[index, 0], "d_y": model_d.forces_xy[index, 1],
            "collective_x": residual[index, 0], "collective_y": residual[index, 1],
        })
    convergence_rows = []
    for order, d_result, b_result, residual_value, dd, db, dr in history:
        convergence_rows.append({
            "case_id": configuration.case_id, "lmax": order,
            "rms_d": rms_vector_magnitude(d_result.forces_xy),
            "rms_matched_pairwise": rms_vector_magnitude(b_result.forces_xy),
            "rms_collective_residual": rms_vector_magnitude(residual_value),
            "delta_d": dd[0], "delta_d_applicable": str(dd[1]).lower(),
            "delta_b": db[0], "delta_b_applicable": str(db[1]).lower(),
            "delta_r": dr[0], "delta_r_applicable": str(dr[1]).lower(),
            "physical_residual": max(d_result.solution.residual_relative, b_result.maximum_residual),
            "balanced_condition": max(d_result.solution.condition_number, b_result.maximum_balanced_condition),
            "raw_condition": max(d_result.solution.physical_condition_number, b_result.maximum_raw_condition),
        })
    return case, force_rows, convergence_rows


def _format_rows(rows):
    return [{key: format(value, ".17g") if isinstance(value, float) else value for key, value in row.items()} for row in rows]


def _write(path, rows):
    rows = _format_rows(rows)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _evaluate_group(configurations):
    """Evaluate four contrasts for one geometry/distance with shared caches."""
    return [_evaluate(configuration) for configuration in configurations]


def _run(configurations, workers):
    if workers == 1:
        return [_evaluate(configuration) for configuration in configurations]
    groups = []
    for start in range(0, len(configurations), 24):
        geometry = configurations[start:start + 24]
        for distance_index in range(6):
            groups.append(tuple(geometry[contrast * 6 + distance_index] for contrast in range(4)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        grouped_results = list(pool.map(_evaluate_group, groups, chunksize=1))
    by_case = {result[0]["case_id"]: result for group in grouped_results for result in group}
    return [by_case[configuration.case_id] for configuration in configurations]

def _audit():
    configurations = {item.case_id: item for item in enumerate_transferability_configurations()}
    selected = [
        "n2_pair_f0.1_d10.0", "n4_irregular_f1.0_d2.1",
        "n10_linear_f0.1_d10.0", "n10_linear_f1.0_d2.1",
        "n10_compact_f0.1_d10.0", "n10_compact_f1.0_d2.1",
        "n10_irregular_f0.1_d10.0", "n10_irregular_f1.0_d2.1",
    ]
    existing = {row["case_id"]: row for row in csv.DictReader((DATA / "t08_cases.csv").read_text(encoding="utf-8").splitlines())}
    for case_id in selected:
        observed, _, _ = _evaluate(configurations[case_id])
        expected = existing[case_id]
        for field in ("reference_lmax", "total_converged", "joint_converged"):
            if str(observed[field]) != expected[field]:
                raise RuntimeError(f"audit mismatch for {case_id}: {field}")
        for field in ("epsilon_a", "epsilon_b", "eta", "lambda_max", "rho_l1", "rms_d"):
            if not np.isclose(float(observed[field]), float(expected[field]), rtol=3e-12, atol=3e-13):
                raise RuntimeError(f"audit mismatch for {case_id}: {field}")
    print(f"T08 stratified audit passed for {len(selected)} cases")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-existing", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    if args.audit_existing:
        _audit(); return
    configurations = enumerate_transferability_configurations()
    results = _run(configurations, args.workers)
    cases = [result[0] for result in results]
    forces = [row for result in results for row in result[1]]
    convergence = [row for result in results for row in result[2]]
    _write(DATA / "t08_cases.csv", cases)
    _write(DATA / "t08_forces.csv", forces)
    _write(DATA / "t08_convergence.csv", convergence)
    print(f"T08 generated {len(cases)} cases, {len(forces)} force rows, {len(convergence)} convergence rows")


if __name__ == "__main__":
    main()
