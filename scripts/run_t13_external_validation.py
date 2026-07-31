#!/usr/bin/env python3
"""Run exactly the published 24-case T13 Model-E holdout campaign."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from acoustic_ms import (
    EXPECTED_CASE_IDS,
    cluster_family,
    maximum_geometric_coupling,
    nodal_pair_force_on_probe,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    solve_multipolar_nodal_interaction_forces,
    spectral_radius_l1,
)
from acoustic_ms.external_validation import (
    canonical_coordinate_hash,
    minimum_two_step_confirmation,
    successive_change,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
MANIFEST = DATA / "t13_holdout_manifest.csv"
PREDICTIONS = DATA / "t13_frozen_predictions.csv"
PROTOCOL = DATA / "t13_frozen_protocol.csv"
RAW_PATH = DATA / "t13_model_e_convergence.csv"
T08_CASES = DATA / "t08_cases.csv"
T08_FORCES = DATA / "t08_forces.csv"
CACHE = Path("/tmp/silva_bruus_t13_cache")
RADIUS = 1.0
K = 0.1
KA = 0.1
ENERGY = 1.0
F0 = 0.0
TOLERANCE = 1.0e-5
DIAGNOSTIC_TOLERANCE = 1.0e-12
CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)
PHASE_A_PATHS = (MANIFEST, PREDICTIONS, PROTOCOL)
AUDIT_ORDERS = (1, 6, 11, 8, 13, 18, 23, 24)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _serialize(vectors: np.ndarray) -> str:
    values = np.asarray(vectors, dtype=float)
    return ";".join(
        ":".join(format(float(value), ".17g") for value in row)
        for row in values
    )


def _deserialize(value: str) -> np.ndarray:
    result = np.asarray([
        [float(item) for item in row.split(":")] for row in value.split(";")
    ])
    if result.ndim != 2 or result.shape[1] != 3 or not np.all(np.isfinite(result)):
        raise ValueError("serialized vector field is invalid")
    return result


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot publish an empty T13 campaign")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise ValueError("campaign rows have inconsistent fields")
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(value) for key, value in row.items()})
    temporary.replace(path)


def _verify_phase_a_against_head() -> None:
    for path in PHASE_A_PATHS:
        relative = str(path.relative_to(ROOT))
        expected = subprocess.run(
            ["git", "show", f"HEAD:{relative}"], cwd=ROOT, check=True,
            capture_output=True,
        ).stdout
        if sha256(expected).digest() != sha256(path.read_bytes()).digest():
            raise RuntimeError(f"published phase-A artifact changed: {relative}")


def _model_a(positions: np.ndarray, f1: float) -> np.ndarray:
    forces = np.zeros((len(positions), 2), dtype=float)
    for first, second in combinations(range(len(positions)), 2):
        forces[first] += nodal_pair_force_on_probe(
            positions[first, :2], positions[second, :2], K, RADIUS, ENERGY, f1
        )
        forces[second] += nodal_pair_force_on_probe(
            positions[second, :2], positions[first, :2], K, RADIUS, ENERGY, f1
        )
    return forces


def _frozen_t08() -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    cases = {row["case_id"]: row for row in _read(T08_CASES)}
    forces: dict[str, list[dict[str, str]]] = {}
    for row in _read(T08_FORCES):
        forces.setdefault(row["case_id"], []).append(row)
    for rows in forces.values():
        rows.sort(key=lambda row: int(row["particle_index"]))
    return cases, forces


def _positions(row: dict[str, str]) -> np.ndarray:
    return cluster_family(
        int(row["particle_count"]), row["family"], float(row["distance_ratio"])
    )


def audit_manifest_and_models() -> None:
    """Audit frozen geometry, predictors, and T08 A/D vectors before E."""

    manifest = _read(MANIFEST)
    if len(manifest) != 24 or tuple(row["case_id"] for row in manifest) != EXPECTED_CASE_IDS:
        raise RuntimeError("published T13 manifest changed")
    cases, forces = _frozen_t08()
    for row in manifest:
        case_id = row["case_id"]
        positions = _positions(row)
        if canonical_coordinate_hash(positions) != row["coordinate_sha256"]:
            raise RuntimeError(f"coordinate hash failed for {case_id}")
        calculated_lambda = maximum_geometric_coupling(
            positions, RADIUS, float(row["f1"])
        )
        if not np.isclose(calculated_lambda, float(row["lambda_max"]), rtol=3e-14, atol=3e-15):
            raise RuntimeError(f"Lambda_max audit failed for {case_id}")
        l1 = solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY, F0, float(row["f1"]), 1
        )
        calculated_rho = spectral_radius_l1(l1)
        if not np.isclose(calculated_rho, float(row["rho_l1"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"rho1 audit failed for {case_id}")
        frozen_rows = forces[case_id]
        frozen_positions = np.asarray([
            [float(item[field]) for field in ("x", "y", "z")] for item in frozen_rows
        ])
        if not np.allclose(frozen_positions, positions, rtol=0.0, atol=5e-15):
            raise RuntimeError(f"T08 positions changed for {case_id}")
        frozen_a = np.asarray([
            [float(item["a_x"]), float(item["a_y"])] for item in frozen_rows
        ])
        frozen_d = np.asarray([
            [float(item["d_x"]), float(item["d_y"])] for item in frozen_rows
        ])
        recalculated_a = _model_a(positions, float(row["f1"]))
        recalculated_d = solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY, F0, float(row["f1"]),
            int(row["reference_lmax_d"]),
        ).forces_xy
        if not np.allclose(recalculated_a, frozen_a, rtol=5e-12, atol=5e-14):
            raise RuntimeError(f"Model A audit failed for {case_id}")
        if not np.allclose(recalculated_d, frozen_d, rtol=5e-12, atol=5e-14):
            raise RuntimeError(f"Model D audit failed for {case_id}")
        if cases[case_id]["total_converged"] != "true":
            raise RuntimeError(f"T08 D reference is not confirmed for {case_id}")


def _confirmation(results: list[object], attribute: str) -> int:
    orders = [result.lmax for result in results]
    changes: list[float] = []
    applicable: list[bool] = []
    for index, result in enumerate(results):
        if index == 0:
            changes.append(0.0)
            applicable.append(False)
        else:
            change, flag, _ = successive_change(
                getattr(result, attribute), getattr(results[index - 1], attribute)
            )
            changes.append(change)
            applicable.append(flag)
    return minimum_two_step_confirmation(
        changes, applicable, orders, tolerance=TOLERANCE
    )


def _all_confirmed(results: list[object]) -> bool:
    return all(_confirmation(results, attribute) > 0 for _, attribute in CHANNELS)


def _interaction_confirmed(results: list[object]) -> bool:
    return _confirmation(results, "interaction_forces_xyz") > 0


def _evaluate_case(row: dict[str, str]) -> tuple[int, list[dict[str, object]]]:
    positions = _positions(row)
    results = []
    for order in range(2, 22):
        result = solve_model_e_nodal(
            positions, K, RADIUS, ENERGY, F0, float(row["f1"]), order
        )
        results.append(result)
        if order >= 5 and order < 13 and _all_confirmed(results):
            break
        if order == 13 and _interaction_confirmed(results):
            break
        if order > 13 and _interaction_confirmed(results):
            break
    confirmations = {
        short: _confirmation(results, attribute) for short, attribute in CHANNELS
    }
    final_lmax = results[-1].lmax
    campaign_rows: list[dict[str, object]] = []
    for index, result in enumerate(results):
        solution = result.solution
        previous = results[index - 1] if index else None
        force_arrays = [getattr(result, attribute) for _, attribute in CHANNELS]
        maximum_force_rms = max(rms_vector_magnitude_xyz(values) for values in force_arrays)
        max_abs_fz = max(float(np.max(np.abs(values[:, 2]))) for values in force_arrays)
        finite = bool(
            all(np.all(np.isfinite(values)) for values in force_arrays)
            and np.isfinite(solution.balanced_condition_number)
            and np.isfinite(solution.balanced_backward_error)
            and np.isfinite(solution.effective_incident_closure_error)
            and np.isfinite(solution.scattering_closure_error)
            and np.isfinite(result.decomposition_residual)
        )
        expected_dimension = len(positions) * len(solution.active_modes)
        mode_consistent = (
            solution.balanced_system_matrix.shape == (expected_dimension, expected_dimension)
            and len(solution.modes) == (result.lmax + 1) ** 2
        )
        fz_tolerance = 128.0 * np.finfo(float).eps * maximum_force_rms
        diagnostics_pass = bool(
            solution.production_solver == "balanced_sqrt"
            and finite
            and solution.balanced_condition_number < 10.0
            and solution.balanced_backward_error < DIAGNOSTIC_TOLERANCE
            and solution.effective_incident_closure_error < DIAGNOSTIC_TOLERANCE
            and solution.scattering_closure_error < DIAGNOSTIC_TOLERANCE
            and result.decomposition_residual < DIAGNOSTIC_TOLERANCE
            and max_abs_fz <= fz_tolerance
            and mode_consistent
        )
        output: dict[str, object] = {
            "holdout_order": int(row["holdout_order"]),
            "case_id": row["case_id"],
            "particle_count": int(row["particle_count"]),
            "family": row["family"],
            "stratum": row["stratum"],
            "target_level": int(row["target_level"]),
            "f1": float(row["f1"]),
            "distance_ratio": float(row["distance_ratio"]),
            "lambda_max": float(row["lambda_max"]),
            "rho_l1": float(row["rho_l1"]),
            "lmax": result.lmax,
            "final_lmax": final_lmax,
            "standard_lmax_cap": 13,
            "extension_lmax_cap": 21,
            "extended_beyond_13": final_lmax > 13,
            "full_modes_per_particle": len(solution.modes),
            "active_modes_per_particle": len(solution.active_modes),
            "system_dimension": solution.balanced_system_matrix.shape[0],
            "expected_system_dimension": expected_dimension,
            "mode_dimension_consistent": mode_consistent,
            "production_solver": solution.production_solver,
            "physical_residual_relative": solution.residual_relative,
            "balanced_condition_number": solution.balanced_condition_number,
            "balanced_backward_error": solution.balanced_backward_error,
            "effective_incident_closure_error": solution.effective_incident_closure_error,
            "scattering_closure_error": solution.scattering_closure_error,
            "force_decomposition_residual": result.decomposition_residual,
            "max_abs_fz": max_abs_fz,
            "fz_tolerance": fz_tolerance,
            "finite": finite,
            "diagnostics_pass": diagnostics_pass,
            "coordinates_xyz": _serialize(positions),
            "total_forces_xyz": _serialize(result.total_forces_xyz),
            "external_forces_xyz": _serialize(result.external_forces_xyz),
            "interaction_forces_xyz": _serialize(result.interaction_forces_xyz),
            "external_scattered_forces_xyz": _serialize(result.external_scattered_forces_xyz),
            "scattered_scattered_forces_xyz": _serialize(result.scattered_scattered_forces_xyz),
            "campaign_complete": True,
        }
        for short, attribute in CHANNELS:
            current = getattr(result, attribute)
            if previous is None:
                change, applicable, absolute = 0.0, False, 0.0
            else:
                change, applicable, absolute = successive_change(
                    current, getattr(previous, attribute)
                )
            output[f"{short}_rms"] = rms_vector_magnitude_xyz(current)
            output[f"{short}_successive_change"] = change
            output[f"{short}_absolute_change"] = absolute
            output[f"{short}_change_applicable"] = applicable
            output[f"{short}_minimum_confirmed_lmax"] = confirmations[short]
            output[f"{short}_confirmed"] = confirmations[short] > 0
        campaign_rows.append(output)
    return int(row["holdout_order"]), campaign_rows


def _cache_path(order: int) -> Path:
    return CACHE / f"case_{order:02d}.json"


def _save_cache(order: int, rows: list[dict[str, object]]) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = _cache_path(order)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=CACHE, delete=False
    ) as stream:
        temporary = Path(stream.name)
        json.dump(rows, stream, sort_keys=False, separators=(",", ":"))
    temporary.replace(path)


def _load_cache(order: int) -> list[dict[str, object]] | None:
    path = _cache_path(order)
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as stream:
        rows = json.load(stream)
    if not rows or any(int(row["holdout_order"]) != order for row in rows):
        raise RuntimeError(f"invalid temporary cache for holdout order {order}")
    return rows


def run_campaign(workers: int) -> None:
    if RAW_PATH.exists():
        raise RuntimeError("official T13 raw campaign already exists; use --audit-existing")
    _verify_phase_a_against_head()
    audit_manifest_and_models()
    manifest = _read(MANIFEST)
    completed: dict[int, list[dict[str, object]]] = {}
    pending = []
    for row in manifest:
        order = int(row["holdout_order"])
        cached = _load_cache(order)
        if cached is None:
            pending.append(row)
        else:
            completed[order] = cached
            print(f"T13 cache restored {order:02d}/24: {row['case_id']}", flush=True)
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1:
        for row in pending:
            order, rows = _evaluate_case(row)
            completed[order] = rows
            _save_cache(order, rows)
            print(f"T13 Model E completed {order:02d}/24: {row['case_id']} L={rows[-1]['final_lmax']}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_evaluate_case, row): row for row in pending}
            for future in as_completed(futures):
                row = futures[future]
                order, rows = future.result()
                completed[order] = rows
                _save_cache(order, rows)
                print(f"T13 Model E completed {order:02d}/24: {row['case_id']} L={rows[-1]['final_lmax']}", flush=True)
    if set(completed) != set(range(1, 25)):
        raise RuntimeError("T13 did not complete all 24 frozen cases")
    rows = [item for order in range(1, 25) for item in completed[order]]
    _atomic_write(RAW_PATH, rows)
    print(f"T13 raw campaign published atomically: {len(rows)} order rows")


def audit_existing() -> None:
    _verify_phase_a_against_head()
    if not RAW_PATH.exists():
        raise RuntimeError("cannot audit a missing T13 raw campaign")
    raw = _read(RAW_PATH)
    grouped: dict[int, list[dict[str, str]]] = {}
    for row in raw:
        grouped.setdefault(int(row["holdout_order"]), []).append(row)
    manifest = {int(row["holdout_order"]): row for row in _read(MANIFEST)}
    for order in AUDIT_ORDERS:
        expected = grouped[order][-1]
        source = manifest[order]
        result = solve_model_e_nodal(
            _positions(source), K, RADIUS, ENERGY, F0, float(source["f1"]),
            int(expected["final_lmax"]),
        )
        for field, attribute in CHANNELS:
            observed = getattr(result, attribute)
            frozen = _deserialize(expected[f"{field}_forces_xyz"])
            if not np.allclose(observed, frozen, rtol=3e-12, atol=3e-13):
                raise RuntimeError(f"audit mismatch {source['case_id']} {field}")
        if not np.isclose(
            result.solution.balanced_condition_number,
            float(expected["balanced_condition_number"]), rtol=3e-12, atol=3e-13,
        ):
            raise RuntimeError(f"audit condition mismatch {source['case_id']}")
        print(f"T13 audit passed {order:02d}/24: {source['case_id']}", flush=True)
    print(f"T13 stratified audit passed for {len(AUDIT_ORDERS)} cases")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--audit-existing", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    arguments = parser.parse_args()
    if arguments.audit_existing and arguments.analyze_only:
        parser.error("--audit-existing and --analyze-only are mutually exclusive")
    if arguments.analyze_only:
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "analyze_t13_external_validation.py")],
            cwd=ROOT,
            check=True,
        )
    elif arguments.audit_existing:
        audit_existing()
    else:
        run_campaign(arguments.workers)


if __name__ == "__main__":
    main()
