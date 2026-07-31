#!/usr/bin/env python3
"""Run the single official 24-case T14 Model-E campaign."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np

from acoustic_ms import (
    EXPECTED_SCALE_OUT_CASE_IDS,
    build_scale_out_cases,
    canonical_coordinate_hash,
    maximum_geometric_coupling,
    nodal_pair_force_on_probe,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    solve_multipolar_nodal_interaction_forces,
    spectral_radius_l1,
)
from acoustic_ms.external_validation import minimum_two_step_confirmation, successive_change


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
MANIFEST = DATA / "t14_scale_manifest.csv"
PREDICTIONS = DATA / "t14_frozen_predictions.csv"
PROTOCOL = DATA / "t14_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_prior_artifact_hashes.csv"
RAW = DATA / "t14_model_e_convergence.csv"
CACHE = Path("/tmp/silva_bruus_t14_cache")
PHASE_A_PATHS = (MANIFEST, PREDICTIONS, PROTOCOL, PRIOR_HASHES)
RADIUS = 1.0
K = 0.1
ENERGY = 1.0
F0 = 0.0
F1 = 0.8
TOLERANCE = 1.0e-5
DIAGNOSTIC_TOLERANCE = 1.0e-12
CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)
AUDIT_IDS = (
    "t14_n15_linear_level1",
    "t14_n15_compact_level3",
    "t14_n15_irregular_level4",
    "t14_n28_linear_level2",
    "t14_n28_compact_level4",
    "t14_n28_irregular_level1",
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _serialize(values: np.ndarray) -> str:
    array = np.asarray(values, dtype=float)
    return ";".join(
        ":".join(format(float(value), ".17g") for value in row)
        for row in array
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
        raise ValueError("cannot publish an empty T14 campaign")
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


def _verify_prior_hashes() -> None:
    for row in _read(PRIOR_HASHES):
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"prior artifact size changed: {row['path']}")
        if sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise RuntimeError(f"prior artifact hash changed: {row['path']}")


def _manifest() -> list[dict[str, str]]:
    rows = _read(MANIFEST)
    if len(rows) != 24 or tuple(row["case_id"] for row in rows) != EXPECTED_SCALE_OUT_CASE_IDS:
        raise RuntimeError("published T14 manifest identity failed")
    generated = {case.case_id: case for case in build_scale_out_cases(f1=F1)}
    for row in rows:
        case = generated[row["case_id"]]
        stored = _deserialize(row["coordinates_xyz"])
        if not np.array_equal(stored, case.positions_xyz):
            raise RuntimeError(f"coordinates changed for {case.case_id}")
        if canonical_coordinate_hash(stored) != row["coordinate_sha256"]:
            raise RuntimeError(f"coordinate hash failed for {case.case_id}")
        lam = maximum_geometric_coupling(stored, RADIUS, F1)
        if not np.isclose(lam, float(row["lambda_max"]), rtol=5e-13, atol=5e-15):
            raise RuntimeError(f"Lambda_max failed for {case.case_id}")
        rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(
            stored, K, RADIUS, ENERGY, F0, F1, 1
        ))
        if not np.isclose(rho, float(row["rho_l1"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"rho1 failed for {case.case_id}")
    return rows


def _model_a(positions: np.ndarray) -> np.ndarray:
    result = np.zeros((len(positions), 3), dtype=float)
    for first, second in combinations(range(len(positions)), 2):
        result[first, :2] += nodal_pair_force_on_probe(
            positions[first, :2], positions[second, :2], K, RADIUS, ENERGY, F1
        )
        result[second, :2] += nodal_pair_force_on_probe(
            positions[second, :2], positions[first, :2], K, RADIUS, ENERGY, F1
        )
    return result


def _confirmation(results: list[object], attribute: str) -> int:
    changes = [0.0]
    applicable = [False]
    orders = [result.lmax for result in results]
    for index in range(1, len(results)):
        change, flag, _ = successive_change(
            getattr(results[index], attribute), getattr(results[index - 1], attribute)
        )
        changes.append(change)
        applicable.append(flag)
    return minimum_two_step_confirmation(changes, applicable, orders, tolerance=TOLERANCE)


def _all_confirmed(results: list[object]) -> bool:
    return all(_confirmation(results, attribute) > 0 for _, attribute in CHANNELS)


def _evaluate_case(row: dict[str, str]) -> tuple[int, list[dict[str, object]]]:
    positions = _deserialize(row["coordinates_xyz"])
    model_a = _model_a(positions)
    results = []
    timings = []
    accumulated = 0.0
    peak_memory = 0
    for order in range(2, 14):
        started = time.perf_counter()
        result = solve_model_e_nodal(positions, K, RADIUS, ENERGY, F0, F1, order)
        wall = time.perf_counter() - started
        accumulated += wall
        peak_memory = max(peak_memory, int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        results.append(result)
        timings.append((wall, accumulated, peak_memory))
        if order >= 5 and order < 13 and _all_confirmed(results):
            break
    confirmations = {short: _confirmation(results, attribute) for short, attribute in CHANNELS}
    final_lmax = results[-1].lmax
    stop_reason = "all_channels_confirmed" if _all_confirmed(results) else "maximum_lmax_13"
    rows: list[dict[str, object]] = []
    for index, (result, timing) in enumerate(zip(results, timings)):
        solution = result.solution
        previous = results[index - 1] if index else None
        force_arrays = [getattr(result, attribute) for _, attribute in CHANNELS]
        maximum_force_rms = max(rms_vector_magnitude_xyz(values) for values in force_arrays)
        max_abs_fz = max(float(np.max(np.abs(values[:, 2]))) for values in force_arrays)
        finite = bool(
            all(np.all(np.isfinite(values)) for values in force_arrays)
            and np.all(np.isfinite(solution.scattered_coefficients))
            and np.isfinite(solution.balanced_condition_number)
            and np.isfinite(solution.balanced_backward_error)
            and np.isfinite(solution.effective_incident_closure_error)
            and np.isfinite(solution.scattering_closure_error)
            and np.isfinite(result.decomposition_residual)
        )
        expected_dimension = len(positions) * len(solution.active_modes)
        mode_consistent = bool(
            solution.balanced_system_matrix.shape == (expected_dimension, expected_dimension)
            and len(solution.modes) == (result.lmax + 1) ** 2
        )
        fz_tolerance = 128.0 * np.finfo(float).eps * maximum_force_rms
        planar = max_abs_fz <= fz_tolerance
        diagnostics = bool(
            solution.production_solver == "balanced_sqrt"
            and finite
            and solution.balanced_condition_number < 10.0
            and solution.balanced_backward_error < DIAGNOSTIC_TOLERANCE
            and solution.effective_incident_closure_error < DIAGNOSTIC_TOLERANCE
            and solution.scattering_closure_error < DIAGNOSTIC_TOLERANCE
            and result.decomposition_residual < DIAGNOSTIC_TOLERANCE
            and mode_consistent
            and planar
        )
        output: dict[str, object] = {
            "scale_order": int(row["scale_order"]),
            "case_id": row["case_id"],
            "particle_count": int(row["particle_count"]),
            "family": row["family"],
            "target_level": int(row["target_level"]),
            "f1": F1,
            "distance_ratio": float(row["distance_ratio"]),
            "lambda_max": float(row["lambda_max"]),
            "rho_l1": float(row["rho_l1"]),
            "lmax": result.lmax,
            "final_lmax": final_lmax,
            "maximum_allowed_lmax": 13,
            "stop_reason": stop_reason,
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
            "planar_symmetry_pass": planar,
            "finite": finite,
            "diagnostics_pass": diagnostics,
            "order_wall_seconds": timing[0],
            "case_accumulated_seconds": timing[1],
            "process_peak_memory_kib": timing[2],
            "coordinates_xyz": row["coordinates_xyz"],
            "model_a_forces_xyz": _serialize(model_a),
            "total_forces_xyz": _serialize(result.total_forces_xyz),
            "external_forces_xyz": _serialize(result.external_forces_xyz),
            "interaction_forces_xyz": _serialize(result.interaction_forces_xyz),
            "external_scattered_forces_xyz": _serialize(result.external_scattered_forces_xyz),
            "scattered_scattered_forces_xyz": _serialize(result.scattered_scattered_forces_xyz),
            "campaign_complete": True,
        }
        for short, attribute in CHANNELS:
            if previous is None:
                change, applicable, absolute = 0.0, False, 0.0
            else:
                change, applicable, absolute = successive_change(
                    getattr(result, attribute), getattr(previous, attribute)
                )
            output[f"{short}_rms"] = rms_vector_magnitude_xyz(getattr(result, attribute))
            output[f"{short}_successive_change"] = change
            output[f"{short}_absolute_change"] = absolute
            output[f"{short}_change_applicable"] = applicable
            output[f"{short}_minimum_confirmed_lmax"] = confirmations[short]
            output[f"{short}_confirmed"] = confirmations[short] > 0
        rows.append(output)
    return int(row["scale_order"]), rows


def _cache_path(order: int) -> Path:
    return CACHE / f"case_{order:02d}.json"


def _save_cache(order: int, rows: list[dict[str, object]]) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = _cache_path(order)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=CACHE, delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(
            rows, stream, separators=(",", ":"),
            default=lambda value: value.item(),
        )
    temporary.replace(path)


def _load_cache(order: int) -> list[dict[str, object]] | None:
    path = _cache_path(order)
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as stream:
        rows = json.load(stream)
    if not rows or any(int(row["scale_order"]) != order for row in rows):
        raise RuntimeError(f"invalid cache for scale order {order}")
    return rows


def run_campaign(workers: int) -> None:
    if RAW.exists():
        raise RuntimeError("official T14 raw campaign exists; use --audit-existing")
    _verify_phase_a_against_head()
    _verify_prior_hashes()
    manifest = _manifest()
    completed: dict[int, list[dict[str, object]]] = {}
    pending = []
    for row in manifest:
        order = int(row["scale_order"])
        cached = _load_cache(order)
        if cached is None:
            pending.append(row)
        else:
            completed[order] = cached
            print(f"T14 cache restored {order:02d}/24: {row['case_id']}", flush=True)
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1:
        for row in pending:
            order, rows = _evaluate_case(row)
            completed[order] = rows
            _save_cache(order, rows)
            print(f"T14 Model E completed {order:02d}/24: {row['case_id']} L={rows[-1]['final_lmax']}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_evaluate_case, row): row for row in pending}
            for future in as_completed(futures):
                row = futures[future]
                order, rows = future.result()
                completed[order] = rows
                _save_cache(order, rows)
                print(f"T14 Model E completed {order:02d}/24: {row['case_id']} L={rows[-1]['final_lmax']}", flush=True)
    if set(completed) != set(range(1, 25)):
        raise RuntimeError("T14 did not complete all 24 cases")
    rows = [item for order in range(1, 25) for item in completed[order]]
    _atomic_write(RAW, rows)
    print(f"T14 raw campaign published atomically: {len(rows)} order rows")


def audit_existing() -> None:
    _verify_phase_a_against_head()
    _verify_prior_hashes()
    if not RAW.exists():
        raise RuntimeError("cannot audit a missing T14 campaign")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read(RAW):
        grouped.setdefault(row["case_id"], []).append(row)
    manifest = {row["case_id"]: row for row in _manifest()}
    for case_id in AUDIT_IDS:
        expected = sorted(grouped[case_id], key=lambda row: int(row["lmax"]))[-1]
        source = manifest[case_id]
        positions = _deserialize(source["coordinates_xyz"])
        result = solve_model_e_nodal(
            positions, K, RADIUS, ENERGY, F0, F1, int(expected["final_lmax"])
        )
        for short, attribute in CHANNELS:
            if not np.allclose(
                getattr(result, attribute), _deserialize(expected[f"{short}_forces_xyz"]),
                rtol=3e-12, atol=3e-13,
            ):
                raise RuntimeError(f"audit mismatch {case_id} {short}")
        if not np.isclose(
            result.solution.balanced_condition_number,
            float(expected["balanced_condition_number"]), rtol=3e-12, atol=3e-13,
        ):
            raise RuntimeError(f"audit condition mismatch {case_id}")
        lam = maximum_geometric_coupling(positions, RADIUS, F1)
        rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY, F0, F1, 1
        ))
        if not np.isclose(lam, float(source["lambda_max"]), rtol=5e-13, atol=5e-15):
            raise RuntimeError(f"audit Lambda mismatch {case_id}")
        if not np.isclose(rho, float(source["rho_l1"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"audit rho mismatch {case_id}")
        print(f"T14 audit passed: {case_id}", flush=True)
    print(f"T14 stratified audit passed for {len(AUDIT_IDS)} cases")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--audit-existing", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    arguments = parser.parse_args()
    if arguments.audit_existing and arguments.analyze_only:
        parser.error("--audit-existing and --analyze-only are mutually exclusive")
    if arguments.analyze_only:
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "analyze_t14_scale_out.py")],
            cwd=ROOT, check=True,
        )
    elif arguments.audit_existing:
        audit_existing()
    else:
        run_campaign(arguments.workers)


if __name__ == "__main__":
    main()
