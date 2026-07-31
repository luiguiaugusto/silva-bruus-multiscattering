#!/usr/bin/env python3
"""Run or audit the single official T14.1 large-N Model-E campaign."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
from itertools import combinations
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np

from acoustic_ms import (
    EXPECTED_LARGE_N_CASE_IDS,
    build_large_n_cases,
    canonical_coordinate_hash,
    local_coupling_statistics,
    local_geometric_coupling,
    maximum_geometric_coupling,
    nodal_pair_force_on_probe,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    solve_multipolar_nodal,
)
from acoustic_ms.external_validation import minimum_two_step_confirmation, successive_change


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
MANIFEST = DATA / "t14_1_large_n_manifest.csv"
LOCAL = DATA / "t14_1_local_coupling.csv"
PREDICTIONS = DATA / "t14_1_frozen_predictions.csv"
PROTOCOL = DATA / "t14_1_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_1_prior_artifact_hashes.csv"
RAW = DATA / "t14_1_model_e_convergence.csv"
CACHE = Path("/tmp/silva_bruus_t14_1_cache")
PHASE_A = (MANIFEST, LOCAL, PREDICTIONS, PROTOCOL, PRIOR_HASHES)
RADIUS = 1.0
K = 0.1
ENERGY = 1.0
F0 = 0.0
F1 = 0.8
TOLERANCE = 1.0e-5
DIAGNOSTIC_TOLERANCE = 1.0e-12
MEMORY_FACTOR = 12.0
MEMORY_FRACTION = 0.75
CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)
AUDIT_IDS = (
    "t14_1_n45_linear_level1",
    "t14_1_n45_compact_level3",
    "t14_1_n45_irregular_level4",
    "t14_1_n105_linear_level2",
    "t14_1_n105_compact_level4",
    "t14_1_n105_irregular_level1",
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
    return ";".join(":".join(format(float(value), ".17g") for value in row) for row in array)


def _deserialize(value: str) -> np.ndarray:
    array = np.asarray([[float(item) for item in row.split(":")] for row in value.split(";")])
    if array.ndim != 2 or array.shape[1] != 3 or not np.all(np.isfinite(array)):
        raise ValueError("serialized vector field is invalid")
    return array


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot publish an empty T14.1 campaign")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise ValueError("campaign rows have inconsistent fields")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(value) for key, value in row.items()})
    temporary.replace(path)


def _phase_a_hash() -> str:
    digest = sha256()
    for path in PHASE_A:
        digest.update(path.name.encode("ascii"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _verify_phase_a_against_head() -> None:
    for path in PHASE_A:
        relative = str(path.relative_to(ROOT))
        expected = subprocess.run(["git", "show", f"HEAD:{relative}"], cwd=ROOT, check=True, capture_output=True).stdout
        if sha256(expected).digest() != sha256(path.read_bytes()).digest():
            raise RuntimeError(f"published phase-A artifact changed: {relative}")
    protocol = _read(PROTOCOL)
    for row in protocol:
        if row["category"] == "phase_a_code_sha256":
            path = ROOT / row["key"]
            if sha256(path.read_bytes()).hexdigest() != row["value"]:
                raise RuntimeError(f"published phase-A code changed: {row['key']}")


def _verify_prior_hashes() -> None:
    for row in _read(PRIOR_HASHES):
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"prior artifact size changed: {row['path']}")
        if sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise RuntimeError(f"prior artifact hash changed: {row['path']}")


def _rho_l1(positions: np.ndarray) -> float:
    solution = solve_multipolar_nodal(positions, K, RADIUS, F0, F1, 1)
    matrix = solution.system_matrix
    coupling = np.eye(len(matrix), dtype=complex) - matrix
    return 0.0 if matrix.size == 0 else float(np.max(np.abs(np.linalg.eigvals(coupling))))


def _manifest() -> list[dict[str, str]]:
    rows = _read(MANIFEST)
    if len(rows) != 24 or tuple(row["case_id"] for row in rows) != EXPECTED_LARGE_N_CASE_IDS:
        raise RuntimeError("published T14.1 manifest identity failed")
    generated = {case.case_id: case for case in build_large_n_cases(f1=F1)}
    for row in rows:
        case = generated[row["case_id"]]
        stored = _deserialize(row["coordinates_xyz"])
        if not np.array_equal(stored, case.positions_xyz):
            raise RuntimeError(f"coordinates changed for {case.case_id}")
        if canonical_coordinate_hash(stored) != row["coordinate_sha256"]:
            raise RuntimeError(f"coordinate hash failed for {case.case_id}")
        local = local_geometric_coupling(stored, RADIUS, F1)
        stats = local_coupling_statistics(local)
        if stats["local_coupling_sha256"] != row["local_coupling_sha256"]:
            raise RuntimeError(f"local coupling changed for {case.case_id}")
        lam = maximum_geometric_coupling(stored, RADIUS, F1)
        if not np.isclose(lam, float(row["lambda_max"]), rtol=5e-13, atol=5e-15):
            raise RuntimeError(f"Lambda_max failed for {case.case_id}")
        rho = _rho_l1(stored)
        if not np.isclose(rho, float(row["rho_l1"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"rho1 failed for {case.case_id}")
    return rows


def _model_a(positions: np.ndarray) -> np.ndarray:
    result = np.zeros((len(positions), 3), dtype=float)
    for first, second in combinations(range(len(positions)), 2):
        result[first, :2] += nodal_pair_force_on_probe(positions[first, :2], positions[second, :2], K, RADIUS, ENERGY, F1)
        result[second, :2] += nodal_pair_force_on_probe(positions[second, :2], positions[first, :2], K, RADIUS, ENERGY, F1)
    return result


def _confirmation(rows: list[dict[str, object]], channel: str) -> int:
    orders = [int(row["lmax"]) for row in rows]
    changes = [float(row[f"{channel}_successive_change"]) for row in rows]
    applicable = [bool(row[f"{channel}_change_applicable"]) for row in rows]
    return minimum_two_step_confirmation(changes, applicable, orders, tolerance=TOLERANCE)


def _all_confirmed(rows: list[dict[str, object]]) -> bool:
    return len(rows) >= 3 and all(_confirmation(rows, channel) > 0 for channel, _ in CHANNELS)


def _available_memory() -> int:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return 0
    for line in meminfo.read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    return 0


def _resource_estimate(particle_count: int, lmax: int) -> tuple[int, int, int, bool]:
    active = lmax * (lmax + 1) // 2
    dimension = particle_count * active
    estimated = int(MEMORY_FACTOR * 16 * dimension * dimension)
    available = _available_memory()
    allowed = available == 0 or estimated <= int(MEMORY_FRACTION * available)
    return active, dimension, estimated, allowed


def _signature(row: dict[str, str]) -> str:
    payload = "|".join((row["case_id"], row["coordinate_sha256"], _phase_a_hash(), sys.version, np.__version__))
    return sha256(payload.encode("utf-8")).hexdigest()


def _cache_path(order: int) -> Path:
    return CACHE / f"case_{order:02d}.json"


def _save_cache(order: int, signature: str, rows: list[dict[str, object]], resource_failed: bool = False) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = _cache_path(order)
    document = {"signature": signature, "resource_failed": resource_failed, "rows": rows}
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=CACHE, delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(document, stream, separators=(",", ":"), default=lambda value: value.item())
    temporary.replace(path)


def _load_cache(order: int, signature: str, resume: bool) -> tuple[list[dict[str, object]], bool]:
    path = _cache_path(order)
    if not path.exists():
        return [], False
    if not resume:
        raise RuntimeError(f"cache exists for case {order}; rerun with --resume")
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("signature") != signature:
        raise RuntimeError(f"cache signature mismatch for case {order}")
    rows = document.get("rows", [])
    if any(int(row["scale_order"]) != order for row in rows):
        raise RuntimeError(f"invalid cache rows for case {order}")
    return rows, bool(document.get("resource_failed", False))


def _order_row(
    source: dict[str, str], model_a: np.ndarray, result: object,
    previous: dict[str, object] | None, wall: float, accumulated: float,
    peak_memory: int, active_expected: int, dimension_expected: int,
    estimated_bytes: int, available_bytes: int, blas_threads: int,
) -> dict[str, object]:
    solution = result.solution
    arrays = [getattr(result, attribute) for _, attribute in CHANNELS]
    maximum_force_rms = max(rms_vector_magnitude_xyz(values) for values in arrays)
    max_abs_fz = max(float(np.max(np.abs(values[:, 2]))) for values in arrays)
    finite = bool(
        all(np.all(np.isfinite(values)) for values in arrays)
        and np.all(np.isfinite(solution.scattered_coefficients))
        and all(np.isfinite(value) for value in (
            solution.balanced_condition_number, solution.balanced_backward_error,
            solution.effective_incident_closure_error, solution.scattering_closure_error,
            result.decomposition_residual,
        ))
    )
    observed_dimension = solution.balanced_system_matrix.shape[0]
    mode_consistent = bool(
        len(solution.active_modes) == active_expected
        and observed_dimension == dimension_expected
        and len(solution.modes) == (result.lmax + 1) ** 2
    )
    fz_tolerance = 128.0 * np.finfo(float).eps * maximum_force_rms
    planar = max_abs_fz <= fz_tolerance
    diagnostics = bool(
        solution.production_solver == "balanced_sqrt" and finite
        and solution.balanced_condition_number < 10.0
        and solution.balanced_backward_error < DIAGNOSTIC_TOLERANCE
        and solution.effective_incident_closure_error < DIAGNOSTIC_TOLERANCE
        and solution.scattering_closure_error < DIAGNOSTIC_TOLERANCE
        and result.decomposition_residual < DIAGNOSTIC_TOLERANCE
        and mode_consistent and planar
    )
    output: dict[str, object] = {
        "scale_order": int(source["scale_order"]), "case_id": source["case_id"],
        "particle_count": int(source["particle_count"]), "family": source["family"],
        "target_level": int(source["target_level"]), "f1": F1,
        "distance_ratio": float(source["distance_ratio"]),
        "lambda_max": float(source["lambda_max"]), "rho_l1": float(source["rho_l1"]),
        "lmax": result.lmax, "final_lmax": 0, "maximum_allowed_lmax": 13,
        "stop_reason": "in_progress", "full_modes_per_particle": len(solution.modes),
        "active_modes_per_particle": len(solution.active_modes),
        "predicted_active_modes_per_particle": active_expected,
        "system_dimension": observed_dimension, "predicted_system_dimension": dimension_expected,
        "estimated_memory_bytes": estimated_bytes, "available_memory_bytes": available_bytes,
        "memory_safety_fraction": MEMORY_FRACTION, "mode_dimension_consistent": mode_consistent,
        "production_solver": solution.production_solver,
        "physical_residual_relative": solution.residual_relative,
        "balanced_condition_number": solution.balanced_condition_number,
        "balanced_backward_error": solution.balanced_backward_error,
        "effective_incident_closure_error": solution.effective_incident_closure_error,
        "scattering_closure_error": solution.scattering_closure_error,
        "force_decomposition_residual": result.decomposition_residual,
        "max_abs_fz": max_abs_fz, "fz_tolerance": fz_tolerance,
        "planar_symmetry_pass": planar, "finite": finite, "diagnostics_pass": diagnostics,
        "assembly_seconds": 0.0, "solve_seconds": 0.0, "postprocess_seconds": 0.0,
        "timing_breakdown_available": False, "order_wall_seconds": wall,
        "case_accumulated_seconds": accumulated, "process_peak_memory_kib": peak_memory,
        "workers": 1, "blas_threads": blas_threads, "resource_precheck_pass": True,
        "coordinates_xyz": source["coordinates_xyz"], "coordinate_sha256": source["coordinate_sha256"],
        "local_coupling_sha256": source["local_coupling_sha256"],
        "model_a_forces_xyz": _serialize(model_a),
        "total_forces_xyz": _serialize(result.total_forces_xyz),
        "external_forces_xyz": _serialize(result.external_forces_xyz),
        "interaction_forces_xyz": _serialize(result.interaction_forces_xyz),
        "external_scattered_forces_xyz": _serialize(result.external_scattered_forces_xyz),
        "scattered_scattered_forces_xyz": _serialize(result.scattered_scattered_forces_xyz),
        "campaign_complete": False,
    }
    for channel, attribute in CHANNELS:
        values = getattr(result, attribute)
        if previous is None:
            change, applicable, absolute = 0.0, False, 0.0
        else:
            change, applicable, absolute = successive_change(values, _deserialize(str(previous[f"{channel}_forces_xyz"])))
        output[f"{channel}_rms"] = rms_vector_magnitude_xyz(values)
        output[f"{channel}_successive_change"] = change
        output[f"{channel}_absolute_change"] = absolute
        output[f"{channel}_change_applicable"] = applicable
        output[f"{channel}_minimum_confirmed_lmax"] = 0
        output[f"{channel}_confirmed"] = False
    return output


def _finalize(rows: list[dict[str, object]], resource_failed: bool) -> None:
    confirmations = {channel: _confirmation(rows, channel) for channel, _ in CHANNELS}
    final_lmax = int(rows[-1]["lmax"])
    all_confirmed = all(value > 0 for value in confirmations.values())
    stop = "all_channels_confirmed" if all_confirmed else ("resource_precheck_failed" if resource_failed else "unconfirmed_at_13")
    for row in rows:
        row["final_lmax"] = final_lmax
        row["stop_reason"] = stop
        row["campaign_complete"] = True
        for channel, _ in CHANNELS:
            row[f"{channel}_minimum_confirmed_lmax"] = confirmations[channel]
            row[f"{channel}_confirmed"] = confirmations[channel] > 0


def _evaluate_case(source: dict[str, str], resume: bool, blas_threads: int) -> tuple[list[dict[str, object]], bool]:
    signature = _signature(source)
    rows, resource_failed = _load_cache(int(source["scale_order"]), signature, resume)
    if resource_failed:
        return rows, True
    positions = _deserialize(source["coordinates_xyz"])
    model_a = _model_a(positions)
    accumulated = sum(float(row["order_wall_seconds"]) for row in rows)
    start_order = int(rows[-1]["lmax"]) + 1 if rows else 2
    if rows and int(rows[-1]["lmax"]) >= 5 and _all_confirmed(rows):
        _finalize(rows, False)
        return rows, False
    for order in range(start_order, 14):
        active, dimension, estimated, allowed = _resource_estimate(len(positions), order)
        available = _available_memory()
        print(f"T14.1 precheck {source['case_id']} L={order}: dim={dimension} estimate={estimated} available={available} pass={allowed}", flush=True)
        if not allowed:
            resource_failed = True
            _save_cache(int(source["scale_order"]), signature, rows, True)
            break
        started = time.perf_counter()
        result = solve_model_e_nodal(positions, K, RADIUS, ENERGY, F0, F1, order)
        wall = time.perf_counter() - started
        accumulated += wall
        previous = rows[-1] if rows else None
        rows.append(_order_row(
            source, model_a, result, previous, wall, accumulated,
            int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss), active,
            dimension, estimated, available, blas_threads,
        ))
        _save_cache(int(source["scale_order"]), signature, rows)
        if order >= 5 and order < 13 and _all_confirmed(rows):
            break
    if not rows:
        raise RuntimeError(f"resource precheck prevented every order for {source['case_id']}")
    _finalize(rows, resource_failed)
    _save_cache(int(source["scale_order"]), signature, rows, resource_failed)
    return rows, resource_failed


def _thread_environment(blas_threads: int) -> None:
    expected = str(blas_threads)
    keys = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    if any(os.environ.get(key) != expected for key in keys):
        environment = os.environ.copy()
        for key in keys:
            environment[key] = expected
        os.execve(sys.executable, [sys.executable, *sys.argv], environment)


def run_campaign(*, resume: bool, workers: int, blas_threads: int) -> None:
    if workers != 1:
        raise ValueError("the official resource-safe T14.1 campaign requires --workers 1")
    _thread_environment(blas_threads)
    if RAW.exists():
        raise RuntimeError("official T14.1 raw campaign exists; use --audit-existing")
    _verify_phase_a_against_head()
    _verify_prior_hashes()
    manifest = _manifest()
    completed: dict[int, list[dict[str, object]]] = {}
    resource_limit = False
    for source in manifest:
        order = int(source["scale_order"])
        rows, failed = _evaluate_case(source, resume, blas_threads)
        completed[order] = rows
        resource_limit |= failed
        print(f"T14.1 completed {order:02d}/24: {source['case_id']} L={rows[-1]['final_lmax']} stop={rows[-1]['stop_reason']}", flush=True)
    official = [item for order in range(1, 25) for item in completed[order]]
    _atomic_write(RAW, official)
    print(f"T14.1 raw campaign published atomically: {len(official)} rows; resource_limit={resource_limit}")


def audit_existing(blas_threads: int) -> None:
    _thread_environment(blas_threads)
    _verify_phase_a_against_head()
    _verify_prior_hashes()
    if not RAW.exists():
        raise RuntimeError("cannot audit a missing T14.1 campaign")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read(RAW):
        grouped.setdefault(row["case_id"], []).append(row)
    manifest = {row["case_id"]: row for row in _manifest()}
    for case_id in AUDIT_IDS:
        expected = sorted(grouped[case_id], key=lambda row: int(row["lmax"]))[-1]
        source = manifest[case_id]
        positions = _deserialize(source["coordinates_xyz"])
        result = solve_model_e_nodal(positions, K, RADIUS, ENERGY, F0, F1, int(expected["final_lmax"]))
        for channel, attribute in CHANNELS:
            if not np.allclose(getattr(result, attribute), _deserialize(expected[f"{channel}_forces_xyz"]), rtol=3e-12, atol=3e-13):
                raise RuntimeError(f"audit mismatch {case_id} {channel}")
        local = local_geometric_coupling(positions, RADIUS, F1)
        if local_coupling_statistics(local)["local_coupling_sha256"] != source["local_coupling_sha256"]:
            raise RuntimeError(f"audit local coupling mismatch {case_id}")
        lam = maximum_geometric_coupling(positions, RADIUS, F1)
        rho = _rho_l1(positions)
        if not np.isclose(lam, float(source["lambda_max"]), rtol=5e-13, atol=5e-15):
            raise RuntimeError(f"audit Lambda mismatch {case_id}")
        if not np.isclose(rho, float(source["rho_l1"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"audit rho mismatch {case_id}")
        if not np.isclose(result.solution.balanced_condition_number, float(expected["balanced_condition_number"]), rtol=3e-12, atol=3e-13):
            raise RuntimeError(f"audit condition mismatch {case_id}")
        print(f"T14.1 audit passed: {case_id}", flush=True)
    print(f"T14.1 stratified audit passed for {len(AUDIT_IDS)} cases")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--blas-threads", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--audit-existing", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    arguments = parser.parse_args()
    if sum((arguments.audit_existing, arguments.analyze_only)) > 1:
        parser.error("--audit-existing and --analyze-only are mutually exclusive")
    if arguments.blas_threads < 1:
        parser.error("--blas-threads must be positive")
    if arguments.analyze_only:
        subprocess.run([sys.executable, str(ROOT / "scripts" / "analyze_t14_1_large_n.py")], cwd=ROOT, check=True)
    elif arguments.audit_existing:
        audit_existing(arguments.blas_threads)
    else:
        run_campaign(resume=arguments.resume, workers=arguments.workers, blas_threads=arguments.blas_threads)


if __name__ == "__main__":
    main()
