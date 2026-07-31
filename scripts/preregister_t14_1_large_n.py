#!/usr/bin/env python3
"""Publish the response-blind T14.1 large-N preregistration."""

from __future__ import annotations

import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import tempfile

import matplotlib
import numpy as np
import scipy

from acoustic_ms import (
    EXPECTED_LARGE_N_CASE_IDS,
    build_large_n_cases,
    canonical_coordinate_hash,
    frozen_external_predictions,
    local_coupling_statistics,
    maximum_geometric_coupling,
    solve_multipolar_nodal,
)
from acoustic_ms.external_validation import (
    LAMBDA_TARGETS,
    LAMBDA_THRESHOLDS,
    M1_EXPONENT,
    M1_PREFACTOR,
    M1_SAFETY_FACTOR,
    P3_EXPONENT,
    P3_PREFACTOR,
    P3_SAFETY_FACTOR,
    RHO_THRESHOLDS,
    TOLERANCES,
)
from acoustic_ms.scale_out_validation import IRREGULAR_AMPLITUDE


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
BASE_COMMIT = "7c6329ffa6208ac517602f4e411172128fc43465"
MANIFEST = DATA / "t14_1_large_n_manifest.csv"
LOCAL = DATA / "t14_1_local_coupling.csv"
PREDICTIONS = DATA / "t14_1_frozen_predictions.csv"
PROTOCOL = DATA / "t14_1_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_1_prior_artifact_hashes.csv"
PHASE_B = (
    DATA / "t14_1_model_e_convergence.csv",
    DATA / "t14_1_forces.csv",
    DATA / "t14_1_case_summary.csv",
    DATA / "t14_1_large_n_predictions.csv",
    DATA / "t14_1_metrics.csv",
    DATA / "t14_1_threshold_audit.csv",
    DATA / "t14_1_matched_large_n_pairs.csv",
    DATA / "t14_1_combined_scale_sequence.csv",
    DATA / "t14_1_performance.csv",
    DATA / "t14_1_gate.csv",
    FIGURES / "t14_1_large_n_validation.png",
)
CODE_PATHS = (
    "src/acoustic_ms/large_n_validation.py",
    "scripts/preregister_t14_1_large_n.py",
    "scripts/run_t14_1_large_n.py",
    "scripts/analyze_t14_1_large_n.py",
    "tests/test_t14_1_preregistration.py",
)
RADIUS = 1.0
K = 0.1
KA = 0.1
ENERGY = 1.0
F0 = 0.0
F1 = 0.8


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot publish empty table {path.name}")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise ValueError("rows must share exact field order")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(value) for key, value in row.items()})
    temporary.replace(path)


def _serialize(array: np.ndarray) -> str:
    return ";".join(":".join(format(float(value), ".17g") for value in row) for row in array)


def _git_bytes(commit: str, relative: str) -> bytes:
    return subprocess.run(["git", "show", f"{commit}:{relative}"], cwd=ROOT, check=True, capture_output=True).stdout


def _prior_rows() -> list[dict[str, object]]:
    paths = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", BASE_COMMIT, "results/data", "results/figures"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    rows = []
    for relative in paths:
        payload = _git_bytes(BASE_COMMIT, relative)
        rows.append({"path": relative, "size_bytes": len(payload), "sha256": sha256(payload).hexdigest(), "source_commit": BASE_COMMIT})
    if len(rows) != 95:
        raise RuntimeError(f"expected 95 prior artifacts, found {len(rows)}")
    return rows


def _hardware() -> dict[str, str]:
    memory = "unknown"
    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        memory = next((line.split(":", 1)[1].strip() for line in meminfo.read_text().splitlines() if line.startswith("MemTotal:")), "unknown")
    cpu = platform.processor() or platform.machine()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next((line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines() if line.startswith("model name")), cpu)
    return {"cpu": cpu, "logical_cpus": str(os.cpu_count()), "memory_total": memory, "rlimit_as": str(resource.getrlimit(resource.RLIMIT_AS))}


def _rho_l1(positions: np.ndarray) -> float:
    solution = solve_multipolar_nodal(positions, K, RADIUS, F0, F1, 1)
    matrix = solution.system_matrix
    coupling = np.eye(len(matrix), dtype=complex) - matrix
    return 0.0 if matrix.size == 0 else float(np.max(np.abs(np.linalg.eigvals(coupling))))


def preregister() -> None:
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    if head != BASE_COMMIT:
        raise RuntimeError("T14.1 preregistration must start at the frozen T14 commit")
    existing = [str(path.relative_to(ROOT)) for path in PHASE_B if path.exists()]
    if existing:
        raise RuntimeError(f"response-bearing T14.1 files exist before preregistration: {existing}")
    cases = build_large_n_cases(f1=F1)
    if tuple(case.case_id for case in cases) != EXPECTED_LARGE_N_CASE_IDS:
        raise RuntimeError("the exact T14.1 case checksum failed")
    manifests: list[dict[str, object]] = []
    locals_: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    for case in cases:
        if case.minimum_distance <= 2.0 * RADIUS:
            raise RuntimeError(f"overlap margin failed for {case.case_id}")
        audited_lambda = maximum_geometric_coupling(case.positions_xyz, RADIUS, F1)
        if not np.isclose(audited_lambda, case.lambda_target, rtol=5e-13, atol=5e-15):
            raise RuntimeError(f"Lambda_max audit failed for {case.case_id}")
        stats = local_coupling_statistics(case.local_coupling)
        rho_l1 = _rho_l1(case.positions_xyz)
        manifests.append({
            "scale_order": case.scale_order, "case_id": case.case_id,
            "particle_count": case.particle_count, "family": case.family,
            "target_level": case.target_level, "triangular_rows": case.triangular_rows,
            "radius": RADIUS, "k": K, "ka": KA, "energy_density": ENERGY,
            "f0": F0, "f1": F1, "geometry_factor": case.geometry_factor,
            "lambda_target": case.lambda_target, "distance_ratio": case.distance_ratio,
            "minimum_distance": case.minimum_distance,
            "overlap_margin": case.minimum_distance - 2.0 * RADIUS,
            "cluster_diameter": case.cluster_diameter,
            "k_cluster_diameter": K * case.cluster_diameter,
            "coordinates_xyz": _serialize(case.positions_xyz),
            "coordinate_sha256": canonical_coordinate_hash(case.positions_xyz),
            **stats, "rho_l1": rho_l1,
            "provenance": "pre_registered_large_n_uniform_dilation",
            "model_e_response_consulted": False,
        })
        for index, value in enumerate(case.local_coupling):
            locals_.append({
                "scale_order": case.scale_order, "case_id": case.case_id,
                "particle_count": case.particle_count, "family": case.family,
                "target_level": case.target_level, "particle_index": index,
                "lambda_i": value, "lambda_i_over_max": value / audited_lambda,
                "is_first_argmax": index == int(stats["first_argmax"]),
                "local_coupling_sha256": stats["local_coupling_sha256"],
            })
        for prediction in frozen_external_predictions(case.case_id, audited_lambda, rho_l1):
            predictions.append({
                "scale_order": case.scale_order, "case_id": case.case_id,
                "particle_count": case.particle_count, "family": case.family,
                "target_level": case.target_level, "model": prediction.model,
                "predictor_name": "lambda_max" if prediction.model == "M1" else "rho_l1",
                "predictor_value": audited_lambda if prediction.model == "M1" else rho_l1,
                "point_prediction": prediction.point_prediction,
                "safety_factor": prediction.safety_factor,
                "conservative_prediction": prediction.conservative_prediction,
                "safe_1pct": prediction.safe_1pct, "safe_5pct": prediction.safe_5pct,
                "safe_10pct": prediction.safe_10pct,
                "model_e_response_column_present": False,
            })
    m1 = [row for row in predictions if row["model"] == "M1"]
    counts = tuple(sum(bool(row[field]) for row in m1) for field in ("safe_1pct", "safe_5pct", "safe_10pct"))
    if counts != (6, 12, 18):
        raise RuntimeError(f"blind M1 safe-count checksum changed: {counts}")
    _atomic_write(MANIFEST, manifests)
    _atomic_write(LOCAL, locals_)
    _atomic_write(PREDICTIONS, predictions)
    _atomic_write(PRIOR_HASHES, _prior_rows())

    protocol: list[dict[str, object]] = []
    def add(category: str, key: str, value: object, source: str = "T14.1 preregistration") -> None:
        protocol.append({"category": category, "key": key, "value": value, "source": source})
    frozen = (
        ("base_commit", BASE_COMMIT), ("m1_prefactor", M1_PREFACTOR),
        ("m1_exponent", M1_EXPONENT), ("m1_safety_factor", M1_SAFETY_FACTOR),
        ("p3_prefactor", P3_PREFACTOR), ("p3_exponent", P3_EXPONENT),
        ("p3_safety_factor", P3_SAFETY_FACTOR),
        ("comparison", "conservative_prediction < tolerance and observed_error < tolerance"),
        ("radius", RADIUS), ("energy_density", ENERGY), ("k", K), ("ka", KA),
        ("f0", F0), ("f1", F1), ("convergence_tolerance", 1e-5),
        ("minimum_lmax", 5), ("maximum_lmax", 13), ("workers_default", 1),
        ("blas_threads_default", 1), ("required_eligible_all", 20),
        ("required_eligible_each_n", 10), ("required_eligible_each_family", 6),
        ("required_eligible_each_level", 5), ("required_global_factor2", 0.80),
        ("required_each_n_factor2", 0.75), ("required_spearman", 0.90),
        ("required_max_rmse_log", float(np.log(2.0))),
        ("trend_no_deterioration_median_max", 1.10),
        ("trend_no_deterioration_p90_max", 1.25),
        ("trend_systematic_median_min", 1.25),
        ("trend_systematic_count_over_1_10", 9),
        ("trend_minimum_applicable_pairs", 10),
        ("m2_status", "excluded_UNSTABLE_COLLINEARITY"),
        ("phase_a_t14_1_model_e_solves", 0),
    )
    for key, value in frozen:
        add("frozen_rule", key, value)
    for index, (tolerance, lambda_threshold, rho_threshold) in enumerate(zip(TOLERANCES, LAMBDA_THRESHOLDS, RHO_THRESHOLDS), start=1):
        add("threshold", f"tolerance_{index}", tolerance)
        add("threshold", f"lambda_max_{index}", lambda_threshold)
        add("threshold", f"rho_l1_{index}", rho_threshold)
    for index, target in enumerate(LAMBDA_TARGETS, start=1):
        add("selection", f"lambda_target_{index}", target)
    add("selection", "case_ids", ";".join(EXPECTED_LARGE_N_CASE_IDS))
    add("selection", "case_count", 24)
    add("geometry", "particle_counts", "45;105")
    add("geometry", "families", "linear;compact;irregular")
    add("geometry", "compact_rows", "N45:q9;N105:q14")
    add("geometry", "irregular_amplitude", IRREGULAR_AMPLITUDE)
    add("geometry", "irregular_rule", "T14 deterministic perturbation; recenter; renormalize")
    add("resource", "memory_estimate", "16*n_active^2 bytes per dense complex matrix times conservative factor 12")
    for key, value in _hardware().items():
        add("environment", key, value)
    add("environment", "python", platform.python_version())
    add("environment", "numpy", np.__version__)
    add("environment", "scipy", scipy.__version__)
    add("environment", "matplotlib", matplotlib.__version__)
    add("environment", "blas_lapack", json.dumps(np.__config__.CONFIG.get("Build Dependencies", {}), sort_keys=True, default=str))
    for relative in CODE_PATHS:
        add("phase_a_code_sha256", relative, sha256((ROOT / relative).read_bytes()).hexdigest(), relative)
    for relative in (
        "results/data/t12_3_logo_coefficients.csv", "results/data/t12_3_nested_safety_factors.csv",
        "results/data/t12_3_gate.csv", "results/data/t13_frozen_protocol.csv",
        "results/data/t13_frozen_predictions.csv", "results/data/t13_gate.csv",
        "results/data/t14_frozen_protocol.csv", "results/data/t14_frozen_predictions.csv",
        "results/data/t14_gate.csv",
    ):
        add("provenance_sha256", relative, sha256(_git_bytes(BASE_COMMIT, relative)).hexdigest(), relative)
    _atomic_write(PROTOCOL, protocol)
    print(f"T14.1 preregistered {len(cases)} cases; blind M1 safe counts={counts}")


if __name__ == "__main__":
    preregister()
