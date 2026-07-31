#!/usr/bin/env python3
"""Publish the response-blind T14 scale-out preregistration."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path
import platform
import subprocess
import tempfile

import matplotlib
import numpy as np
import scipy

from acoustic_ms import (
    EXPECTED_SCALE_OUT_CASE_IDS,
    build_scale_out_cases,
    canonical_coordinate_hash,
    frozen_external_predictions,
    maximum_geometric_coupling,
    solve_multipolar_nodal_interaction_forces,
    spectral_radius_l1,
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
BASE_COMMIT = "a0e6d97c45915e5974f10d5f5ee9523070d9876f"
MANIFEST = DATA / "t14_scale_manifest.csv"
PREDICTIONS = DATA / "t14_frozen_predictions.csv"
PROTOCOL = DATA / "t14_frozen_protocol.csv"
PRIOR_HASHES = DATA / "t14_prior_artifact_hashes.csv"
PHASE_B_PATHS = (
    DATA / "t14_model_e_convergence.csv",
    DATA / "t14_forces.csv",
    DATA / "t14_case_summary.csv",
    DATA / "t14_scale_predictions.csv",
    DATA / "t14_metrics.csv",
    DATA / "t14_threshold_audit.csv",
    DATA / "t14_matched_scale_pairs.csv",
    DATA / "t14_performance.csv",
    DATA / "t14_gate.csv",
    FIGURES / "t14_scale_out_validation.png",
)
RADIUS = 1.0
K = 0.1
KA = 0.1
ENERGY_DENSITY = 1.0
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
        raise ValueError("preregistration rows must share exact field order")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(value) for key, value in row.items()})
    temporary.replace(path)


def _serialize(positions: np.ndarray) -> str:
    return ";".join(
        ":".join(format(float(value), ".17g") for value in row)
        for row in positions
    )


def _git_bytes(commit: str, relative: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{commit}:{relative}"], cwd=ROOT, check=True,
        capture_output=True,
    ).stdout


def _prior_artifact_rows() -> list[dict[str, object]]:
    paths = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", BASE_COMMIT,
         "results/data", "results/figures"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    rows = []
    for relative in paths:
        payload = _git_bytes(BASE_COMMIT, relative)
        rows.append({
            "path": relative,
            "size_bytes": len(payload),
            "sha256": sha256(payload).hexdigest(),
            "source_commit": BASE_COMMIT,
        })
    if len(rows) != 81:
        raise RuntimeError(f"expected 81 prior artifacts, found {len(rows)}")
    return rows


def preregister() -> None:
    if subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip() != BASE_COMMIT:
        raise RuntimeError("T14 preregistration must start at the frozen T13 commit")
    existing = [str(path.relative_to(ROOT)) for path in PHASE_B_PATHS if path.exists()]
    if existing:
        raise RuntimeError(f"response-bearing T14 files exist before preregistration: {existing}")
    cases = build_scale_out_cases(f1=F1)
    if tuple(case.case_id for case in cases) != EXPECTED_SCALE_OUT_CASE_IDS:
        raise RuntimeError("the exact T14 24-case checksum failed")
    manifest_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for case in cases:
        positions = case.positions_xyz
        if case.minimum_distance <= 2.0 * RADIUS:
            raise RuntimeError(f"overlap margin failed for {case.case_id}")
        audited_lambda = maximum_geometric_coupling(positions, RADIUS, F1)
        if not np.isclose(
            audited_lambda, case.lambda_target, rtol=5.0e-13, atol=5.0e-15
        ):
            raise RuntimeError(f"Lambda_max audit failed for {case.case_id}")
        l1 = solve_multipolar_nodal_interaction_forces(
            positions, K, RADIUS, ENERGY_DENSITY, F0, F1, 1
        )
        rho_l1 = spectral_radius_l1(l1)
        manifest_rows.append({
            "scale_order": case.scale_order,
            "case_id": case.case_id,
            "particle_count": case.particle_count,
            "family": case.family,
            "target_level": case.target_level,
            "radius": RADIUS,
            "k": K,
            "ka": KA,
            "energy_density": ENERGY_DENSITY,
            "f0": F0,
            "f1": F1,
            "geometry_factor": case.geometry_factor,
            "lambda_target": case.lambda_target,
            "distance_ratio": case.distance_ratio,
            "minimum_distance": case.minimum_distance,
            "overlap_margin": case.minimum_distance - 2.0 * RADIUS,
            "cluster_diameter": case.cluster_diameter,
            "k_cluster_diameter": K * case.cluster_diameter,
            "coordinates_xyz": _serialize(positions),
            "coordinate_sha256": canonical_coordinate_hash(positions),
            "lambda_max": audited_lambda,
            "rho_l1": rho_l1,
            "provenance": "pre_registered_analytic_scale_out_geometry",
            "model_e_response_consulted": False,
        })
        for prediction in frozen_external_predictions(
            case.case_id, audited_lambda, rho_l1
        ):
            prediction_rows.append({
                "scale_order": case.scale_order,
                "case_id": case.case_id,
                "particle_count": case.particle_count,
                "family": case.family,
                "target_level": case.target_level,
                "model": prediction.model,
                "predictor_name": "lambda_max" if prediction.model == "M1" else "rho_l1",
                "predictor_value": audited_lambda if prediction.model == "M1" else rho_l1,
                "point_prediction": prediction.point_prediction,
                "safety_factor": prediction.safety_factor,
                "conservative_prediction": prediction.conservative_prediction,
                "safe_1pct": prediction.safe_1pct,
                "safe_5pct": prediction.safe_5pct,
                "safe_10pct": prediction.safe_10pct,
                "model_e_response_column_present": False,
            })
    m1 = [row for row in prediction_rows if row["model"] == "M1"]
    counts = tuple(
        sum(bool(row[field]) for row in m1)
        for field in ("safe_1pct", "safe_5pct", "safe_10pct")
    )
    if counts != (6, 12, 18):
        raise RuntimeError(f"blind M1 safe-count checksum changed: {counts}")
    expected_points = (
        0.0014863937048093015, 0.0087007541172098232,
        0.027514200371470676, 0.10186146774330709,
    )
    expected_safe = (
        0.0038199876936037326, 0.02236067977499789,
        0.070710678118654738, 0.26178094805761298,
    )
    for level in range(1, 5):
        rows = [row for row in m1 if row["target_level"] == level]
        if len(rows) != 6 or any(
            not np.isclose(float(row["point_prediction"]), expected_points[level - 1], rtol=3e-15, atol=3e-17)
            or not np.isclose(float(row["conservative_prediction"]), expected_safe[level - 1], rtol=3e-15, atol=3e-17)
            for row in rows
        ):
            raise RuntimeError(f"blind prediction checksum failed at level {level}")
    _atomic_write(MANIFEST, manifest_rows)
    _atomic_write(PREDICTIONS, prediction_rows)
    _atomic_write(PRIOR_HASHES, _prior_artifact_rows())

    protocol_rows: list[dict[str, object]] = []
    def add(category: str, key: str, value: object, source: str = "T14 preregistration") -> None:
        protocol_rows.append({"category": category, "key": key, "value": value, "source": source})
    for key, value in (
        ("base_commit", BASE_COMMIT),
        ("m1_prefactor", M1_PREFACTOR), ("m1_exponent", M1_EXPONENT),
        ("m1_safety_factor", M1_SAFETY_FACTOR), ("p3_prefactor", P3_PREFACTOR),
        ("p3_exponent", P3_EXPONENT), ("p3_safety_factor", P3_SAFETY_FACTOR),
        ("comparison", "conservative_prediction < tolerance and observed_error < tolerance"),
        ("radius", RADIUS), ("energy_density", ENERGY_DENSITY), ("k", K),
        ("ka", KA), ("f0", F0), ("f1", F1),
        ("convergence_tolerance", 1.0e-5), ("minimum_lmax", 5),
        ("maximum_lmax", 13), ("workers_default", 1),
        ("required_eligible_all", 20), ("required_eligible_each_n", 10),
        ("required_eligible_each_family", 6), ("required_eligible_each_level", 5),
        ("required_global_factor2", 0.80), ("required_each_n_factor2", 0.75),
        ("required_spearman", 0.90), ("required_max_rmse_log", float(np.log(2.0))),
        ("m2_status", "excluded_UNSTABLE_COLLINEARITY"),
        ("phase_a_t14_model_e_solves", 0),
    ):
        add("frozen_rule", key, value)
    for index, (tolerance, lambda_threshold, rho_threshold) in enumerate(
        zip(TOLERANCES, LAMBDA_THRESHOLDS, RHO_THRESHOLDS), start=1
    ):
        add("threshold", f"tolerance_{index}", tolerance)
        add("threshold", f"lambda_max_{index}", lambda_threshold)
        add("threshold", f"rho_l1_{index}", rho_threshold)
    for index, target in enumerate(LAMBDA_TARGETS, start=1):
        add("selection", f"lambda_target_{index}", target)
    add("geometry", "particle_counts", "15;28")
    add("geometry", "families", "linear;compact;irregular")
    add("geometry", "compact", "triangular rows q=5 or q=7; centroid zero; minimum distance one")
    add("geometry", "irregular_amplitude", IRREGULAR_AMPLITUDE)
    add("geometry", "irregular_rule", "explicit frac((i+1)*sqrt(2/3)) deterministic perturbation; recenter; renormalize")
    add("selection", "case_ids", ";".join(EXPECTED_SCALE_OUT_CASE_IDS))
    add("selection", "case_count", 24)
    add("environment", "python", platform.python_version())
    add("environment", "numpy", np.__version__)
    add("environment", "scipy", scipy.__version__)
    add("environment", "matplotlib", matplotlib.__version__)
    for relative in (
        "results/data/t12_3_logo_coefficients.csv",
        "results/data/t12_3_nested_safety_factors.csv",
        "results/data/t12_3_gate.csv",
        "results/data/t13_frozen_protocol.csv",
        "results/data/t13_frozen_predictions.csv",
        "results/data/t13_gate.csv",
        "TAREFA_T13_VALIDACAO_EXTERNA_LAMBDA_MAX.md",
    ):
        add("provenance_sha256", relative, sha256(_git_bytes(BASE_COMMIT, relative)).hexdigest(), relative)
    _atomic_write(PROTOCOL, protocol_rows)
    print(f"T14 preregistered {len(cases)} cases; blind M1 safe counts={counts}")


if __name__ == "__main__":
    preregister()
