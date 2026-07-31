#!/usr/bin/env python3
"""Materialize the response-blind T13 external-validation preregistration."""

from __future__ import annotations

import csv
from hashlib import sha256
import platform
from pathlib import Path
import tempfile

import matplotlib
import numpy as np
import scipy

from acoustic_ms import cluster_family, maximum_geometric_coupling
from acoustic_ms.external_validation import (
    EXPECTED_CASE_IDS,
    EXTERNAL_STRATA,
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
    canonical_coordinate_hash,
    frozen_external_predictions,
    select_external_validation_cases,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
T08_CASES = DATA / "t08_cases.csv"
MANIFEST = DATA / "t13_holdout_manifest.csv"
PREDICTIONS = DATA / "t13_frozen_predictions.csv"
PROTOCOL = DATA / "t13_frozen_protocol.csv"
PHASE_B_PATHS = (
    DATA / "t13_model_e_convergence.csv",
    DATA / "t13_forces.csv",
    DATA / "t13_case_summary.csv",
    DATA / "t13_external_predictions.csv",
    DATA / "t13_metrics.csv",
    DATA / "t13_threshold_audit.csv",
    DATA / "t13_gate.csv",
    ROOT / "results" / "figures" / "t13_external_validation.png",
)
RADIUS = 1.0
K = 0.1
KA = 0.1
ENERGY_DENSITY = 1.0
F0 = 0.0
METADATA_FIELDS = (
    "case_id", "split", "particle_count", "family", "radius", "k", "ka",
    "energy_density", "f0", "f1", "distance_ratio", "reference_lmax",
    "total_converged", "lambda_max", "rho_l1",
)


def _format(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    return value


def _atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot publish empty table {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise ValueError("all preregistration rows must share exact field order")
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(value) for key, value in row.items()})
    temporary.replace(path)


def _read_response_blind_metadata(path: Path = T08_CASES) -> list[dict[str, str]]:
    """Load only the explicitly allowed T08 metadata/predictor columns."""

    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        if len(header) != len(set(header)):
            raise RuntimeError("T08 header contains duplicate fields")
        missing = [field for field in METADATA_FIELDS if field not in header]
        if missing:
            raise RuntimeError(f"T08 metadata fields missing: {missing}")
        indices = {field: header.index(field) for field in METADATA_FIELDS}
        rows = [
            {field: values[index] for field, index in indices.items()}
            for values in reader
        ]
    if len(rows) != 312:
        raise RuntimeError("the frozen T08 table must contain 312 metadata rows")
    return rows


def _serialize(positions: np.ndarray) -> str:
    return ";".join(
        ":".join(format(float(value), ".17g") for value in row)
        for row in positions
    )


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _validate_physics(rows: list[dict[str, str]]) -> None:
    for row in rows:
        if row["split"] != "holdout":
            continue
        checks = (
            float(row["radius"]) == RADIUS,
            float(row["k"]) == K,
            float(row["ka"]) == KA,
            float(row["energy_density"]) == ENERGY_DENSITY,
            float(row["f0"]) == F0,
            int(row["particle_count"]) in (6, 10),
            row["family"] in ("linear", "compact", "irregular"),
        )
        if not all(checks):
            raise RuntimeError(f"frozen physical metadata changed for {row['case_id']}")


def preregister() -> None:
    if any(path.exists() for path in PHASE_B_PATHS):
        existing = [path.name for path in PHASE_B_PATHS if path.exists()]
        raise RuntimeError(f"response-bearing T13 files exist before preregistration: {existing}")
    metadata = _read_response_blind_metadata()
    _validate_physics(metadata)
    cases = select_external_validation_cases(metadata)
    if tuple(case.case_id for case in cases) != EXPECTED_CASE_IDS:
        raise RuntimeError("the exact 24-case checksum failed")
    manifest_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for order, case in enumerate(cases, start=1):
        positions = cluster_family(case.particle_count, case.family, case.distance_ratio)
        audited_lambda = maximum_geometric_coupling(positions, RADIUS, case.f1)
        if not np.isclose(audited_lambda, case.lambda_max, rtol=3e-14, atol=3e-15):
            raise RuntimeError(f"Lambda_max audit failed for {case.case_id}")
        manifest_rows.append({
            "holdout_order": order,
            "case_id": case.case_id,
            "particle_count": case.particle_count,
            "family": case.family,
            "stratum": case.stratum,
            "target_level": case.target_level,
            "lambda_target": case.lambda_target,
            "f1": case.f1,
            "distance_ratio": case.distance_ratio,
            "radius": RADIUS,
            "k": K,
            "ka": KA,
            "energy_density": ENERGY_DENSITY,
            "f0": F0,
            "coordinates_xyz": _serialize(positions),
            "coordinate_sha256": canonical_coordinate_hash(positions),
            "lambda_max": case.lambda_max,
            "lambda_public_audit": audited_lambda,
            "rho_l1": case.rho_l1,
            "reference_lmax_d": case.reference_lmax,
            "reference_d_confirmed": True,
            "provenance": "pre_registered_t08_holdout_metadata",
            "model_e_response_consulted": False,
        })
        for prediction in frozen_external_predictions(
            case.case_id, case.lambda_max, case.rho_l1
        ):
            prediction_rows.append({
                "holdout_order": order,
                "case_id": case.case_id,
                "particle_count": case.particle_count,
                "family": case.family,
                "stratum": case.stratum,
                "target_level": case.target_level,
                "model": prediction.model,
                "predictor_name": "lambda_max" if prediction.model == "M1" else "rho_l1",
                "predictor_value": case.lambda_max if prediction.model == "M1" else case.rho_l1,
                "point_prediction": prediction.point_prediction,
                "safety_factor": prediction.safety_factor,
                "conservative_prediction": prediction.conservative_prediction,
                "safe_1pct": prediction.safe_1pct,
                "safe_5pct": prediction.safe_5pct,
                "safe_10pct": prediction.safe_10pct,
                "model_e_response_column_present": False,
            })
    m1 = [row for row in prediction_rows if row["model"] == "M1"]
    counts = tuple(sum(bool(row[field]) for row in m1) for field in ("safe_1pct", "safe_5pct", "safe_10pct"))
    if counts != (6, 12, 18):
        raise RuntimeError(f"blind M1 safety checksum changed: {counts}")
    _atomic_write(MANIFEST, manifest_rows)
    _atomic_write(PREDICTIONS, prediction_rows)

    provenance_paths = (
        DATA / "t12_3_logo_coefficients.csv",
        DATA / "t12_3_nested_safety_factors.csv",
        DATA / "t12_3_gate.csv",
        ROOT / "TAREFA_T12_3_CRITERIO_MECANISTICO_VALIDACAO_AGRUPADA.md",
        T08_CASES,
    )
    protocol_rows: list[dict[str, object]] = []
    def add(category: str, key: str, value: object, source: str = "T13 preregistration") -> None:
        protocol_rows.append({"category": category, "key": key, "value": value, "source": source})
    for key, value in (
        ("m1_prefactor", M1_PREFACTOR), ("m1_exponent", M1_EXPONENT),
        ("m1_safety_factor", M1_SAFETY_FACTOR), ("p3_prefactor", P3_PREFACTOR),
        ("p3_exponent", P3_EXPONENT), ("p3_safety_factor", P3_SAFETY_FACTOR),
        ("comparison", "conservative_prediction < tolerance and observed_error < tolerance"),
        ("convergence_tolerance", 1e-5), ("minimum_lmax", 5),
        ("standard_lmax_cap", 13), ("extension_lmax_cap", 21),
        ("required_eligible_all", 20), ("required_eligible_n6", 10),
        ("required_eligible_n10", 10), ("required_global_factor2", 0.80),
        ("required_each_n_factor2", 0.75), ("required_spearman", 0.90),
        ("required_max_rmse_log", float(np.log(2.0))),
        ("m2_status", "excluded_UNSTABLE_COLLINEARITY"),
        ("phase_a_model_e_solves", 0),
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
    add("selection", "strata", ";".join(EXTERNAL_STRATA))
    add("selection", "case_count", len(cases))
    add("selection", "assignment", "minimum total absolute log-distance; lexicographic tie at 1e-15")
    add("environment", "python", platform.python_version())
    add("environment", "numpy", np.__version__)
    add("environment", "scipy", scipy.__version__)
    add("environment", "matplotlib", matplotlib.__version__)
    for path in provenance_paths:
        add("provenance_sha256", str(path.relative_to(ROOT)), _hash(path), str(path.relative_to(ROOT)))
    _atomic_write(PROTOCOL, protocol_rows)
    print(f"T13 preregistered {len(cases)} cases; blind M1 safe counts={counts}")


if __name__ == "__main__":
    preregister()
