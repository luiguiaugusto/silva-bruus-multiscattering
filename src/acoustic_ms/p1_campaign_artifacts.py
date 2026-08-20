"""Pure, deterministic P1.6 artifact generation and frozen G1 evaluation.

This module deliberately imports no solver or other :mod:`acoustic_ms` module.
It can therefore regenerate every campaign derivative from checkpoints in a
fresh process without making scientific code importable, let alone callable.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence


G1_BUDGET = 1.0e-12
CAMPAIGN_ID = "p1_dimer_confirmatory"
ARTIFACT_PATHS = {
    "data_raw.csv": "campaigns/p1/data_raw.csv",
    "data_derived.csv": "campaigns/p1/data_derived.csv",
    "data_plot.csv": "campaigns/p1/data_plot.csv",
    "failures.csv": "campaigns/p1/failures.csv",
    "performance.csv": "campaigns/p1/performance.csv",
}
CHANNELS = (
    "total",
    "interaction",
    "external_scattered",
    "scattered_scattered",
)


class CampaignArtifactError(RuntimeError):
    """Raised when checkpoint data cannot produce trustworthy artifacts."""


@dataclass(frozen=True)
class G1Result:
    """Frozen P1 confirmatory gate result."""

    gate_status: str
    decision: str
    attempted_count: int
    eligible_count: int
    eligible_primary_count: int
    eligible_audit_count: int
    covered_strata: tuple[tuple[float, str], ...]
    missing_strata: tuple[tuple[float, str], ...]
    eligible_audit_twin_pairs: int
    identity_error_max: float | None
    rotation_error_max: float | None
    reasons: tuple[str, ...]


def _format(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".17g")
    return value


def _csv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    field_tuple = tuple(fields)
    row_list = list(rows)
    expected = set(field_tuple)
    for row in row_list:
        if set(row) != expected:
            raise CampaignArtifactError(
                "campaign CSV row differs from its frozen field contract"
            )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=field_tuple,
        lineterminator="\n",
    )
    writer.writeheader()
    for row in row_list:
        writer.writerow({field: _format(row[field]) for field in field_tuple})
    return stream.getvalue().encode("utf-8")


def _vector_norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(float(component) ** 2 for component in vector))


def _vector_error(
    observed: Sequence[Sequence[float]],
    reference: Sequence[Sequence[float]],
) -> float:
    """Return the frozen scale-safe maximum vector discrepancy.

    Forces are already normalized by ``a^2 E0``.  The denominator is the
    larger of one and all compared vector magnitudes, so the gate is absolute
    below unit scale and relative above it.
    """

    if len(observed) != len(reference) or not observed:
        raise CampaignArtifactError("force arrays must have equal nonzero length")
    if any(len(vector) != 3 for vector in (*observed, *reference)):
        raise CampaignArtifactError("force vectors must have three components")
    residuals = [
        _vector_norm([float(left) - float(right) for left, right in zip(a, b)])
        for a, b in zip(observed, reference)
    ]
    scale = max(
        1.0,
        *(_vector_norm(vector) for vector in observed),
        *(_vector_norm(vector) for vector in reference),
    )
    return max(residuals) / scale


def _rotate_planar(
    vectors: Sequence[Sequence[float]], angle: float
) -> list[list[float]]:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return [
        [
            cosine * float(vector[0]) - sine * float(vector[1]),
            sine * float(vector[0]) + cosine * float(vector[1]),
            float(vector[2]),
        ]
        for vector in vectors
    ]


def _completed(record: Mapping[str, Any]) -> bool:
    return record.get("state") == "completed" and isinstance(
        record.get("outcome"), Mapping
    )


def _eligible(record: Mapping[str, Any]) -> bool:
    return _completed(record) and bool(record["outcome"].get("eligible"))


def _identity_error(record: Mapping[str, Any]) -> float | None:
    if not _eligible(record):
        return None
    outcome = record["outcome"]
    return _vector_error(
        outcome["model_be_forces_xyz"],
        outcome["model_e_forces_xyz"],
    )


def _rotation_error(
    audit: Mapping[str, Any], twin: Mapping[str, Any]
) -> float | None:
    if not (_eligible(audit) and _eligible(twin)):
        return None
    audit_parameters = audit["parameters"]
    twin_parameters = twin["parameters"]
    angle = float(audit_parameters["theta_rad"]) - float(
        twin_parameters["theta_rad"]
    )
    expected = _rotate_planar(
        twin["outcome"]["model_be_forces_xyz"], angle
    )
    return _vector_error(audit["outcome"]["model_be_forces_xyz"], expected)


def evaluate_g1(
    manifest: Mapping[str, Any],
    records: Sequence[Mapping[str, Any]],
    *,
    contract_valid: bool = True,
) -> G1Result:
    """Evaluate the response-blind G1 gate frozen by P1.6A."""

    cases = list(manifest["cases"])
    by_record_id = {str(record["case_id"]): record for record in records}
    if len(by_record_id) != len(records):
        contract_valid = False
    expected_identity = [
        (int(case["case_order"]), case["case_id"]) for case in cases
    ]
    observed_identity = [
        (int(record["case_order"]), record["case_id"]) for record in records
    ]
    if observed_identity != expected_identity:
        contract_valid = False
    for case, record in zip(cases, records):
        if record.get("parameters") != case.get("parameters"):
            contract_valid = False
        if int(record.get("attempt_count", -1)) not in (0, 1):
            contract_valid = False
        outcome = record.get("outcome")
        if isinstance(outcome, Mapping):
            evaluated = list(outcome.get("evaluated_lmax", []))
            order_values = [
                int(order.get("lmax", -1)) for order in outcome.get("orders", [])
            ]
            if (
                order_values != evaluated
                or len(evaluated) != len(set(evaluated))
                or int(outcome.get("model_e_solve_count", -1)) != len(evaluated)
                or any(
                    set(order.get("channels", {})) != set(CHANNELS)
                    for order in outcome.get("orders", [])
                )
            ):
                contract_valid = False
    attempted = [
        record for record in records if int(record.get("attempt_count", 0)) == 1
    ]
    eligible = [record for record in records if _eligible(record)]
    primaries = [
        record
        for record in eligible
        if record["parameters"]["evidence_role"] == "primary"
    ]
    audits = [
        record
        for record in records
        if record["parameters"]["evidence_role"] == "rotational_audit"
    ]

    required_strata = sorted(
        {
            (float(case["parameters"]["ka"]), case["parameters"]["material_id"])
            for case in cases
            if case["parameters"]["evidence_role"] == "primary"
        }
    )
    covered = sorted(
        {
            (
                float(record["parameters"]["ka"]),
                record["parameters"]["material_id"],
            )
            for record in primaries
        }
    )
    missing = sorted(set(required_strata) - set(covered))

    identity_errors = [
        error
        for error in (_identity_error(record) for record in eligible)
        if error is not None
    ]
    rotation_errors: list[float] = []
    eligible_pairs = 0
    for audit in audits:
        twin_id = audit["parameters"]["twin_case_id"]
        twin = by_record_id.get(twin_id)
        if twin is None:
            continue
        error = _rotation_error(audit, twin)
        if error is not None:
            eligible_pairs += 1
            rotation_errors.append(error)

    identity_max = max(identity_errors) if identity_errors else None
    rotation_max = max(rotation_errors) if rotation_errors else None
    symmetry_or_contract_reasons: list[str] = []
    if not contract_valid:
        symmetry_or_contract_reasons.append("campaign_contract_invalid")
    if any(record.get("failure_stage") == "contract" for record in records):
        symmetry_or_contract_reasons.append("case_contract_invalid")
    if identity_max is not None and identity_max > G1_BUDGET:
        symmetry_or_contract_reasons.append("be_e_identity_budget_exceeded")
    if rotation_max is not None and rotation_max > G1_BUDGET:
        symmetry_or_contract_reasons.append("rotational_covariance_budget_exceeded")

    coverage_reasons: list[str] = []
    if len(attempted) != len(cases) or len(cases) != 102:
        coverage_reasons.append("not_all_102_cases_attempted")
    if missing:
        coverage_reasons.append("eligible_primary_stratum_coverage_incomplete")
    if eligible_pairs != 6:
        coverage_reasons.append("eligible_audit_twin_coverage_incomplete")
    if any(
        record.get("failure_stage") in {"timeout", "memory", "global_limit"}
        for record in records
    ):
        coverage_reasons.append("resource_exhaustion")

    if symmetry_or_contract_reasons:
        gate_status = "FAIL_G1"
        decision = "NO_GO_P2"
        reasons = symmetry_or_contract_reasons + coverage_reasons
    elif coverage_reasons:
        gate_status = "INCONCLUSIVE_G1"
        decision = "INCONCLUSIVE_P1"
        reasons = coverage_reasons
    else:
        gate_status = "PASS_G1"
        decision = "GO_P2"
        reasons = ["all_frozen_g1_gates_passed"]

    return G1Result(
        gate_status=gate_status,
        decision=decision,
        attempted_count=len(attempted),
        eligible_count=len(eligible),
        eligible_primary_count=len(primaries),
        eligible_audit_count=sum(_eligible(record) for record in audits),
        covered_strata=tuple(covered),
        missing_strata=tuple(missing),
        eligible_audit_twin_pairs=eligible_pairs,
        identity_error_max=identity_max,
        rotation_error_max=rotation_max,
        reasons=tuple(reasons),
    )


def load_checkpoint_records(state_directory: str | Path) -> tuple[dict[str, Any], ...]:
    """Load ordered case checkpoints without importing a solver."""

    directory = Path(state_directory)
    ledger_path = directory / "campaign_ledger.json"
    try:
        ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignArtifactError(f"cannot load campaign ledger: {exc}") from exc
    records: list[dict[str, Any]] = []
    for entry in ledger.get("cases", []):
        checkpoint = directory / "cases" / f"{int(entry['case_order']):03d}.json"
        if checkpoint.exists():
            record = json.loads(checkpoint.read_text(encoding="utf-8"))
        else:
            record = dict(entry)
            record["parameters"] = entry["parameters"]
            record["outcome"] = None
        records.append(record)
    return tuple(records)


_RAW_FIELDS = (
    "schema_version", "campaign_id", "manifest_sha256", "git_commit",
    "case_order", "case_id", "classification", "evidence_role",
    "include_in_scientific_tables", "attempt_state", "attempt_count",
    "started_utc", "completed_utc", "radius_m", "k_rad_m", "ka",
    "energy_density_j_m3", "material_id", "material_model", "f0",
    "f0_applicable", "f1", "distance_ratio", "theta_rad",
    "particle_count", "position_x_m", "position_y_m", "position_z_m",
    "worker_count", "blas_threads",
    "model", "force_channel", "lmax", "particle_index",
    "order_wall_seconds", "order_peak_rss_bytes", "diagnostics_json",
    "force_x_over_a2e0", "force_y_over_a2e0", "force_z_over_a2e0",
    "successive_change", "absolute_change", "change_applicable",
    "confirmed", "eligible", "failure_stage", "ineligibility_reason",
)
_DERIVED_FIELDS = (
    "schema_version", "campaign_id", "manifest_sha256", "case_order",
    "case_id", "classification", "evidence_role", "eligible",
    "include_in_scientific_tables", "metric", "value", "unit",
    "applicable", "reason",
)
_PLOT_FIELDS = (
    "schema_version", "campaign_id", "case_order", "case_id", "series_id",
    "point_order", "x_name", "x_value", "x_unit", "y_name", "y_value",
    "y_unit", "eligible", "annotation",
)
_FAILURE_FIELDS = (
    "schema_version", "campaign_id", "manifest_sha256", "case_order",
    "case_id", "attempt_state", "attempt_count", "eligible",
    "failure_stage", "failure_reason", "include_in_scientific_tables",
)
_PERFORMANCE_FIELDS = (
    "schema_version", "campaign_id", "manifest_sha256", "case_order",
    "case_id", "attempt_state", "attempt_count", "started_utc",
    "completed_utc", "attempted_lmax", "evaluated_lmax", "final_lmax",
    "case_wall_seconds", "peak_rss_bytes", "model_e_solve_count",
    "wall_seconds_per_case", "peak_rss_bytes_per_case",
    "wall_seconds_campaign", "worker_count", "blas_threads",
    "final_diagnostics_json",
    "converged", "eligible", "failure_stage", "failure_reason",
)


def _base(record: Mapping[str, Any], manifest: Mapping[str, Any]) -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "campaign_id": manifest["campaign_id"],
        "manifest_sha256": manifest["provenance"]["manifest_sha256"],
        "case_order": int(record["case_order"]),
        "case_id": record["case_id"],
    }


def _positions(
    parameters: Mapping[str, Any], radius: float
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    distance = radius * float(parameters["distance_ratio"])
    theta = float(parameters["theta_rad"])
    displacement = (distance * math.cos(theta), distance * math.sin(theta), 0.0)
    first = tuple(-0.5 * value for value in displacement)
    second = tuple(0.5 * value for value in displacement)
    return first, second


def build_campaign_artifacts(
    manifest: Mapping[str, Any],
    records: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, bytes], G1Result]:
    """Build all five deterministic P1.6 tables from checkpoint records."""

    ordered = sorted(records, key=lambda record: int(record["case_order"]))
    g1 = evaluate_g1(manifest, ordered)
    raw_rows: list[dict[str, object]] = []
    derived_rows: list[dict[str, object]] = []
    plot_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    performance_rows: list[dict[str, object]] = []

    for record in ordered:
        base = _base(record, manifest)
        parameters = record["parameters"]
        outcome = record.get("outcome") or {}
        eligible = _eligible(record)
        common = {
            **base,
            "git_commit": manifest["provenance"]["git_commit"],
            "classification": manifest["classification"],
            "evidence_role": parameters["evidence_role"],
            "include_in_scientific_tables": bool(
                parameters["include_in_scientific_tables"]
            ),
            "attempt_state": record["state"],
            "attempt_count": int(record.get("attempt_count", 0)),
            "started_utc": record.get("started_utc"),
            "completed_utc": record.get("completed_utc"),
            "radius_m": float(manifest["physical"]["radius_m"]),
            "k_rad_m": float(parameters["k_rad_m"]),
            "ka": float(parameters["ka"]),
            "energy_density_j_m3": float(
                manifest["physical"]["energy_density_j_m3"]
            ),
            "material_id": parameters["material_id"],
            "material_model": parameters["material_model"],
            "f0": float(parameters["f0"]),
            "f0_applicable": bool(parameters["f0_applicable"]),
            "f1": float(parameters["f1"]),
            "distance_ratio": float(parameters["distance_ratio"]),
            "theta_rad": float(parameters["theta_rad"]),
            "particle_count": 2,
            "worker_count": int(manifest["resources"]["worker_count"]),
            "blas_threads": int(manifest["resources"]["blas_threads"]),
        }
        positions = _positions(parameters, float(manifest["physical"]["radius_m"]))
        for order in outcome.get("orders", []):
            for channel in CHANNELS:
                channel_record = order["channels"][channel]
                for particle_index, force in enumerate(channel_record["forces_xyz"]):
                    raw_rows.append(
                        {
                            **common,
                            "model": "E",
                            "force_channel": channel,
                            "lmax": int(order["lmax"]),
                            "particle_index": particle_index,
                            "position_x_m": positions[particle_index][0],
                            "position_y_m": positions[particle_index][1],
                            "position_z_m": positions[particle_index][2],
                            "order_wall_seconds": float(order["wall_seconds"]),
                            "order_peak_rss_bytes": int(order["peak_rss_bytes"]),
                            "diagnostics_json": json.dumps(
                                order["diagnostics"],
                                sort_keys=True,
                                separators=(",", ":"),
                                ensure_ascii=True,
                                allow_nan=False,
                            ),
                            "force_x_over_a2e0": float(force[0]),
                            "force_y_over_a2e0": float(force[1]),
                            "force_z_over_a2e0": float(force[2]),
                            "successive_change": float(
                                channel_record["successive_change"]
                            ),
                            "absolute_change": float(
                                channel_record["absolute_change"]
                            ),
                            "change_applicable": bool(
                                channel_record["applicable"]
                            ),
                            "confirmed": bool(channel_record["confirmed"]),
                            "eligible": eligible,
                            "failure_stage": record.get("failure_stage"),
                            "ineligibility_reason": record.get("failure_reason"),
                        }
                    )
        if outcome:
            for model, key in (
                ("A", "model_a_forces_xyz"),
                ("B_E", "model_be_forces_xyz"),
            ):
                forces = outcome.get(key)
                if forces is None:
                    continue
                for particle_index, force in enumerate(forces):
                    raw_rows.append(
                        {
                            **common,
                            "model": model,
                            "force_channel": "interaction",
                            "lmax": outcome.get("final_lmax"),
                            "particle_index": particle_index,
                            "position_x_m": positions[particle_index][0],
                            "position_y_m": positions[particle_index][1],
                            "position_z_m": positions[particle_index][2],
                            "order_wall_seconds": None,
                            "order_peak_rss_bytes": None,
                            "diagnostics_json": "{}",
                            "force_x_over_a2e0": float(force[0]),
                            "force_y_over_a2e0": float(force[1]),
                            "force_z_over_a2e0": float(force[2]),
                            "successive_change": None,
                            "absolute_change": None,
                            "change_applicable": False,
                            "confirmed": bool(outcome.get("converged")),
                            "eligible": eligible,
                            "failure_stage": record.get("failure_stage"),
                            "ineligibility_reason": record.get("failure_reason"),
                        }
                    )

        be_minus_a: float | None = None
        identity = _identity_error(record)
        if eligible:
            be = outcome["model_be_forces_xyz"]
            model_a = outcome["model_a_forces_xyz"]
            differences = [
                [float(x) - float(y) for x, y in zip(left, right)]
                for left, right in zip(be, model_a)
            ]
            be_minus_a = math.sqrt(
                sum(_vector_norm(vector) ** 2 for vector in differences)
                / len(differences)
            )
        for metric, value, unit, applicable, reason in (
            (
                "be_minus_a_rms",
                be_minus_a,
                "a2e0",
                eligible,
                "scientific_magnitude_not_a_gate" if eligible else "case_ineligible",
            ),
            (
                "be_e_identity_error",
                identity,
                "1",
                identity is not None,
                "g1_budget_1e-12" if identity is not None else "case_ineligible",
            ),
        ):
            derived_rows.append(
                {
                    **base,
                    "classification": manifest["classification"],
                    "evidence_role": parameters["evidence_role"],
                    "eligible": eligible,
                    "include_in_scientific_tables": bool(
                        parameters["include_in_scientific_tables"]
                    ),
                    "metric": metric,
                    "value": value,
                    "unit": unit,
                    "applicable": applicable,
                    "reason": reason,
                }
            )
        if eligible and parameters["evidence_role"] == "primary":
            plot_rows.append(
                {
                    "schema_version": "1.0.0",
                    "campaign_id": manifest["campaign_id"],
                    "case_order": int(record["case_order"]),
                    "case_id": record["case_id"],
                    "series_id": (
                        f"ka={format(float(parameters['ka']), '.17g')};"
                        f"material={parameters['material_id']}"
                    ),
                    "point_order": int(record["case_order"]),
                    "x_name": "distance_ratio",
                    "x_value": float(parameters["distance_ratio"]),
                    "x_unit": "1",
                    "y_name": "be_minus_a_rms",
                    "y_value": be_minus_a,
                    "y_unit": "a2e0",
                    "eligible": True,
                    "annotation": "eligible_primary_only",
                }
            )
        if not eligible:
            failure_rows.append(
                {
                    **base,
                    "attempt_state": record["state"],
                    "attempt_count": int(record.get("attempt_count", 0)),
                    "eligible": False,
                    "failure_stage": record.get("failure_stage"),
                    "failure_reason": record.get("failure_reason"),
                    "include_in_scientific_tables": bool(
                        parameters["include_in_scientific_tables"]
                    ),
                }
            )
        performance_rows.append(
            {
                **base,
                "attempt_state": record["state"],
                "attempt_count": int(record.get("attempt_count", 0)),
                "started_utc": record.get("started_utc"),
                "completed_utc": record.get("completed_utc"),
                "attempted_lmax": ";".join(
                    str(value) for value in outcome.get("attempted_lmax", [])
                ),
                "evaluated_lmax": ";".join(
                    str(value) for value in outcome.get("evaluated_lmax", [])
                ),
                "final_lmax": outcome.get("final_lmax"),
                "case_wall_seconds": outcome.get("wall_seconds"),
                "peak_rss_bytes": outcome.get("peak_rss_bytes"),
                "model_e_solve_count": outcome.get("model_e_solve_count"),
                "wall_seconds_per_case": manifest["resources"][
                    "wall_seconds_per_case"
                ],
                "peak_rss_bytes_per_case": manifest["resources"][
                    "peak_rss_bytes_per_case"
                ],
                "wall_seconds_campaign": manifest["resources"][
                    "wall_seconds_campaign"
                ],
                "worker_count": manifest["resources"]["worker_count"],
                "blas_threads": manifest["resources"]["blas_threads"],
                "final_diagnostics_json": json.dumps(
                    outcome.get("orders", [{}])[-1].get("diagnostics", {})
                    if outcome.get("orders")
                    else {},
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                    allow_nan=False,
                ),
                "converged": bool(outcome.get("converged", False)),
                "eligible": eligible,
                "failure_stage": record.get("failure_stage"),
                "failure_reason": record.get("failure_reason"),
            }
        )

    by_id = {record["case_id"]: record for record in ordered}
    for audit in ordered:
        parameters = audit["parameters"]
        if parameters["evidence_role"] != "rotational_audit":
            continue
        twin = by_id.get(parameters["twin_case_id"])
        rotation_error = None if twin is None else _rotation_error(audit, twin)
        derived_rows.append(
            {
                **_base(audit, manifest),
                "classification": manifest["classification"],
                "evidence_role": "rotational_audit",
                "eligible": _eligible(audit) and twin is not None and _eligible(twin),
                "include_in_scientific_tables": False,
                "metric": "rotational_covariance_error",
                "value": rotation_error,
                "unit": "1",
                "applicable": rotation_error is not None,
                "reason": (
                    "g1_budget_1e-12"
                    if rotation_error is not None
                    else "audit_or_twin_ineligible"
                ),
            }
        )

    artifacts = {
        "data_raw.csv": _csv_bytes(_RAW_FIELDS, raw_rows),
        "data_derived.csv": _csv_bytes(_DERIVED_FIELDS, derived_rows),
        "data_plot.csv": _csv_bytes(_PLOT_FIELDS, plot_rows),
        "failures.csv": _csv_bytes(_FAILURE_FIELDS, failure_rows),
        "performance.csv": _csv_bytes(_PERFORMANCE_FIELDS, performance_rows),
    }
    return artifacts, g1


def artifact_sha256(artifacts: Mapping[str, bytes]) -> dict[str, str]:
    """Return deterministic hashes keyed by frozen repository path."""

    if set(artifacts) != set(ARTIFACT_PATHS):
        raise CampaignArtifactError("campaign artifact set is incomplete")
    return {
        ARTIFACT_PATHS[name]: sha256(artifacts[name]).hexdigest()
        for name in sorted(artifacts)
    }


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def publish_campaign_artifacts(
    root: str | Path, artifacts: Mapping[str, bytes]
) -> dict[str, str]:
    """Publish the frozen artifact set atomically per file, without overwrite."""

    repository = Path(root)
    if set(artifacts) != set(ARTIFACT_PATHS):
        raise CampaignArtifactError("campaign artifact set is incomplete")
    destinations = {
        name: repository / relative for name, relative in ARTIFACT_PATHS.items()
    }
    existing = [str(path) for path in destinations.values() if path.exists()]
    if existing:
        raise FileExistsError(
            "campaign outputs exist; overwrite and second publication are forbidden: "
            + ", ".join(existing)
        )
    for name in sorted(destinations):
        _atomic_write(destinations[name], artifacts[name])
    return artifact_sha256(artifacts)
