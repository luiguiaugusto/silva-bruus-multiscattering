"""Lightweight manifest validation for the confirmatory paper pipeline.

The campaign examples are JSON documents stored with a ``.yaml`` suffix.  JSON
is a strict subset of YAML, so the standard library is sufficient for P0 and
no runtime YAML or JSON-Schema dependency is introduced.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import math
from pathlib import Path
import re
import sys
from types import MappingProxyType
from typing import Any, Mapping


class ManifestValidationError(ValueError):
    """Raised when a manifest does not satisfy its versioned paper schema."""


_CAMPAIGN_SCHEMA_BY_VERSION = {
    "1.0.0": "campaign_manifest.schema.json",
    "1.1.0": "campaign_manifest.v1.1.schema.json",
}
_FIGURE_SCHEMA = "figure_manifest.schema.json"
_MANIFEST_KINDS = ("campaign", "figure")
_ZERO_MANIFEST_SHA256 = "0" * 64
_MANIFEST_SHA256_PATTERN = re.compile(
    rb'(?P<prefix>"manifest_sha256"[ \t\r\n]*:[ \t\r\n]*")'
    rb"(?P<value>TBD|[0-9a-f]{64})"
    rb'(?P<suffix>")'
)

P1_FROZEN_MANIFEST_SHA256 = MappingProxyType(
    {
        "p1_dimer_confirmatory": (
            "9d360de6e61d901cff3f84c477f367773251103db12386dbb8156bd1ec2addca"
        ),
        "p1_dimer_resource_pilot": (
            "d8f56ce20f6f0821d84fd6f36e1f76c855f63f55d809ba9a7201ba52097a43bf"
        ),
    }
)


def load_json_yaml(path: str | Path) -> dict[str, Any]:
    """Load a JSON-compatible YAML mapping using only the standard library."""

    manifest_path = Path(path)
    try:
        value = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestValidationError(
            f"cannot load JSON-compatible YAML from {manifest_path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise ManifestValidationError("manifest root must be an object")
    return value


def _normalized_manifest_bytes(
    exact_bytes: bytes,
) -> tuple[bytes, str, dict[str, Any]]:
    """Return exact bytes with only the embedded hash value normalized."""

    if not isinstance(exact_bytes, bytes):
        raise TypeError("exact_bytes must be bytes")
    try:
        decoded = exact_bytes.decode("utf-8")
        document = json.loads(decoded)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ManifestValidationError(
            f"manifest bytes must be exact UTF-8 JSON-compatible YAML: {exc}"
        ) from exc
    if not isinstance(document, dict):
        raise ManifestValidationError("manifest root must be an object")
    try:
        stored = document["provenance"]["manifest_sha256"]
    except (KeyError, TypeError) as exc:
        raise ManifestValidationError(
            "manifest provenance.manifest_sha256 is required for hashing"
        ) from exc
    if not isinstance(stored, str):
        raise ManifestValidationError(
            "manifest provenance.manifest_sha256 must be a string"
        )

    matches = list(_MANIFEST_SHA256_PATTERN.finditer(exact_bytes))
    if len(matches) != 1:
        raise ManifestValidationError(
            "manifest must contain exactly one JSON manifest_sha256 field"
        )
    match = matches[0]
    raw_value = match.group("value").decode("ascii")
    if raw_value != stored:
        raise ManifestValidationError(
            "raw manifest_sha256 bytes do not match the parsed provenance value"
        )
    normalized = (
        exact_bytes[: match.start("value")]
        + _ZERO_MANIFEST_SHA256.encode("ascii")
        + exact_bytes[match.end("value") :]
    )
    return normalized, stored, document


def manifest_sha256(exact_bytes: bytes) -> str:
    """Hash exact UTF-8 bytes after zeroing only the embedded hash value.

    The 64 ASCII zeroes avoid self-reference. Every other byte, including
    whitespace and the final newline, remains part of the SHA-256 input.
    """

    normalized, _, _ = _normalized_manifest_bytes(exact_bytes)
    return hashlib.sha256(normalized).hexdigest()


def manifest_file_sha256(path: str | Path) -> str:
    """Return :func:`manifest_sha256` for the exact bytes of one file."""

    manifest_path = Path(path)
    try:
        exact_bytes = manifest_path.read_bytes()
    except OSError as exc:
        raise ManifestValidationError(
            f"cannot read manifest bytes from {manifest_path}: {exc}"
        ) from exc
    return manifest_sha256(exact_bytes)


def verify_manifest_sha256(
    exact_bytes: bytes, *, expected_sha256: str | None = None
) -> str:
    """Verify the stored self-hash and an optional external frozen digest."""

    normalized, stored, _ = _normalized_manifest_bytes(exact_bytes)
    if stored == "TBD":
        raise ManifestValidationError(
            "provenance.manifest_sha256: TBD cannot authorize execution"
        )
    digest = hashlib.sha256(normalized).hexdigest()
    if not hmac.compare_digest(stored, digest):
        raise ManifestValidationError(
            "provenance.manifest_sha256: stored hash does not match exact bytes"
        )
    if expected_sha256 is not None and not hmac.compare_digest(
        digest, expected_sha256
    ):
        raise ManifestValidationError(
            "provenance.manifest_sha256: digest differs from frozen public lock"
        )
    return digest


def _matches_type(value: object, expected: str) -> bool:
    if expected == "object":
        return isinstance(value, Mapping)
    if expected == "array":
        return isinstance(value, list)
    if expected == "string":
        return isinstance(value, str)
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected == "null":
        return value is None
    raise ManifestValidationError(f"unsupported schema type {expected!r}")


def _validate(value: object, schema: Mapping[str, Any], path: str) -> list[str]:
    errors: list[str] = []
    expected = schema.get("type")
    if expected is not None:
        expected_types = [expected] if isinstance(expected, str) else list(expected)
        if not any(_matches_type(value, item) for item in expected_types):
            return [f"{path}: expected {' or '.join(expected_types)}"]

    if "const" in schema and value != schema["const"]:
        errors.append(f"{path}: expected constant {schema['const']!r}")
    if "enum" in schema and value not in schema["enum"]:
        errors.append(f"{path}: value {value!r} is not in the allowed set")

    if isinstance(value, str):
        if len(value) < schema.get("minLength", 0):
            errors.append(f"{path}: string is shorter than minLength")
        pattern = schema.get("pattern")
        if pattern is not None and re.fullmatch(pattern, value) is None:
            errors.append(f"{path}: string does not match {pattern!r}")

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            errors.append(f"{path}: value is below minimum {schema['minimum']}")
        if "exclusiveMinimum" in schema and value <= schema["exclusiveMinimum"]:
            errors.append(
                f"{path}: value is not above {schema['exclusiveMinimum']}"
            )

    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            errors.append(f"{path}: array has fewer than minItems entries")
        if schema.get("uniqueItems"):
            serialized = [json.dumps(item, sort_keys=True) for item in value]
            if len(serialized) != len(set(serialized)):
                errors.append(f"{path}: array entries must be unique")
        item_schema = schema.get("items")
        if item_schema is not None:
            for index, item in enumerate(value):
                errors.extend(_validate(item, item_schema, f"{path}[{index}]"))

    if isinstance(value, Mapping):
        properties = schema.get("properties", {})
        for key in schema.get("required", []):
            if key not in value:
                errors.append(f"{path}: missing required property {key!r}")
        if schema.get("additionalProperties") is False:
            for key in value:
                if key not in properties:
                    errors.append(f"{path}: unexpected property {key!r}")
        for key, item in value.items():
            if key in properties:
                errors.extend(_validate(item, properties[key], f"{path}.{key}"))
    return errors


def validate_manifest(
    manifest: Mapping[str, Any], schema: Mapping[str, Any]
) -> None:
    """Validate one mapping against the supported JSON-Schema subset."""

    errors = _validate(manifest, schema, "$")
    schema_id = str(schema.get("$id", ""))
    if not errors and schema_id.endswith("/campaign_manifest.schema.json"):
        physical = manifest["physical"]
        expected_ka = float(physical["radius_m"]) * float(physical["k_rad_m"])
        if not math.isclose(
            float(physical["ka"]),
            expected_ka,
            rel_tol=64.0 * sys.float_info.epsilon,
            abs_tol=0.0,
        ):
            errors.append("$.physical.ka: must equal radius_m * k_rad_m")
        errors.extend(_campaign_common_errors(manifest))
    if not errors and schema_id.endswith(
        "/campaign_manifest.v1.1.schema.json"
    ):
        errors.extend(_campaign_common_errors(manifest))
        errors.extend(_campaign_v1_1_errors(manifest))
    if not errors and schema_id.endswith("/figure_manifest.schema.json"):
        table_ids = {table["table_id"] for table in manifest["source_tables"]}
        for index, panel in enumerate(manifest["panels"]):
            if panel["source_table"] not in table_ids:
                errors.append(
                    f"$.panels[{index}].source_table: unknown source table"
                )
    if errors:
        raise ManifestValidationError("manifest validation failed:\n- " + "\n- ".join(errors))


def _tbd_paths(value: object, path: str = "$") -> list[str]:
    paths: list[str] = []
    if isinstance(value, str) and "TBD" in value:
        paths.append(path)
    elif isinstance(value, Mapping):
        for key, item in value.items():
            paths.extend(_tbd_paths(item, f"{path}.{key}"))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            paths.extend(_tbd_paths(item, f"{path}[{index}]"))
    return paths


def _campaign_common_errors(manifest: Mapping[str, Any]) -> list[str]:
    """Return semantic errors shared by campaign schema versions."""

    errors: list[str] = []
    numerical = manifest["numerical"]
    if numerical["lmax_min"] > numerical["lmax_max"]:
        errors.append("$.numerical: lmax_min must not exceed lmax_max")
    case_ids = [case["case_id"] for case in manifest["cases"]]
    if len(case_ids) != len(set(case_ids)):
        errors.append("$.cases: case_id values must be unique")
    if manifest["status"] != "planned":
        errors.extend(
            f"{path}: TBD is allowed only for planned status"
            for path in _tbd_paths(manifest)
        )
    return errors


def _campaign_v1_1_errors(manifest: Mapping[str, Any]) -> list[str]:
    """Return per-case physical and ordering errors introduced in version 1.1."""

    errors: list[str] = []
    cases = list(manifest["cases"])
    radius = float(manifest["physical"]["radius_m"])
    numerical = manifest["numerical"]
    minimum_stop = int(numerical["minimum_stop_lmax"])
    if not (
        int(numerical["lmax_min"])
        <= minimum_stop
        <= int(numerical["lmax_max"])
    ):
        errors.append(
            "$.numerical.minimum_stop_lmax: must lie within the lmax range"
        )

    orders = [case["case_order"] for case in cases]
    expected_orders = list(range(1, len(cases) + 1))
    if orders != expected_orders:
        errors.append(
            "$.cases: case_order values must be contiguous and match list order"
        )
    if manifest["status"] == "planned" and any(case["enabled"] for case in cases):
        errors.append("$.cases: planned campaigns must keep every case disabled")

    by_id = {case["case_id"]: case for case in cases}
    finite_fields = ("ka", "k_rad_m", "f0", "f1", "distance_ratio", "theta_rad")
    for index, case in enumerate(cases):
        path = f"$.cases[{index}]"
        parameters = case["parameters"]
        if any(
            not math.isfinite(float(parameters[field]))
            for field in finite_fields
        ):
            errors.append(f"{path}.parameters: physical values must be finite")
            continue

        expected_ka = radius * float(parameters["k_rad_m"])
        if not math.isclose(
            float(parameters["ka"]),
            expected_ka,
            rel_tol=64.0 * sys.float_info.epsilon,
            abs_tol=0.0,
        ):
            errors.append(
                f"{path}.parameters.ka: must equal physical.radius_m * k_rad_m"
            )

        f0 = float(parameters["f0"])
        f1 = float(parameters["f1"])
        material_model = parameters["material_model"]
        f0_applicable = parameters["f0_applicable"]
        if f0 >= 1.0:
            errors.append(f"{path}.parameters.f0: API value must be smaller than 1")
        if material_model == "fluid":
            if not f0_applicable:
                errors.append(
                    f"{path}.parameters.f0_applicable: fluid material requires true"
                )
            if not -2.0 < f1 < 1.0:
                errors.append(
                    f"{path}.parameters.f1: fluid material requires -2 < f1 < 1"
                )
        if material_model == "rigid":
            if f1 != 1.0:
                errors.append(
                    f"{path}.parameters.f1: rigid material requires f1=1"
                )
            if f0_applicable:
                errors.append(
                    f"{path}.parameters.f0_applicable: rigid material requires false"
                )
            if f0 != 0.0:
                errors.append(
                    f"{path}.parameters.f0: rigid API sentinel must be zero"
                )

        role = parameters["evidence_role"]
        twin_id = parameters["twin_case_id"]
        include = parameters["include_in_scientific_tables"]
        if role == "primary" and twin_id is not None:
            errors.append(f"{path}.parameters.twin_case_id: primary case requires null")
        if role in ("rotational_audit", "resource_pilot") and include:
            errors.append(
                f"{path}.parameters.include_in_scientific_tables: "
                f"{role} must be false"
            )
        if role == "resource_pilot" and twin_id is not None:
            errors.append(
                f"{path}.parameters.twin_case_id: resource pilot requires null"
            )
        if role == "rotational_audit":
            if twin_id is None:
                errors.append(
                    f"{path}.parameters.twin_case_id: rotational audit requires a twin"
                )
                continue
            twin = by_id.get(twin_id)
            if twin is None:
                errors.append(
                    f"{path}.parameters.twin_case_id: unknown case {twin_id!r}"
                )
                continue
            if twin_id == case["case_id"]:
                errors.append(f"{path}.parameters.twin_case_id: case cannot twin itself")
                continue
            twin_parameters = twin["parameters"]
            if twin_parameters["evidence_role"] != "primary":
                errors.append(
                    f"{path}.parameters.twin_case_id: twin must be a primary case"
                )
            matching_fields = (
                "ka", "k_rad_m", "material_id", "material_model", "f0",
                "f0_applicable", "f1", "distance_ratio",
            )
            if any(
                parameters[field] != twin_parameters[field]
                for field in matching_fields
            ):
                errors.append(
                    f"{path}.parameters.twin_case_id: twin physical parameters differ"
                )
            if not math.isclose(
                float(twin_parameters["theta_rad"]),
                0.0,
                rel_tol=0.0,
                abs_tol=0.0,
            ):
                errors.append(
                    f"{path}.parameters.twin_case_id: twin theta_rad must be zero"
                )
            if math.isclose(
                float(parameters["theta_rad"]),
                float(twin_parameters["theta_rad"]),
                rel_tol=0.0,
                abs_tol=0.0,
            ):
                errors.append(
                    f"{path}.parameters.theta_rad: audit must rotate its twin"
                )
    return errors


def validate_manifest_file(
    manifest_path: str | Path,
    *,
    kind: str,
    schema_directory: str | Path | None = None,
) -> dict[str, Any]:
    """Load and validate a campaign or figure manifest, returning its mapping."""

    manifest_path = Path(manifest_path)
    if kind not in _MANIFEST_KINDS:
        raise ValueError(f"kind must be one of {list(_MANIFEST_KINDS)}")
    if schema_directory is None:
        schema_directory = Path(__file__).resolve().parents[2] / "campaigns" / "schemas"
    manifest = load_json_yaml(manifest_path)
    if kind == "campaign":
        version = str(manifest.get("schema_version", ""))
        schema_name = _CAMPAIGN_SCHEMA_BY_VERSION.get(version)
        if schema_name is None:
            raise ManifestValidationError(
                f"unsupported campaign schema_version {version!r}"
            )
    else:
        schema_name = _FIGURE_SCHEMA
    schema = load_json_yaml(Path(schema_directory) / schema_name)
    validate_manifest(manifest, schema)
    if (
        kind == "campaign"
        and manifest["provenance"]["manifest_sha256"] != "TBD"
    ):
        try:
            exact_bytes = manifest_path.read_bytes()
        except OSError as exc:
            raise ManifestValidationError(
                f"cannot read manifest bytes from {manifest_path}: {exc}"
            ) from exc
        expected = P1_FROZEN_MANIFEST_SHA256.get(str(manifest["campaign_id"]))
        verify_manifest_sha256(exact_bytes, expected_sha256=expected)
    return manifest


def validate_executable_manifest_file(
    manifest_path: str | Path,
    *,
    expected_campaign_id: str,
    schema_directory: str | Path | None = None,
) -> dict[str, Any]:
    """Validate an immutable preregistration before any execution.

    Callers must name the expected campaign explicitly. A valid frozen hash is
    necessary but not sufficient: at least one case must be enabled.
    """

    manifest = validate_manifest_file(
        manifest_path,
        kind="campaign",
        schema_directory=schema_directory,
    )
    if manifest["campaign_id"] != expected_campaign_id:
        raise ManifestValidationError(
            "campaign_id does not match the explicitly authorized campaign"
        )
    if manifest["status"] != "preregistered":
        raise ManifestValidationError(
            "manifest status must be preregistered before execution"
        )
    if manifest["provenance"]["manifest_sha256"] == "TBD":
        raise ManifestValidationError("TBD hash cannot authorize execution")
    if not any(case["enabled"] for case in manifest["cases"]):
        raise ManifestValidationError("manifest has no enabled cases")
    return manifest
