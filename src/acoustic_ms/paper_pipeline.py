"""Lightweight manifest validation for the confirmatory paper pipeline.

The campaign examples are JSON documents stored with a ``.yaml`` suffix.  JSON
is a strict subset of YAML, so the standard library is sufficient for P0 and
no runtime YAML or JSON-Schema dependency is introduced.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Mapping


class ManifestValidationError(ValueError):
    """Raised when a manifest does not satisfy its frozen P0 schema."""


_SCHEMA_BY_KIND = {
    "campaign": "campaign_manifest.schema.json",
    "figure": "figure_manifest.schema.json",
}


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
    """Validate one mapping against the supported P0 JSON-Schema subset."""

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
        numerical = manifest["numerical"]
        if numerical["lmax_min"] > numerical["lmax_max"]:
            errors.append("$.numerical: lmax_min must not exceed lmax_max")
        case_ids = [case["case_id"] for case in manifest["cases"]]
        if len(case_ids) != len(set(case_ids)):
            errors.append("$.cases: case_id values must be unique")
        if (
            manifest["status"] != "planned"
            and manifest["provenance"]["manifest_sha256"] == "TBD"
        ):
            errors.append(
                "$.provenance.manifest_sha256: TBD is allowed only for planned status"
            )
    if not errors and schema_id.endswith("/figure_manifest.schema.json"):
        table_ids = {table["table_id"] for table in manifest["source_tables"]}
        for index, panel in enumerate(manifest["panels"]):
            if panel["source_table"] not in table_ids:
                errors.append(
                    f"$.panels[{index}].source_table: unknown source table"
                )
    if errors:
        raise ManifestValidationError("manifest validation failed:\n- " + "\n- ".join(errors))


def validate_manifest_file(
    manifest_path: str | Path,
    *,
    kind: str,
    schema_directory: str | Path | None = None,
) -> dict[str, Any]:
    """Load and validate a campaign or figure manifest, returning its mapping."""

    if kind not in _SCHEMA_BY_KIND:
        raise ValueError(f"kind must be one of {sorted(_SCHEMA_BY_KIND)}")
    if schema_directory is None:
        schema_directory = Path(__file__).resolve().parents[2] / "campaigns" / "schemas"
    schema = load_json_yaml(Path(schema_directory) / _SCHEMA_BY_KIND[kind])
    manifest = load_json_yaml(manifest_path)
    validate_manifest(manifest, schema)
    return manifest
