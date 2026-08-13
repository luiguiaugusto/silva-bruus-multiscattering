# Confirmatory campaign manifests

This directory is the entry point for paper campaigns P1–P6. A planned
manifest remains non-executable: its cases must be disabled and its final
hash is deferred until preregistration.

The `.yaml` examples deliberately use JSON syntax. JSON is valid YAML and lets
the lightweight validator use Python's standard library without adding a YAML
or JSON-Schema runtime dependency.

Validate an example with:

```python
from acoustic_ms.paper_pipeline import validate_manifest_file

validate_manifest_file(
    "campaigns/templates/campaign_manifest.example.yaml",
    kind="campaign",
)
```

Campaign schema `1.0.0` remains available at
`campaigns/schemas/campaign_manifest.schema.json` for existing single-physical
manifests. Schema `1.1.0` is additive and selected automatically by
`validate_manifest_file`; it keeps radius, energy, geometry convention,
numerical policy and resource policy global while requiring `ka`, `k_rad_m`,
material, `f0`, `f0_applicable`, `f1`, `distance_ratio` and `theta_rad` for
every ordered case.

The examples cover both versions:

- `campaign_manifest.example.yaml`: unchanged `1.0.0`;
- `campaign_manifest.multicase.example.yaml`: disabled per-case `1.1.0`.

P1 uses two disabled `1.1.0` manifests. `campaign_manifest.yaml` contains the
102-case confirmatory design; `pilot_manifest.yaml` contains one separate
`development` resource pilot that is forbidden from entering P1.6 scientific
tables. A schema change must increment `schema_version`; legacy manifests and
CSVs remain immutable.
