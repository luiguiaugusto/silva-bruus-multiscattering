# Confirmatory campaign manifests

This directory is the entry point for paper campaigns P1–P6. P0 contains no
real campaign and the example case is disabled. Copy a template into a new
campaign directory only after its scientific `TBD` items have been decided and
the protocol has been preregistered.

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

Schemas are versioned at `1.0.0`. A schema change must increment
`schema_version`; legacy CSVs remain immutable and are adapted only in the
derived layer described by `docs/PAPER_DATA_CONTRACT.md`.
