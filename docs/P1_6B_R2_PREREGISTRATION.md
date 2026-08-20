# P1.6B-R2 response-blind replacement preregistration

Status: **frozen before replacement execution**.

P1.6B-R2 replaces only the invalid infrastructure execution documented in
`docs/P1_6B_INCIDENT.md`. It is not a scientific retuning. No force or
transient response from the invalid run is recoverable or used here.

## Frozen identity and namespaces

- campaign ID: `p1_dimer_confirmatory_r2`;
- manifest: `campaigns/p1/campaign_manifest_r2.yaml`;
- checkpoint namespace: `campaigns/p1/.p1_6b_r2_checkpoint/`;
- raw: `campaigns/p1/p1_6b_r2/data_raw.csv`;
- derived: `campaigns/p1/p1_6b_r2/data_derived.csv`;
- plot: `campaigns/p1/p1_6b_r2/data_plot.csv`;
- failures: `campaigns/p1/p1_6b_r2/failures.csv`;
- performance: `campaigns/p1/p1_6b_r2/performance.csv`;
- response-blind lock:
  `a041e07ae93e9a858bad809427039bf593641ad1f9e341ed89b9d91f648f297d`.

The P1.6A lock
`3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe`,
P1.4 historical lock
`9d360de6e61d901cff3f84c477f367773251103db12386dbb8156bd1ec2addca`,
and P1.5 lock
`d8f56ce20f6f0821d84fd6f36e1f76c855f63f55d809ba9a7201ba52097a43bf`
remain historical and unchanged.

## Scientific invariants preserved exactly

The R2 manifest is identical to the P1.6A manifest for schema,
classification, status, physical constants, geometry, numerical policy,
resources and the complete ordered case list. Therefore it retains:

- the same 102 case IDs in order 1--102;
- 96 primaries and six rotational audits;
- the same materials, contrasts, separations, orientations and audit twins;
- Model E at `L=2..21`, no stop before 5, tolerance `1e-5`, two consecutive
  confirmations and all four applicable channels;
- one worker and one BLAS thread;
- 1800 s and 4 GiB per case, with 64800 s total;
- the frozen metrics and G1 gate, including the `1e-12` identity/rotation
  budget.

Only campaign identity, provenance and checkpoint/output namespaces differ.

## Infrastructure gates

The JSON boundary converts only NumPy boolean, integer and real scalars through
`.item()`; complex and unknown types are rejected. Diagnostics must be fully
JSON-native before checkpointing. Unexpected serialization, contract or
infrastructure failures stop the campaign immediately, persist
`INVALID_P1.6B_R2_INFRASTRUCTURE` and propagate a nonzero process status.
They cannot be represented as scientific inelegibility, `unconfirmed_at_21`,
timeout or memory exhaustion.

Artifact publication requires a closed ledger, at least one completed case and
a non-null G1 decision persisted before publication. The read-only status gate
checks service state, decision, ledger counts and SHA-256 of all five outputs.
An entirely interrupted campaign cannot be
`READY_FOR_POSTPROCESSING`; explicit infrastructure invalidity is
`FAILED_INVALID`.

There is no retry, retuning, automatic restart or response-informed change.
