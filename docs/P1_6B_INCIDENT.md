# P1.6B infrastructure incident

Status: **INVALID_P1.6B_INFRASTRUCTURE**.

This record is factual and preserves the failed execution as non-scientific
infrastructure evidence. The original worktree
`/tmp/silva-bruus-p1-6b-confirmatory-execution` and detached launch directory
`/tmp/p1_6b_detached.xdiF9B` remain byte-immutable and are not copied into any
valid scientific result path.

## Root cause

`evaluate_model_e_numerical_diagnostics` returned
`planar_symmetry_pass` as `numpy.bool_`. The dataclass `asdict` operation
and the shallow diagnostics copy in `_normalize_outcome` preserved that NumPy
scalar at
`$.orders[0].diagnostics.planar_symmetry_pass`. The strict
`_json_bytes/json.dumps` boundary then raised
`TypeError: Object of type bool is not JSON serializable`.

The solver returned before this serialization failure. The in-memory outcomes
were never checkpointed and must not be reconstructed from timing.

## Scope

- attempted: 102;
- completed: 0;
- interrupted: 102;
- all 102 checkpoints have `attempt_count=1`, `state=interrupted` and
  `outcome=null`;
- recoverable scientific data: **none**;
- classification: **INVALID_P1.6B_INFRASTRUCTURE**.

The runner treated every unexpected serialization exception as a local
interruption and continued through all cases. It then set `closed=true`
because no case remained `never_started`, despite `campaign_decision=null`.
The CLI published empty/inapplicable CSVs based only on `closed`. The detached
status script checked service success, valid JSON and `closed=true`, but did
not require a decision, completed cases or output hashes. That combination
caused the false `READY_FOR_POSTPROCESSING`.

## Preserved evidence hashes

| Evidence | SHA-256 |
|---|---|
| `campaign_ledger.json` | `86ef17f76199ae3f3fd73960c1e66b65e68c1aa191d812ca810c21c870a656b7` |
| `data_raw.csv` | `79b438c500dbf0fa6b2c8c63cca2e19490238dee086842486e09dbe5d2d93383` |
| `data_derived.csv` | `77cbe29ed3d2b915541321dbb609f1b0ba65839ed103e6821a36e2974422efcb` |
| `data_plot.csv` | `853a5438b201c596885f969c412e054eed7c3743b0f59b530c5aa45c4c7171f7` |
| `failures.csv` | `3240fcf936d74153fcedf9b92ef19e41ba3b71a7daee882ece3b820fd7f33a96` |
| `performance.csv` | `a9a40f644b34597d2c024724e9dee81a721d467b9fa5dbd17ee912903476f4c8` |

The complete ledger/checkpoint and detached-launch hashes were recorded in the
forensic audit preceding this repair. No original evidence is staged by R2.
