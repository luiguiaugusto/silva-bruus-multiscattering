# Paper data and figure contract

Schema version: `1.0.0`

## Immutable levels

```text
raw solver output → derived analysis data → plot-ready data
```

1. **Raw solver output** records what was attempted and returned, including
   failures. Once a campaign is closed, raw bytes are immutable.
2. **Derived analysis data** is a deterministic, versioned transformation of
   named raw files. It may add aliases, metrics and eligibility, but never
   replace raw values.
3. **Plot-ready data** contains every value needed to draw a panel. Plot code
   must not solve physics, refit a model or infer an undocumented exclusion.

A correction creates a new file/schema version and provenance edge. It never
edits a closed upstream layer.

## Stable identifiers

`case_id` is lowercase ASCII matching
`^[a-z0-9][a-z0-9._-]*$`. It is assigned in the response-blind campaign
manifest and must not depend on a force, error, convergence result or plot
order. It never changes between raw, derived and plot-ready layers.

Recommended construction:

```text
<campaign>_<particle-count>_<family>_<design-level>_<replicate>
```

`campaign_id` identifies the campaign; `solve_id` identifies an attempt or
rerun; `analysis_id` identifies a transformation; `figure_id` and
`panel_id` identify graphical consumers. Compound uniqueness is:

- raw: (`campaign_id`, `case_id`, `solve_id`, `lmax`,
  `particle_index`, `model`, `force_channel`);
- derived metric: (`analysis_id`, `case_id`, `model`, `metric`);
- plot ready: (`figure_id`, `panel_id`, `series_id`, `point_order`).

Legacy files without `case_id` remain immutable. An adapter creates a
deterministic key in the derived layer and stores `legacy_source_path` plus
`legacy_row_number`.

## Minimum solve metadata

Every raw row or its lossless one-to-one case table must provide:

- identifiers: `schema_version`, `campaign_id`, `case_id`, `solve_id`;
- classification: `exploratory`, `development`,
  `legacy_validation`, or `confirmatory_new`;
- physical: `radius_m`, `k_rad_m`, `ka`,
  `energy_density_j_m3`, `f0`, `f1`,
  `temporal_convention=exp(-i omega t)`;
- geometry: `particle_count`, `particle_index`,
  `position_x_m`, `position_y_m`, `position_z_m`,
  `minimum_distance_m`, `distance_ratio`, `family`,
  `coordinate_sha256`;
- numerical: `model`, `force_channel`, `lmax`,
  `full_modes_per_particle`, `active_modes_per_particle`,
  `system_dimension`, `solver_name`, worker/thread counts;
- force: `force_x_n`, `force_y_n`, `force_z_n` and, when used,
  `force_x_over_a2e0`, `force_y_over_a2e0`,
  `force_z_over_a2e0`;
- convergence: per-channel successive/absolute change, applicability,
  confirmation order, confirmed flag and stop reason;
- quality: finite, mode-dimension, planar-symmetry and resource-precheck
  flags; balanced condition/backward error; incident/scattering closure and
  force-decomposition residual;
- eligibility: `attempted`, `eligible`, `failure_stage` and
  `ineligibility_reason`;
- provenance: `git_commit`, `campaign_manifest_path`,
  `campaign_manifest_sha256`, `schema_version` and `created_utc`.

`ka` and (`radius_m`, `k_rad_m`) must be consistent; sweeps may not treat
\(ka\), \(kd\), and \(d/a\) as independent (`docs/CONVENTIONS.md`). A future
runner must validate this before solving.

## Forces and normalizations

Physical forces are stored in newtons. Normalized values use only the explicit
scale \(a^2E_0\) and carry `_over_a2e0` in the name. A normalized value never
occupies a `_n` column.

Model E's paper reference is
`ModelENodalResult.interaction_forces_xyz`
(`src/acoustic_ms/model_e.py::ModelENodalResult`). Total, external,
external–scattered and scattered–scattered remain for diagnostics. A/B/C/D
planar vectors receive explicit zero z only in a documented derived step.

`model` uses `A`, `B_E`, `C`, `D`, `E`. Historical `B` and `B_L`
remain legacy values and must never be silently mapped to `B_E`.

## Convergence, quality, and exclusions

Missing, inapplicable, unconfirmed and failed are distinct:

- `*_applicable=false`: no resolved denominator; numeric field blank.
- `*_confirmed=false`: direct convergence rule did not pass.
- `finite=false`: numerical output was not finite.
- `eligible=false`: requires a controlled nonempty reason.
- A missing response is never imputed or represented by fabricated zero.

All attempted cases survive into the eligibility ledger. Scientific filtering
is recorded in `inclusion_rule`. Counts and IDs by exclusion reason accompany
every summary.

## CSV serialization

New CSVs use UTF-8, LF, a header, RFC-4180 quoting, lowercase
`true`/`false`, empty missing fields, and no thousands separators. Floats use
round-trip-safe decimal serialization (equivalent to Python `.17g` for
binary64). Display rounding occurs only in rendering. Units occur in names or
an explicit `unit` column.

Vectors use long rows or separate component columns. Legacy JSON-serialized
arrays are not the canonical new representation.

### `data_raw.csv`

One row per case/order/particle/model/channel:

```text
schema_version,campaign_id,case_id,solve_id,classification,created_utc,
git_commit,campaign_manifest_path,campaign_manifest_sha256,
radius_m,k_rad_m,ka,energy_density_j_m3,f0,f1,temporal_convention,
particle_count,particle_index,family,position_x_m,position_y_m,position_z_m,
minimum_distance_m,distance_ratio,coordinate_sha256,
model,force_channel,lmax,full_modes_per_particle,active_modes_per_particle,
system_dimension,solver_name,force_x_n,force_y_n,force_z_n,
force_x_over_a2e0,force_y_over_a2e0,force_z_over_a2e0,
successive_change,absolute_change,change_applicable,confirmation_lmax,
confirmed,finite,diagnostics_pass,eligible,failure_stage,
ineligibility_reason,stop_reason
```

Detailed numerical diagnostics may be additional columns but may not replace
the minimum quality fields.

### `data_derived.csv`

Long metric table:

```text
schema_version,analysis_id,created_utc,git_commit,source_path,source_sha256,
campaign_id,case_id,classification,eligible,inclusion_rule,
model,reference_model,metric,value,unit,applicable,reason
```

Signed force identities require component-level rows in addition to RMS.
`metric` names are versioned; examples: `epsilon_a_e`, `epsilon_be_e`,
`lambda_max`, `rho_l1` and `connected_y3_e`.

### `data_plot.csv`

One row per graphical mark:

```text
schema_version,figure_id,panel_id,series_id,point_order,case_id,
x_name,x_value,x_unit,y_name,y_value,y_unit,
xerr_low,xerr_high,yerr_low,yerr_high,
marker,color,linestyle,label,eligible,annotation
```

Every filter/transformation is already applied and repeated in the figure
manifest. A panel must rebuild from `data_plot.csv` plus
`figure_manifest.yaml` without reading an image or solver table.

### `fit_parameters.csv`

One row per fit/fold/version:

```text
schema_version,fit_id,model_name,evidence_role,training_split,held_out_group,
predictors,response,transform,weighting,point_count,group_count,
intercept,prefactor,exponent_lambda,exponent_rho,safety_factor,
rmse_log,mae_log,created_utc,git_commit,source_sha256,applicable,reason
```

Unused coefficients are blank, not zero. Frozen/descriptive fits use distinct
`evidence_role` values.

## `figure_manifest.yaml`

The schema is `campaigns/schemas/figure_manifest.schema.json`. Required:

- `figure_id`, title, UTC time and output stem;
- commit, campaign manifest and generator;
- `style_profile` and PDF/SVG/PNG formats;
- source path, immutable SHA-256 and data level;
- every panel's ID/label/title, table, x/y columns, filters and ordered
  transformations.

The example `campaigns/templates/figure_manifest.example.yaml` validates via
`src/acoustic_ms/paper_pipeline.py::validate_manifest_file`.

## Provenance and immutability

`created_utc` is RFC-3339 UTC (`YYYY-MM-DDTHH:MM:SSZ`). `git_commit` is the
40-character source commit. The manifest is hashed from exact bytes after
response-blind fields freeze; `TBD` is allowed only while `status=planned`
and cannot start a campaign.

Each layer records upstream path and SHA-256. Final archives contain:

1. manifest/schema;
2. raw bytes/hash;
3. analysis code/commit and derived bytes/hash;
4. plot-ready bytes/hash;
5. figure manifest, generator and exports;
6. exclusions/failure ledger and environment record.

## Legacy compatibility

No file under `results/` is migrated or regenerated by P0. Adapters must:

- read a declared immutable legacy hash;
- map names in a versioned function;
- preserve 2-D/3-D semantics and applicability flags;
- emit a new derived file with source path/hash;
- test row count, key uniqueness, units and representative formulas;
- never overwrite the source.

Thus all T01–T14.1 CSVs coexist with the new P1–P6 pipeline.
