# Paper execution plan

Status: **P1.1 decisions frozen; P1 disabled pending P1.2--P1.4**. No P1
solve is authorized by this file.

P0 was merged through PR #1 at merge commit
`926e639fe2d327eacd09a2542208500891399687`. The P1.1 decision surface is
recorded in `docs/P1_1_DIMER_DECISION_RECORD.md`. The `1.1.0` confirmatory
manifest freezes 102 ordered IDs, while a separate `development` manifest
freezes the one-case resource pilot. Every case remains disabled and both
hashes remain `TBD` until P1.4.

## Dependency chain

```text
P0 contracts
  → P1 canonical dimer / B_E
    → P2 connected Model-E hierarchy
      → P3 criterion development and freeze
        → P4 locked internal verification
          → P5 independent confirmation
            → P6 synthesis and reproducible release
```

The first task after P0 is **P1: canonical dimer benchmark**.

## Cost anchors from existing measurements

Only measured historical numbers are used:

- T14 recorded 2963.8310701187 s over 162 evaluated orders, or 18.30
  s/order arithmetically, with maximum dimension 1848 and peak 1,032,216 KiB
  (`results/data/t14_performance.csv`;
  `scripts/run_t14_scale_out.py::_evaluate_case`).
- T14.1 recorded 23209.431251620874 s over 160 orders, or 145.06 s/order,
  with maximum dimension 6930 and peak 7,393,840 KiB
  (`results/data/t14_1_performance.csv`;
  `scripts/run_t14_1_large_n.py::_resource_estimate`).
- The complete baseline suite took 87.93 s for 464 tests in the P0
  environment.

There is no existing per-order dimer timing table. Therefore P1 cost is not
given a fabricated point estimate. For planning only, \(n_c n_L\times18.30\)
seconds is a measured large-cluster arithmetic reference, not a dimer bound.
P1 freezes 102 confirmatory cases and at most 20 evaluated orders per case,
or 2,040 case-order evaluations. The corresponding 10.37 h arithmetic
comparison is not a dimer estimate. The separate development pilot must record
assembly/solve/postprocess time and peak memory before the confirmatory grid is
enabled.

## P1 — canonical dimer benchmark

Hypothesis: a sum of independently converged complete Model-E dimers provides
a well-defined \(B_E\), and \(B_E-A\) isolates intrapair rescattering while
preserving signs, symmetries and the E interaction-force convention.

Parameters/sampling: 96 primary cases are
\(2\;ka\times6\) material cases \(\times8\) separations at \(\theta=0\);
six \(\theta=\pi/4\) covariance audits give 102 total. Orders are 2--21 with
minimum stop 5 and two consecutive passes in every applicable channel.

Small tasks:

1. **P1.1 decision record** — complete. Scientific choices, 102 IDs/order,
   schema 1.1, failure policy and provisional resource limits are frozen while
   every case remains disabled.
2. **P1.2 B_E API** — add an importable complete-dimer-sum routine under
   `src/acoustic_ms`; preserve historical B.
3. **P1.3 formula/validation tests** — isolated dimer identity, permutation,
   action–reaction where physically applicable, rotation/reflection, limiting
   behavior, non-overlap and invalid inputs; include the deferred common-order
   sensitivity audit for independently converged pairs.
4. **P1.4 blind manifest** — generate final hashes, audit the frozen IDs and
   perform any authorized disabled-to-enabled transition in a dedicated
   response-blind commit.
5. **P1.5 timed pilot** — execute only the separate `development` rigid case
   at \(ka=0.1,d/a=2.1,\theta=0\); its force cannot enter P1.6 scientific
   tables.
6. **P1.6 campaign and analysis** — execute once, retain all cases, run
   deterministic analysis twice and evaluate G1.

Inputs: P0 schemas/contract, existing A
(`src/acoustic_ms/silva_bruus.py::nodal_pair_force_on_probe`) and E
(`src/acoustic_ms/model_e.py::solve_model_e_nodal`).

Outputs: `data_raw.csv`, `data_derived.csv`, `data_plot.csv`,
`fit_parameters.csv` if a fit is preregistered, failure ledger, performance
table, manifest, figure manifest, canonical dimer figure and benchmark table.

Diagnostics/tests: all E channels, balanced solve/closures, componentwise
identity, force applicability, deterministic hashes, schemas and full suite.

Acceptance/stop: G1 passes only with a validated B_E definition and every
preregistered numerical/scientific item. Otherwise stop before P2 and report
`NO_GO_P2` or `INCONCLUSIVE_P1`. Do not tune parameters after response.

Feeds: paper theory/methods; dimer benchmark figure; model-definition table.

## P2 — connected complete-force hierarchy

Hypothesis: \(E-B_E\) is reconstructed by signed connected Model-E terms
through a feasible frozen order, with exact inclusion–exclusion closure.

Parameters/sampling: particle counts (candidate maximum five), geometry
families, cases, common convergence policy and maximum body order are **TBD**.

Small tasks:

1. derive/test generic subset bookkeeping independent of the acoustic solver;
2. define common-order versus independently confirmed-subset policy;
3. add Model-E subset solve/cache and complete failure ledger;
4. benchmark one preregistered small cluster;
5. freeze P2 manifest and resource stop rule;
6. execute/analyze once; evaluate G2.

Inputs: P1 B_E, E and P0 contract.

Outputs: raw subset table, subset-to-parent map, \(\Phi_E^{(n)}\) component
table, reconstruction residuals, cost table, connected-hierarchy figure/table.

Diagnostics/tests: Möbius/inclusion–exclusion identities on synthetic vectors;
pair reduction to B_E; permutation invariance; every subset solve's
convergence/quality; no RMS additivity claim.

Cost: **TBD after P1 timing**. For a five-particle case, evaluating every
subset of size at least two entails 26 subset problems before cache reuse;
the frozen resource precheck must multiply measured P1 timings rather than
extrapolate from T14.1.

Acceptance/stop: G2 requires reconstruction and the declared-order
feasibility. Reduce order only before response and in a new protocol version;
otherwise `INCONCLUSIVE_P2_RESOURCE_LIMIT`.

Feeds: connected-decomposition figure, hierarchy table, C3–C4.

## P3 — criterion development and freeze

Hypothesis: a prespecified collective descriptor yields a leakage-safe
quantitative description of \(\varepsilon_A^E\) across development groups.

Parameters/sampling: admissible P1/P2 development cases, groups, predictors,
fit form, weighting, thresholds, margin and bootstrap policy are **TBD**.

Small tasks:

1. freeze development/holdout boundaries before fitting;
2. implement canonical adapters to derived long data;
3. preregister candidate hierarchy and identifiability checks;
4. generate OOF predictions and sensitivity/influence diagnostics;
5. select/fail according to literal criteria;
6. commit frozen coefficients, safety factor, thresholds and P4/P5 IDs.

Inputs: eligible P1/P2 derived data and declared legacy context.

Outputs: OOF predictions, fold fits, `fit_parameters.csv`, safety audits,
influence table, frozen-prediction manifest and criterion figure.

Diagnostics/tests: held-out-response leakage mutation, group completeness,
positive-domain validation, strict threshold equality, invalid/collinear
candidates, deterministic serialization.

Cost: analysis-only CPU/memory **TBD**; no Model-E solves are allowed in the
fit step. Use wall time recorded by the implementation.

Acceptance/stop: G3 passes only if one preregistered candidate satisfies all
criteria. No second search after failure; stop with `NO_GO_P4`.

Feeds: criterion-development figure/table; frozen P4/P5 protocol.

## P4 — locked internal verification

Hypothesis: the frozen P3 law and A/B_E/E decomposition satisfy prespecified
accuracy, safety and coverage on a response-blind internal set.

Parameters/sampling: locked cases, strata, target levels, tolerances, coverage
and order caps are **TBD**.

Small tasks:

1. generate response-blind coordinates/predictions and commit their hashes;
2. audit A/B_E and all phase-A identities without E responses;
3. resource precheck/pilot under frozen stop rules;
4. execute once and retain every attempted case;
5. run analysis twice without solves;
6. evaluate sufficiency first, then G4.

Inputs: P3 frozen law and P1/P2 validated APIs.

Outputs: blind/revealed manifests, raw forces, eligibility ledger, threshold
audit, internal-validation figure and gate.

Diagnostics/tests: phase-A identity, eligibility mutation tests, exact strict
inequalities, convergence and all numerical gates.

Cost: \(n_c n_L\) times measured P1/P2 unit costs; all terms **TBD** until
sample freeze.

Acceptance/stop: insufficient coverage gives `INCONCLUSIVE_P4`; a literal
scientific failure gives `NO_GO_P5`. No refit in P4.

Feeds: internal validation and methods audit.

## P5 — independent confirmatory validation

Hypothesis: the unchanged criterion has conservative, non-vacuous performance
on a genuinely new domain and exposes rather than hides boundaries.

Parameters/sampling: independent sizes/families, ranges, target levels,
replicates and extrapolation policy are **TBD**.

Small tasks:

1. freeze case generator, IDs, predictions, local descriptors and hashes;
2. verify no P5 response entered P1–P4 selection;
3. commit/push phase A before first solve;
4. execute once under resource controls;
5. independent stratified audit and two analysis-only regenerations;
6. sufficiency gate followed by G5.

Inputs: only frozen P3 rule plus response-blind P5 geometry.

Outputs: raw/derived/plot data, eligibility/failure ledger, threshold and
boundary audits, independent-validation figure/table and gate.

Diagnostics/tests: blind-hash integrity, no-refit assertions, strict safety,
per-stratum coverage, extrapolation flags and numerical convergence.

Cost: computed from frozen case/order counts and P1/P2 timings; **TBD**.

Acceptance/stop: preregistered sufficiency and every scientific criterion must
pass. P3 comparators cannot rescue/reject the primary gate.

Feeds: C5–C7, external-validation figure, limitations table.

## P6 — synthesis and reproducible release

Hypothesis: all accepted claims can be reconstructed from immutable tables
without new science or hidden exclusions.

Parameters/sampling: no new physical sample. Journal target, final dimensions,
archive and license are **TBD**.

Small tasks:

1. lock claim-to-evidence matrix and exclusions ledger;
2. create final plot-ready tables and figure manifests;
3. replace diagnostic style only after journal requirements freeze;
4. rebuild all panels in PDF/SVG/PNG from tables;
5. run fresh-environment reproduction and hash audit;
6. evaluate G6 and archive release.

Inputs: frozen P1–P5 outputs only.

Outputs: final figures/tables, manifest bundle, environment lock, provenance
report, archive identifier and manuscript evidence matrix.

Diagnostics/tests: schema validation, panel-source completeness, no solver
imports in plot generation, byte determinism where backend/environment permit,
links/hashes, full scientific suite.

Cost: analysis/rendering only; **TBD after the final panel count**.

Acceptance/stop: G6 requires every claim and panel to trace to frozen data,
all exclusions visible and all unresolved TBDs either decided or outside the
paper's claim scope.

## P0 recommendation

**GO_P1_WITH_CONDITIONS**:

- resolve P1 scientific TBDs and change the protocol from draft in a separate
  commit;
- implement/test B_E before any cluster solve;
- record a timed dimer pilot because no dimer cost measurement exists;
- preserve the mixed pre-P0 user worktree outside P1 staging;
- do not treat T13–T14.1 as `confirmatory_new`.

## Current P1.1 gate

`HOLD_P1 — P1.2, P1.3 AND P1.4 REQUIRED`. P1.2 through P1.6 remain
unopened. No `B_E` implementation, solver call, force, result artifact,
campaign enablement or final manifest hash is authorized by P1.1.
