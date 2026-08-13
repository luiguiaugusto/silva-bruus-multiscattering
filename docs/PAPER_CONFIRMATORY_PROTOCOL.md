# Paper confirmatory protocol

## DRAFT — NOT FROZEN

P0 records the decision surface; it does not preregister P1, open a new
sample, or generate a scientific response. Every item marked **TBD** must be
resolved in a versioned protocol before the associated campaign starts.

P1.1 status (2026-08-13): `docs/P1_1_DIMER_DECISION_RECORD.md` freezes the
dimer decisions and 102 ordered confirmatory IDs. The compatible schema
`1.1.0` stores physical sweep values per case. The confirmatory manifest and
the separate one-case `development` pilot manifest are valid, planned and
fully disabled; final hashes and enablement remain deferred to P1.4. This does
not authorize a solve.

## Paper question and intended answer

Question: when and why does the pairwise nodal Silva–Bruus interaction cease
to approximate the complete interaction force of a cluster of identical
fluid spheres?

Intended answer to test, not a conclusion: the discrepancy can be separated
exactly into an intrapair complete-dimer correction and a collective
many-body residual; connected contributions organize that residual, while
\(\Lambda_{\max}\) may provide a conservative quantitative validity
criterion within a declared physical domain.

The frozen narrative order is:

```text
Silva–Bruus
→ intrapair rescattering correction
→ connected many-body contributions
→ collective descriptor Lambda_max
→ quantitative validity criterion
→ independent confirmatory validation
```

## Physical domain and limitations

The established implementation uses \(e^{-i\omega t}\), explicit SI
\(E_0\), identical non-overlapping spheres, Cartesian centers in a pressure
nodal plane and an ideal unbounded fluid
(`docs/CONVENTIONS.md`; validation in
`src/acoustic_ms/solver.py::_validate_positions`). Model E represents
lossless fluid spheres with exact isolated-sphere Mie coefficients and a
complete multipolar interaction force
(`src/acoustic_ms/mie_scattering.py::mie_scattering_coefficients_from_contrasts`;
`src/acoustic_ms/model_e.py::solve_model_e_nodal`).

The paper must not generalize to antinodal excitation, arbitrary off-plane
centers, viscosity, streaming, walls, elasticity, absorption, nonspherical
particles, contact, torque, trajectories, or dynamics. P1's \(ka\), material,
separation and orientation grid is frozen in the P1.1 record; corresponding
ranges for P2–P6 remain **TBD**.

## Scientific claims to test

| ID | Proposed claim | Required evidence |
|---|---|---|
| C1 | A/SB is the controlled pairwise baseline in its stated nodal, non-overlap domain. | Analytical convention audit plus P1 limiting and symmetry checks. |
| C2 | \(B_E-A\) isolates the complete isolated-dimer correction under a common Model-E definition. | P1 canonical dimers, exact vector identities and converged E channels. |
| C3 | \(E-A=(B_E-A)+(E-B_E)\) holds componentwise, and \(E-B_E\) is the collective residual. | P2/P4 signed-vector reconstruction at every eligible case. |
| C4 | The collective residual admits connected \(\Phi_E^{(n)}\) terms through the preregistered order, without interpreting RMS amplitudes as additive fractions. | P2 subset inclusion–exclusion and reconstruction tests; target order through five is **TBD pending feasibility**. |
| C5 | \(\Lambda_{\max}\) captures the dominant variation of \(\varepsilon_A^E\) better than the declared comparator set in the development domain. | P3 leakage-safe grouped comparison, with candidates frozen before response use. |
| C6 | A frozen conservative \(\Lambda_{\max}\) rule controls false-safe classifications at declared tolerances with non-vacuous coverage. | P4 internal lock and P5 independent confirmation; tolerances and minima are **TBD**. |
| C7 | The rule transfers across preregistered sizes/geometries only within explicit boundaries; failures and extrapolations remain visible. | P5 external sample plus P6 boundary/sensitivity table and exclusions ledger. |

C1–C7 are hypotheses until their gates pass. Prior T01–T14.1 results may
motivate design and appear as `exploratory`, `development`, or
`legacy_validation`; none is relabeled `confirmatory_new`.

## Editorial model layer

- A/SB calls the existing pair force
  (`src/acoustic_ms/silva_bruus.py::nodal_pair_force_on_probe`).
- Historical B keeps its old Rayleigh meaning
  (`src/acoustic_ms/comparison.py::compare_nodal_force_models`).
- \(B_E\) is reserved for the future sum of independently converged isolated
  Model-E dimer interaction forces. It does not yet exist.
- C is global dipolar Model C
  (`src/acoustic_ms/force.py::solve_rayleigh_nodal_interaction_forces`).
- D is the multipolar leading-Rayleigh bridge
  (`src/acoustic_ms/model_d.py::solve_multipolar_nodal_interaction_forces`).
- E is the complete reference
  (`src/acoustic_ms/model_e.py::solve_model_e_nodal`), with
  `interaction_forces_xyz` as the paper response.

## Campaigns P1–P6

| Campaign | Role and hypothesis | Inputs/outputs | Evidence and stop condition |
|---|---|---|---|
| P1 — canonical dimer benchmark | Build and validate \(B_E\); test C1–C2 before any cluster campaign. | Frozen input: 96 primary dimers plus six rotational audits; separate development pilot excluded from scientific tables. Output raw E channels, A, B_E, convergence, identity and plot-ready dimer table. | Independent limiting/symmetry checks, channel convergence and deterministic artifacts. Stop on G1 pass or explicit redesign; do not begin P2 on silent partial success. |
| P2 — connected complete-force hierarchy | Build common-protocol Model-E subset expansion; test C3–C4. | Particle counts, families and maximum connected order: **TBD**. Output subset ledger, \(\Phi_E^{(n)}\), reconstruction and cost table. | Every subset retained; common convergence policy; exact signed reconstruction. Stop if cost or unresolved subsets make the declared order infeasible. |
| P3 — criterion development and freeze | Compare preregistered predictor candidates and freeze one law; test development part of C5. | Development groups, candidate transforms and sample size: **TBD**. Output OOF predictions, coefficients, margins and frozen manifest. | Leakage-safe grouped validation and identifiability checks. Stop after one prespecified hierarchy; no post-hoc search. |
| P4 — locked internal verification | Verify the frozen law and A/B_E/E decomposition on a response-blind locked set; test C3 and internal C6. | Locked cases/tolerances/coverage minima: **TBD**. Output eligibility ledger, threshold audit and gate. | Sufficiency precedes science. Any refit ends confirmation and returns to a separately named development cycle. |
| P5 — independent confirmatory validation | Test C5–C7 on a new sample not used by P1–P4. | Sizes, families, contrast/frequency ranges and levels: **TBD**. Output raw E data, frozen predictions, safety audit and boundary flags. | Frozen hashes before first response, adequate eligible coverage and preregistered false-safe/antivacuity gate. |
| P6 — synthesis and reproducible package | Assemble claims, limitations, final tables/figures and archive; no model selection. | Only frozen P1–P5 outputs. Output manuscript-ready data, figure manifests, exclusions ledger and reproducibility report. | Every panel rebuilds from tables; claims map to gates; all `TBD` resolved or explicitly excluded. |

## Metrics

Primary:

1. \(\varepsilon_A^E=\mathcal R(F^A-F^E_{\rm int})/
   \mathcal R(F^E_{\rm int})\), with applicability from
   `src/acoustic_ms/model_e_comparison.py::normalized_rms_error_xyz`.
2. \(\varepsilon_{B_E}^E=\mathcal R(F^{B_E}-F^E_{\rm int})/
   \mathcal R(F^E_{\rm int})\), to be implemented/tested in P1.
3. Componentwise residuals of
   \(E-A=(B_E-A)+(E-B_E)\) and connected reconstructions.
4. Strict conservative false-safe count and non-vacuous predicted-safe
   coverage at tolerances **TBD**.

Secondary:

- RMS amplitudes \(Y_n^E=\mathcal R(\Phi_E^{(n)})/
  \mathcal R(F^E_{\rm int})\), never additive scalar fractions.
- Multiplicative log RMSE/MAE, factor coverage, Spearman, bias and worst
  underprediction.
- Channel-specific truncation changes and numerical diagnostics.
- C/D bridge errors and signed mechanism cosines/projections.
- Cost, system dimension and peak memory.

## Predictor candidates

The provisional list is \(\eta\), \(\Lambda_{\max}\), \(\rho_1\), and
explicitly reference-derived bridge diagnostics. Implementations exist for
the first three
(`src/acoustic_ms/scaling.py::coupling_eta`,
`maximum_geometric_coupling`;
`src/acoustic_ms/transferability.py::spectral_radius_l1`).

Before P3, decide **TBD**:

- whether M1 in \(\Lambda_{\max}\) is the sole confirmatory candidate or is
  selected again in the new development protocol;
- allowed transforms/interactions and candidate hierarchy;
- grouped folds and minimum group count;
- whether historical coefficients are comparator-only;
- safety-margin estimator and exact anti-vacuity minima.

No response-driven feature addition, outlier removal, clipping or residual
floor is allowed after the P3 candidate set is frozen.

## Convergence, failures, missingness, and exclusions

Established E convergence requires two consecutive applicable relative
changes no larger than \(10^{-5}\), separately for total, interaction,
external–scattered and scattered–scattered channels
(`docs/CONVENTIONS.md`; channels from
`src/acoustic_ms/model_e.py::solve_model_e_nodal`). P1 evaluates every integer
order 2--21, never stops before 5 and requires the rule in every applicable
channel. Policies for P2–P5 remain **TBD**.

For every attempted case:

- retain one row even if assembly, resource precheck, solve, numerical gate,
  convergence or metric applicability fails;
- store `eligible=false`, controlled `failure_stage` and
  `ineligibility_reason`;
- never label `unconfirmed_at_limit` as divergent;
- never impute a missing force or denominator;
- calculate scientific metrics only on the preregistered eligible set;
- report attempted/eligible/excluded counts and IDs by reason;
- preserve raw rows; corrections require a new derived version.

## Separation of development and confirmation

| Class | May select predictors/fits/thresholds? | May support a new confirmatory claim? |
|---|---:|---:|
| `exploratory` | Yes, descriptively | No |
| `development` | Yes, with leakage controls | No |
| `legacy_validation` | No new selection; contextual evidence only | No under P1–P6 |
| `confirmatory_new` | No refit or post-response selection | Yes, only through its frozen gate |

If a confirmatory gate prompts recalibration, that campaign closes as failed
or inconclusive. Recalibration occurs in a new development version and
requires a new independent confirmatory sample.

## Gates G1–G6

- **G1 / P1:** B_E API, validation tests, dimer convergence, exact identity,
  schemas and deterministic artifacts pass.
- **G2 / P2:** subset identities pass; feasibility/cost support the frozen
  maximum connected order.
- **G3 / P3:** candidate hierarchy, grouped validation, coefficients, margin,
  tolerances, samples and hashes are committed before P4.
- **G4 / P4:** sufficiency precedes science; frozen rule meets internal
  error/safety/coverage criteria.
- **G5 / P5:** phase-A hashes remain fixed; independent sufficiency and all
  literal scientific criteria pass.
- **G6 / P6:** C1–C7 trace to evidence; each panel rebuilds from plot-ready
  data and manifest; exclusions/limitations are complete.

Exact numerical thresholds are **TBD** unless already physical identities.

## Decisions required before freeze

Luigui and ChatGPT Work must decide:

1. P1.1 is resolved: grid, orientations, materials, \(ka\), ordered IDs and
   convergence policy are frozen but disabled until P1.4.
2. P1 uses independently confirmed pair orders; the common-order sensitivity
   audit is assigned to P1.3.
3. P2 counts/families/subset reuse and feasibility of \(\Phi_E^{(5)}\).
4. Which prior data may inform P3 and exact group/split construction.
5. Predictor hierarchy, fit, safety factor, tolerances and coverage minima.
6. P4/P5 independent sample boundaries and extrapolation flags.
7. Numerical diagnostic limits, resource ceilings and stop conditions.
8. Journal style requirements replacing `diagnostic-v1` at P6.
9. Archival target, license and environment/container record.

Until these are resolved and the status changes in a dedicated commit, this
protocol remains **DRAFT — NOT FROZEN**.
