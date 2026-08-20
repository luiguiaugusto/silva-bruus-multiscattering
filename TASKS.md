# Tasks

- [x] T01 -- infrastructure, contrasts, and nodal pairwise Silva--Bruus force.
- [x] T02 -- corrected two-body formula and reproduction benchmark.
- [x] T03 -- multipolar Rayleigh core and coupled Lmax=1 solver.

T03 provides the Rayleigh field solver; T04 adds the validated nodal two-particle interaction-force observable.
- [x] T04 -- Rayleigh Lmax=1 nodal interaction force for one and two particles.
- [x] T04.1 -- closure of T04 force-test coverage and documentation.
- [x] T05 -- first N=3 Model A/B/C comparison and irreducible multibody force.

- [x] T05.1 -- closure of T05 metrics, validation coverage, artifacts, and documentation.

- [x] T05.1a -- documentary closure of T05.1 audit records.
- [x] T05.1b -- final correction of T05.1 audit documentation.
- [x] T05.1c -- final cleanup of T05.1 audit documentation.
- [x] T06 -- N=4 connected body expansion through irreducible four-body forces.
- [x] T06.1 -- Scaling analysis of connected three- and four-body amplitudes.
- [x] T07 -- multipolar Model D and convergence of total and connected cluster forces.
- [x] T08 -- transferability through N=10 and collective-coupling validity diagnosis.
- [x] T09 -- analytical foundation and independent operator audit of rho_1.
- [x] T10 -- exact isolated-sphere Mie coefficients and Rayleigh-error audit.
- [x] T11 -- exact-Mie global Model E with the complete multipolar radiation force.
- [x] T11.1 -- numerical stabilization and high-precision audit of Model E.
- [x] T12 -- preregistered Model E sentinel audit of the frozen rho_1 criterion.
- [x] T12.1 -- Model E convergence extension and diagnosis of the frozen rho_1 failure.
- [x] T12.2 -- controlled LOGO recalibration of rho_1 against Model E.
- [x] T12.3 -- grouped mechanistic validation of the Lambda_max criterion.
- [x] T13 -- external validation of the frozen Lambda_max criterion on N=6 and N=10 clusters.
- [x] T14 -- scale-out validation of the frozen Lambda_max criterion on N=15 and N=28 clusters.
- [x] T14.1 -- frozen Lambda_max confirmation at N=45 and N=105.
- [x] T14/P0 -- methodological freeze and canonical paper-pipeline contracts;
  P1 remains unopened.
- [x] P1.1 -- canonical-dimer decisions, schema 1.1 and 102 ordered IDs
  frozen; confirmatory and development-pilot cases remain disabled, hashes
  remain deferred to P1.4, and P1.2--P1.6 remain unopened.
- [x] P1.2 -- importable independently converged complete-dimer \(B_E\) API,
  deterministic pair ledger, all-or-nothing eligibility and injected-solver
  unit tests; historical B/\(B_L\), disabled manifests and result trees remain
  unchanged, and P1.3--P1.6 remain unopened.

- [x] P1.3 -- physical validation and P1.3a hierarchical audit of `B_E`;
  the original G7 is preserved as strict xfail because it compares
  non-equivalent modal/force objects, while equivalent gates and the
  distinct-common-order fallback pass.
- [x] P1.4 -- exact-byte response-blind SHA-256 locks published for both P1
  manifests; all 102 confirmatory cases remain disabled and only the separate
  P1.5 resource-pilot case is enabled. No pilot, campaign or response ran.
- [x] P1.5 -- one response-blind rigid resource pilot executed from the
  separately pushed runner commit. Orders 2--21 ended
  `unconfirmed_at_21`, `eligible=false`; 494.133 s and 311857152-byte peak RSS
  pass the provisional limits, artifacts/derivations are deterministic and
  excluded from scientific tables. Decision: `GO_P1.6A_BLIND_FREEZE`.
- [x] P1.6A -- exact-byte confirmatory lock updated with all 102 cases enabled,
  96 scientific primaries plus six excluded rotational audits, frozen
  1-worker/1-thread/1800-s/4-GiB/64800-s limits, single-attempt checkpointed
  runner, pure deterministic analysis and pre-response G1. No confirmatory
  case ran and no confirmatory artifact was created. Decision:
  `GO_P1.6B_EXECUTE` pending audit.
- [x] P1.6A.1 -- response-blind runner amendment on the same draft PR: actual
  execution provenance and allowlisted environment are immutable on resume;
  abandoned floating wall reservations debit the global budget; raw/performance
  distinguish execution and manifest commits; `epsilon_a_e` is the main plot
  response while `epsilon_be_e` and absolute `be_minus_a_rms` remain explicit.
  Manifest/locks and all P1.5 bytes are unchanged; no confirmatory case ran.
- [x] P1.6A.2 -- response-blind CLI worktree amendment on draft PR #7: first
  execution requires a clean tree; resume admits only versionable files below
  the exact `.p1_6_checkpoint/` directory and refuses external/staged changes
  or preexisting confirmatory CSVs before any executor. P1.6B must stage the
  ledger, 102 existing case checkpoints, five CSVs and all hashes without
  cleaning or regenerating observed checkpoints. Frozen bytes remain unchanged.
- [x] P1.6B -- invalid infrastructure execution preserved: 102 attempted,
  zero completed, 102 interrupted and no scientific outcome; classified
  `INVALID_P1.6B_INFRASTRUCTURE`.
- [ ] P1.6B-R2 -- infrastructure-only replacement frozen under campaign ID
  `p1_dimer_confirmatory_r2`, new checkpoint/output namespaces and lock
  `a041e07ae93e9a858bad809427039bf593641ad1f9e341ed89b9d91f648f297d`.
  Launch is authorized only after focused/full tests, clean pre-solve commit,
  push and draft PR.
