# P1.1 canonical-dimer decision record

Status: **DECISIONS FROZEN — CASES DISABLED; P1.4 REQUIRED**

Decision date: 2026-08-13

Decision owners: Luigui and ChatGPT Work

This record applies the scientific review of draft PR #2. It freezes the P1
case IDs, their order and the policies needed by later implementation tasks.
It does not implement \(B_E\), authorize a solver call, generate a force,
execute the pilot/campaign or open P1.2--P1.6.

The final manifest hashes and every `enabled=true` transition are explicitly
deferred to P1.4. Until then, both P1 manifests have `status=planned`,
`manifest_sha256=TBD` and all cases disabled.

## Fixed scope

- Two identical, non-overlapping spheres in a pressure nodal plane of an ideal
  unbounded fluid, with \(e^{-i\omega t}\), explicit SI \(E_0\), and the
  Model-E interaction force as the paper response.
- Total, interaction, external--scattered and scattered--scattered channels
  remain recorded and must converge when applicable.
- Historical B and \(B_L\) keep their meanings. \(B_E\) remains unimplemented
  and reserved for a future sum of independently converged Model-E dimers.
- The campaign uses schema `1.1.0`; legacy manifest schema `1.0.0` remains
  supported without reinterpretation.

## Approved physical grid

The 96 primary cases are the exact Cartesian product

\[
\{0.05,0.1\}_{ka}
\times
\{\mathrm{M01},\ldots,\mathrm{M06}\}_{material}
\times
\{2.1,2.25,2.5,3,4,6,8,10\}_{d/a}
\times
\{0\}_{\theta}.
\]

The ordered material cases are:

| Order | `material_id` | Model | \(f_0\) | `f0_applicable` | \(f_1\) |
|---:|---|---|---:|:---:|---:|
| M01 | `fluid_f000_f1010` | fluid | 0 | true | 0.1 |
| M02 | `fluid_f000_f1040` | fluid | 0 | true | 0.4 |
| M03 | `fluid_f000_f1060` | fluid | 0 | true | 0.6 |
| M04 | `fluid_f020_f1060` | fluid | 0.2 | true | 0.6 |
| M05 | `fluid_f000_f1080` | fluid | 0 | true | 0.8 |
| M06 | `rigid_boundary` | rigid | API sentinel 0 | false | 1 |

For M06, `material_model=rigid` and \(f_1=1\) select the exact sound-hard
boundary. The stored \(f_0=0\) is only the finite sentinel required by the
current API before it dispatches to the rigid coefficient routine.
`f0_applicable=false` forbids interpreting that sentinel as a physical
monopole contrast.

## Frozen IDs and order

`campaigns/p1/campaign_manifest.yaml` contains exactly 102 ordered, unique and
disabled cases:

1. cases 1--96: \(ka\) in ascending order; within each \(ka\), materials M01
   through M06; within each material, \(d/a\) in the listed ascending order;
   all at \(\theta=0\);
2. cases 97--102: one \(\theta=\pi/4\) rotational-covariance audit for each
   material M01--M06, balanced over both \(ka\) values and the two separation
   extremes;
3. every audit records the exact corresponding \(\theta=0\) primary
   `twin_case_id`.

The audit order and physical mapping are frozen as follows:

| Case order | Material | \(ka\) | \(d/a\) |
|---:|---|---:|---:|
| 97 | M01 | 0.05 | 2.1 |
| 98 | M02 | 0.05 | 10 |
| 99 | M03 | 0.1 | 2.1 |
| 100 | M04 | 0.1 | 10 |
| 101 | M05 | 0.05 | 2.1 |
| 102 | M06 | 0.1 | 10 |

Primary IDs use
`p1_dimer_<ka-tag>_<material-id>_<distance-tag>_t000`. Audit IDs append
`_t045_audit`. `case_order` is the immutable integer sequence 1--102. The
manifest, not a generated force or plot order, is authoritative.

The 96 primary cases are eligible for future P1 scientific tables subject to
the frozen numerical gates. The six audit cases have
`include_in_scientific_tables=false` and test rotational covariance only.

## Approved numerical and failure policy

- Evaluate every integer \(L_{\max}=2,\ldots,21\).
- Do not stop before \(L_{\max}=5\).
- Require two consecutive applicable normalized RMS changes no larger than
  \(10^{-5}\) for every applicable force channel.
- A null channel follows the established scale-relative applicability rule;
  it is not assigned a fabricated zero change.
- Reuse the established quality gates: finite outputs, balanced condition
  number below 10, balanced backward error, incident closure, scattering
  closure and force-decomposition residual each below \(10^{-12}\), plus
  \(\max|F_z|\le128\epsilon_{mach}F_{max}\).
- Retain every attempt and failure. Never impute a response or retry by
  changing physics, tolerance or solver. Record controlled failure stage and
  reason; unconfirmed at the cap is not called divergent.
- Future \(B_E\) sums use independently confirmed pair orders and are eligible
  only if every pair is eligible. The common-order sensitivity audit is
  assigned to P1.3 and is not implemented here.

## Resources and separate pilot

The provisional resource policy is one worker, one BLAS thread, 4 GiB peak RSS
per case, 30 minutes per case and 12 hours for the campaign. These limits are
recorded as `limits_status=provisional`; the development pilot may support a
new pre-campaign resource decision but may not alter scientific responses.

`campaigns/p1/pilot_manifest.yaml` freezes one separate `development` case:

```text
p1_pilot_rigid_ka010_d0210_t000
material_model=rigid
ka=0.1
d/a=2.1
theta=0
```

The pilot is disabled and has no response. If later executed under its own
authorization, its force and derived metrics must never enter P1.6 scientific
tables. Its permitted outputs are resource/timing, serialization and
numerical-health evidence.

## Reproducibility

P1 analysis must be regenerated twice without importing or calling a solver.
The six frozen \(\pi/4\) cases are the rotational-covariance audits and must
be compared componentwise with their six \(\theta=0\) twins. No full campaign
duplication is approved.

## Cost record

The confirmatory manifest contains 102 cases. Orders 2--21 give at most 20
evaluated orders per case, hence at most 2,040 case-order evaluations; stopping
at order 5 gives a minimum of four evaluated orders per case, or 408. Using
the unrelated T14 arithmetic reference of 18.30 s/order gives 37,332 s
(10.37 h) only as a large-cluster comparison, not a dimer forecast or bound.

For a planar dimer, the balanced system dimension is \(L(L+1)\): 182 at
\(L=13\) and 462 at \(L=21\). One complex128 dense matrix is approximately
0.51 MiB and 3.26 MiB at those orders before translation tables,
factorizations and runtime overhead. Dense cubic solve work at \(L=21\) is
about 135 times the \(L=9\) work. Actual wall time and peak memory remain
unmeasured until the separate pilot is authorized.

## Applied decision matrix

| ID | Existing evidence | Admissible options considered | Recommendation / applied decision | Estimated cost impact |
|---|---|---|---|---|
| D01 | T10 validates through \(ka=0.1\); T11 has dimer precedent at 0.05 and 0.1. | One validated level, both precedent levels, or a broader grid only inside the validated interval. | `APPROVED`: \(ka=\{0.05,0.1\}\). | Doubles the primary grid relative to one \(ka\): 48 additional primary cases. |
| D02 | Model-E campaigns mostly use \(f_0=0\); the T10 oracle covers \(f_0=0.2\). | Legacy-only zero, one paired positive point, or a free factorial. | `APPROVED`: use \(f_0=0\) in four fluid cases, \(f_0=0.2\) in one fluid case, and a nonphysical API sentinel zero for rigid. | Already included in the six material levels; a free factorial was rejected because it would multiply the grid. |
| D03 | T10 covers \(f_1=0.1,0.4,0.8,1\) and the finite point 0.6; exactly 1 selects the rigid routine. | Historical positive values, a reduced pair, or additional negative contrasts after separate validation. | `APPROVED`: paired fluid values 0.1, 0.4, 0.6, 0.6 and 0.8, plus rigid \(f_1=1\); no negative contrasts. | Already included in six material levels; each additional paired material would add 16 primary cases. |
| D04 | The implementation distinguishes finite fluid contrast from the exact sound-hard boundary; no measured-material table is in scope. | Canonical synthetic cases, named measured materials, or a hybrid; paired or crossed contrasts. | `APPROVED`: the five listed fluid pairs plus `material_model=rigid`, `f1=1`, `f0_applicable=false`. The rigid sentinel is not physical contrast. | Six materials produce 96 primary cases over the approved \(ka\)-distance grid; named materials would add unpriced property review. |
| D05 | Non-overlap requires \(d/a\ge2\); earlier campaigns span 2.1--10 and T11 anchors 2.5, 3 and 4. | Sparse anchors, linear/log grids, or a targeted near/far grid; exact contact was not selected. | `APPROVED`: \(d/a=\{2.1,2.25,2.5,3,4,6,8,10\}\). | Each distance contributes 12 primary cases; eight levels give 96. |
| D06 | Rotational covariance is expected, while the historical axis/diagonal cases change other parameters simultaneously. | One angle plus audits, two full angles, or three or more angular levels. | `APPROVED`: all primaries at \(\theta=0\), plus six \(\pi/4\) audits, one per material, balanced over both \(ka\) levels and the separation extremes, each linked to its exact zero-angle twin. | Six additional cases; a full second angle would add 96. |
| D07 | P0 requires response-blind stable IDs and order. | Full factorial, a balanced reduction, or a preregistered sequential design. | `APPROVED`: freeze exactly 96 primary plus six audit IDs in manifest order 1--102. | Between 408 and 2,040 case-order evaluations under D08. |
| D08 | Cap 9 left historical T11 channels unconfirmed; later work uses no stop below 5 and two-pass changes at \(10^{-5}\). | Caps 9, 13 or 21; any later extension would require a new pre-response version. | `APPROVED`: evaluate integer \(L_{\max}=2,\ldots,21\), no stop before 5, and require two consecutive passes in every applicable channel. | Four to 20 evaluations per case; maximum 2,040. Dense \(L=21\) work remains an unmeasured pilot risk. |
| D09 | P0 defines future \(B_E\) as a sum of independently converged isolated dimers; no aggregator exists. | Independent pair orders, recomputation at the largest order, or rejection without one common order. | `APPROVED`: retain each eligible pair's confirmed order; every pair must be eligible. Defer the common-order sensitivity audit to P1.3. | Minimizes production work; the deferred audit is not costed or executed in P1.1. |
| D10 | The established rule is two consecutive normalized RMS changes at \(10^{-5}\), channel by channel. | Interaction-only, all applicable channels, or a stricter sensitivity lane. | `APPROVED`: require total, interaction, external--scattered and scattered--scattered whenever applicable; null channels use the existing applicability rule. | No extra cases, but slow channels may consume the full 20 orders. |
| D11 | T13--T14.1 use finite values, balanced condition number below 10, four residual gates below \(10^{-12}\), and the planar-force tolerance. | Reuse, independently tighten, or weaken only with new numerical evidence. | `APPROVED`: reuse all established gates and later add dimension/mode and action--reaction/rotation identities; never gate on the legacy ill-scaled condition number. | Diagnostic work is negligible relative to solves; later identity audits remain in P1.3. |
| D12 | P0 forbids imputation and requires every attempt; an unconfirmed cap is not divergence. | Stop globally, continue with a complete ledger, or use a preregistered retry rule. | `APPROVED`: continue after case-local failures, retain stage/reason for every attempt, never change physics/tolerance/solver on retry, and report campaign resource exhaustion explicitly. | Ledger cost is negligible; continuing can spend the remaining fixed budget but prevents response selection. |
| D13 | There is no dimer timing; T14/T14.1 used one worker and one BLAS thread but are not dimer cost bounds. | One deterministic worker or bounded parallelism after audit; limits must remain provisional until pilot evidence. | `APPROVED_PROVISIONAL`: one worker, one BLAS thread, 4 GiB/case, 30 min/case and 12 h total. | At most 12 h campaign wall time; actual dimer time and RSS remain unknown. |
| D14 | P0 requires a separate timed pilot before campaign enablement. | Nominal, worst-cost, or multiple pilots, always outside confirmatory responses. | `APPROVED`: one separate disabled `development` rigid pilot at \(ka=0.1,d/a=2.1,\theta=0\); its force is excluded from P1.6 tables. | At most 20 pilot order evaluations, outside the 102 confirmatory cases. |
| D15 | P0 fixes Model-E interaction force, SI storage, optional normalization and evidence classes. | No alternative inside P1; another observable needs a new protocol. | `CONTRACT_FIXED`: preserve the P0 observable, normalization and classification contract. | Serialization only; no added solve. |
| D16 | Previous campaigns regenerate analysis without solves and audit a prespecified response subset. | No repeat, analysis-only repeat, selected audits, or a full duplicate. | `APPROVED`: two no-solver regenerations plus the six rotational-covariance audits; no full campaign duplicate. | Two analysis passes plus the six audits already counted in 102. |

Current gate: **HOLD_P1 — P1.2, P1.3 AND P1.4 REQUIRED**. No solver, pilot,
campaign, force or \(B_E\) operation is authorized by this record.
