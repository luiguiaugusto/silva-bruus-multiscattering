# P1.1 canonical-dimer decision record

Status: **DRAFT — SCIENTIFIC AUDIT REQUIRED**

Prepared: 2026-08-13

Decision owners: Luigui and ChatGPT Work

This record exposes the decisions required before the P1 canonical-dimer
campaign can be preregistered. It does not approve a scientific value, freeze
a case ID, implement \(B_E\), or authorize a solver call. Every row whose
state is `DECISION_REQUIRED` must be resolved explicitly in a later,
response-blind commit. The only P1 manifest in this change is planned,
disabled and non-executable.

## Fixed scope and existing evidence

The following items are constraints, not open choices:

- P1 covers two identical, non-overlapping, lossless spheres in a pressure
  nodal plane of an ideal unbounded fluid, with the
  \(e^{-i\omega t}\) convention and explicit SI \(E_0\)
  (`docs/CONVENTIONS.md`).
- The paper reference is `ModelENodalResult.interaction_forces_xyz`; total,
  external--scattered and scattered--scattered forces remain mandatory
  convergence and quality channels (`docs/PAPER_DATA_CONTRACT.md`).
- Historical B and \(B_L\) keep their meanings. \(B_E\) is reserved for a
  future sum of independently converged isolated Model-E dimers and is not
  implemented in P1.1 (`docs/PAPER_REDESIGN_INVENTORY.md`).
- Valid finite-fluid inputs satisfy \(f_0<1\) and \(-2<f_1<1\); exactly
  \(f_1=1\) selects the analytic rigid boundary. The conversion is
  \(\widetilde\kappa=1-f_0\),
  \(\widetilde\rho=(2+f_1)/[2(1-f_1)]\), and
  \(c_p/c_0=(\widetilde\rho\widetilde\kappa)^{-1/2}\)
  (`src/acoustic_ms/mie_scattering.py`).
- T10 validated exact isolated-sphere coefficients for
  \(10^{-3}\le ka\le0.1\), \(f_0=0\),
  \(f_1\in\{0.1,0.4,0.8,1\}\), and \(\ell\le5\). Its independent tests also
  include the finite-fluid point \((f_0,f_1)=(0.2,0.6)\).
- T11 contains only three historical dimer configurations: axis
  \((ka,f_1,d/a)=(0.1,0.8,2.5)\), diagonal \((0.05,0.4,4)\), and rigid
  \((0.1,1,3)\), all with \(f_0=0\) and orders 2--9. The diagonal case
  confirmed all four channels; the axis and rigid cases left at least one
  channel unconfirmed at the cap. These data are development evidence, not a
  new P1 sample.
- T13--T14.1 used a two-consecutive-change rule at \(10^{-5}\), no early stop
  below order 5, a standard cap of 13, balanced condition number below 10,
  four normalized numerical residuals below \(10^{-12}\), finite outputs and
  the planar-force tolerance \(128\epsilon_{mach}F_{max}\).
- P0 measured no dimer timing. T14 averaged 18.30 s per evaluated order and
  T14.1 averaged 145.06 s per evaluated order on much larger systems. These
  are arithmetic planning references only, not dimer estimates or bounds.

## Cost model used below

Let \(n_c\) be the number of dimer cases and \(n_L\) the number of evaluated
integer orders. A full Cartesian design has

\[
n_c=n_{ka}\,n_m\,n_d\,n_\theta+n_{audit}.
\]

For orders 2--13, \(n_L\le12\). The recommendation below would contain
\(2\times5\times8=80\) primary cases plus six rotation audits, or 86 cases,
and at most 1,032 case-order evaluations. Multiplying by the T14 arithmetic
reference gives 18,885.6 s (5.25 h); this number is deliberately reported
only as a large-cluster comparison. A timed P1 dimer pilot is required before
an actual wall-time or peak-memory budget can be frozen.

For a planar dimer, the balanced dense system dimension is
\(2L(L+1)/2=L(L+1)\). The complex128 storage of one dense matrix is about
0.13 MiB at \(L=9\), 0.51 MiB at \(L=13\), and 3.26 MiB at \(L=21\), before
translation tables, factorizations and Python/SciPy overhead. Cubic dense
solve work at \(L=13\) and \(L=21\) is approximately 8.3 and 135 times the
\(L=9\) work, respectively. These are scaling estimates, not measured peak
resources.

## Decision matrix

| ID | Decision and state | Existing evidence | Admissible options | Recommendation awaiting approval | Estimated cost impact |
|---|---|---|---|---|---|
| D01 | \(ka\) levels — `DECISION_REQUIRED` | T10 covers \(10^{-3}\) to 0.1; T11 dimers use 0.05 and 0.1; later campaigns fix 0.1. | One legacy level (0.1); two validated dimer levels (0.05, 0.1); a broader grid wholly inside \([10^{-3},0.1]\). Values above 0.1 require new domain validation. | Use \(\{0.05,0.1\}\): both have Model-E dimer precedent and expose frequency dependence without leaving the audited interval. | Two levels double cases relative to \(ka=0.1\); recommended share: 43 of 86 cases per level. |
| D02 | \(f_0\) levels — `DECISION_REQUIRED` | All historical Model-E campaigns use \(f_0=0\); T10 unit/oracle coverage includes \(f_0=0.2\); implementation accepts every finite \(f_0<1\). | Legacy-only \(\{0\}\); add one positive monopole point such as 0.2; add negative and positive points; reject values near 1 until separately stress-tested. | Include \(f_0=0\) in four legacy-linked material cases and one finite-fluid case at \(f_0=0.2\); do not make \(f_0\) and \(f_1\) a free Cartesian product. | One added paired material contributes 16 primary cases under the recommended \(2\,ka\times8\,d/a\) grid. |
| D03 | \(f_1\) levels — `DECISION_REQUIRED` | T10 uses 0.1, 0.4, 0.8 and exact rigid 1; Model-E campaigns otherwise emphasize positive 0.8. Domain is \(-2<f_1\le1\). | Preserve the four historical positive/rigid levels; use a reduced low/high pair; add negative contrast only with a separately justified claim domain. | Retain \(\{0.1,0.4,0.8,1\}\) and add 0.6 only as part of the nonzero-\(f_0\) finite-fluid pair. Defer negative contrast. | Five paired material cases yield 80 primary cases; a full 2-by-5 \(f_0,f_1\) factorial would double that before geometry and is not recommended. |
| D04 | Material cases and pairing — `DECISION_REQUIRED` | Model E represents identical lossless fluid spheres; \(f_1=1\) is rigid, not a finite-density fluid. Existing evidence is contrast-defined rather than a table of named experimental materials. | Synthetic contrast-defined cases; named measured material/host pairs with cited properties; hybrid. Material levels may be paired or fully crossed. | Use five paired canonical cases: `(0,0.1)`, `(0,0.4)`, `(0,0.8)`, `(0,1 rigid)`, and `(0.2,0.6 finite fluid)`. Label them mathematical material cases, not experimental substances. | Five cases per \((ka,d/a)\). Named materials would add literature/property uncertainty and review work not priced by solver timing. |
| D05 | Separations \(d/a\) — `DECISION_REQUIRED` | Non-overlap requires \(d/a\ge2\); T05/T06 used 2.1--10; historical Model-E dimers use 2.5, 3 and 4. | Sparse anchors; linear grid; log grid; targeted near/far grid. Exact contact 2 is admissible geometrically but lacks a clearance margin. | Use \(\{2.1,2.25,2.5,3,4,6,8,10\}\), preserving legacy anchors and resolving the strongest near-contact variation without claiming contact physics. | Eight levels contribute 10 primary cases each, or 80 total before rotation audits. Removing/adding one level changes the primary count by 10. |
| D06 | In-plane orientations — `DECISION_REQUIRED` | The nodal excitation and sphere physics imply rotational covariance; full/planar bases agree in an audited dimer, but T11 axis/diagonal cases confound orientation with other parameters. | One canonical axis plus explicit covariance tests; two angles \(0,\pi/4\) for every case; three or more angular levels. | Use \(\theta=0\) for the primary grid and six prespecified \(\pi/4\) audit cases spanning material, \(ka\), and separation extremes. Rotation must compare rotated vectors componentwise. | Six extra cases (7.5% over 80). A full second orientation adds 80 cases; three full angles add 160. |
| D07 | Case count, order and IDs — `DECISION_REQUIRED` | P0 requires response-blind stable IDs; no P1 ID is frozen. IDs may not depend on force, convergence or plot order. | Full factorial; balanced reduced design; sequential adaptive design only if frozen before responses. | Freeze exactly 86 cases after D01--D06 approval. Proposed pattern: `p1_dimer_<material>_ka<level>_d<level>_o<level>` with a separate deterministic `scale_order`; never reuse the disabled placeholder ID. | 86 cases, 344--1,032 order evaluations for early stops at 5 through cap 13. Any count change scales solve cost approximately linearly. |
| D08 | \(L_{max}\) range and stop — `DECISION_REQUIRED` | T11 evaluated 2--9 and left channels unconfirmed; T13--T14.1 use integer 2--13 and no early stop below 5. | Cap 9; standard cap 13; preregistered cap 21; adaptive cap changes only before response in a new version. | Evaluate every integer from 2, never stop below 5, stop only after all required channels pass twice, and cap at 13. A failed timed pilot may motivate a new pre-response version with cap 21; the campaign itself may not extend post hoc. | Maximum 12 orders/case. Relative dense-solve work is about 8.3 times \(L=9\) at \(L=13\); \(L=21\) is about 135 times \(L=9\). |
| D09 | Different confirmed orders in future \(B_E\) sums — `DECISION_REQUIRED` | P0 defines \(B_E\) as a sum of independently converged isolated dimers; no aggregator exists. | Each pair at its own confirmed order; recompute all pairs at the largest confirmed order; require one common order or reject the parent case. | Preserve independent per-pair confirmation and record every pair's final order; a future sum is eligible only if every pair is eligible. Also report a common-order sensitivity audit before P1 passes. | Independent orders minimize work. Recomputing \(P\) pairs at a common cap costs up to \(P\,n_L\); no P1.1 solve is allowed. |
| D10 | Convergence definition — `DECISION_REQUIRED` | Established rule is two consecutive applicable normalized RMS changes \(\le10^{-5}\), channel by channel; T11 shows the small scattered--scattered channel can lag. | Require interaction only; require all four channels; use \(10^{-5}\); add a stricter \(10^{-6}\) sensitivity lane. | Require total, interaction, external--scattered and scattered--scattered at \(10^{-5}\) for G1. Treat a numerically null channel as inapplicable under the existing scale-relative rule, not as zero change. | No extra cases; strict all-channel stopping can consume all 12 orders. A \(10^{-6}\) rule may require a higher cap and is not costed without pilot data. |
| D11 | Numerical quality gates — `DECISION_REQUIRED` | T13--T14.1 used finite outputs, balanced condition number <10, balanced backward/incident/scattering/decomposition residuals <\(10^{-12}\), and \(\max|F_z|\le128\epsilon F_{max}\). | Reuse all established limits; derive dimer-specific limits before response; weaken a gate only with independent numerical evidence. | Reuse the established gates unchanged and add dimension/mode-count and action--reaction/rotation identity checks in P1.3. Do not gate on the ill-scaled legacy condition number. | Diagnostics are negligible beside solves; independent audit solves add the prespecified audit-case cost. |
| D12 | Failures, missingness and exclusions — `DECISION_REQUIRED` | P0 requires every attempt retained; unconfirmed is not divergent; missing values are not zero or imputed. | Stop whole campaign at first numerical failure; continue and retain all; retry only under a preregistered deterministic rule. | Continue after case-local failure, retain one ledger row per attempt/order, set `eligible=false` with controlled stage/reason, and never retry with changed physics, tolerance or solver. Resource ceiling triggers campaign-level `INCONCLUSIVE_P1_RESOURCE_LIMIT`; scientific insufficiency gives `NO_GO_P2` or `INCONCLUSIVE_P1`. | Ledger cost is negligible. Continuing may spend the remaining frozen budget but prevents response-dependent selection. |
| D13 | Resource ceilings and parallelism — `DECISION_REQUIRED` | No dimer performance table exists. T14/T14.1 ran with one worker and one BLAS thread; their much larger peaks were about 1.0 and 7.4 GiB. | One deterministic worker; bounded parallel workers after a determinism audit. Candidate ceilings must cover pilot evidence with a declared margin. | Pilot with one worker/one BLAS thread. Provisional ceilings for approval: 4 GiB peak RSS, 30 min per case, and 12 h campaign wall time. Freeze final values from the P1.5 pilot before enabling the campaign. | One worker maximizes reproducibility but not throughput. Recommended grid's only available comparison is 5.25 h at the unrelated T14 per-order average; actual dimer time remains unknown. |
| D14 | Timed pilot and enable transition — `DECISION_REQUIRED` | P0 explicitly requires a preregistered timed dimer pilot before the full grid. P1.2--P1.4 must precede it. | Pilot a nominal case; pilot the expected worst case; use two pilots. Any scientific change after a response requires a new protocol version. | Use one disabled-now, later preregistered worst-cost case at largest \(ka\), smallest \(d/a\), and the hardest approved material. Record assembly, solve, postprocess and peak RSS. Enable only after P1.2/P1.3 pass, P1.4 freezes hashes, and all D01--D14 are approved. | At most 12 order evaluations. The T14 arithmetic comparison is 219.6 s, explicitly not a dimer forecast. |
| D15 | Observable, normalization and classification — `CONTRACT_FIXED` | P0 fixes Model-E interaction force, SI storage, optional \(F/(a^2E_0)\), and evidence classes. | No alternative inside P1. A different observable or physical domain requires a new task/protocol. | Store interaction as the primary response; retain all four channels; use `confirmatory_new` only after preregistration; preserve legacy evidence labels. | Metadata/serialization only; no additional solve count. |
| D16 | Reproducibility audit — `DECISION_REQUIRED` | Prior campaigns execute once, rerun analysis without solves, and audit a prespecified subset of final-order solves. P1 has no stochastic geometry. | No repeat; analysis-only repeat; selected final-order audit; full duplicate campaign. | Two byte-identical analysis-only runs plus six prespecified rotated/final-order audits. Do not duplicate the full 86-case campaign. | Six audit cases are already included in the recommended 86; a full duplicate would add 86 cases. |

## Approval and freeze gate

P1 remains disabled until one review records an explicit disposition for
D01--D14 and D16. Approval must state the chosen value/set, not merely
"approved as recommended" if any recommendation is changed. The freeze
commit must then:

1. replace the disabled placeholder with the exact ordered case list and
   immutable IDs;
2. record coordinates or a deterministic coordinate generator and hashes;
3. replace `TBD` provenance with the frozen manifest hash;
4. keep every case disabled until the separate P1.4 transition commit;
5. confirm P1.2 and P1.3 tests pass before the P1.5 pilot; and
6. preserve `results/` and `papers/` until an enabled campaign is explicitly
   authorized.

Current decision: **HOLD_P1 — DECISION_REQUIRED**. No \(B_E\), solver output,
force, campaign enablement or scientific freeze is part of P1.1.
