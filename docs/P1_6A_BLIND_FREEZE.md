# P1.6A — response-blind pre-campaign freeze

Date: 2026-08-20
Status: **frozen; `GO_P1.6B_EXECUTE`; no confirmatory case executed**

## Provenance boundary

PR #6 was audited at exactly
`e5f4efebbdd4e9983ef075f3fa038b875dd28bda`: two commits, 17 files and the
five P1.5 hashes recorded in `docs/P1_5_TIMED_PILOT.md`. It was merged into
`main` as `0e4f643ef8161af57af41c6600944eaaf6f8719a` before this branch was
created. The P1.5 manifest and all five pilot artifacts are byte-unchanged.

The P1.4 confirmatory lock remains a historical, non-executable checkpoint:

```text
9d360de6e61d901cff3f84c477f367773251103db12386dbb8156bd1ec2addca
```

P1.6A changes only confirmatory provenance, case `enabled` flags, resource
status/total wall limit and the self-hash. After normalizing only the
manifest's own hash field to 64 ASCII zeroes, the new exact-byte lock is:

```text
3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe
```

The title still contains the historical word `blocked` because title bytes
were explicitly outside the authorized P1.6A diff. It is not an execution
flag. The normative machine checks are the public lock, `status`, 102
`enabled=true` values, frozen resources and the explicit P1.6B branch guard.
P1.6A itself performs no execution.

## Frozen campaign and resource budget

- cases/order: exactly 102 IDs at `case_order=1..102`;
- evidence: 96 primaries included in scientific tables and six rotational
  audits excluded from them;
- design: two `ka` values by six materials by eight primary separations, plus
  six frozen audit–twin links;
- numerical policy: `L=2..21`, no stop before 5, tolerance `1e-5`, and the two
  latest applicable changes must pass simultaneously in all applicable
  `total`, `interaction`, `external_scattered` and `scattered_scattered`
  channels;
- concurrency: one worker and one BLAS thread;
- local limits: 1800 s and 4 GiB per case;
- global limit: 64800 s (18 h), with `limits_status=frozen`.

The maximum design contains `102 × 20 = 2040` case-order evaluations. The sum
of all local wall ceilings is `102 × 1800 = 183600 s = 51 h`; the independent
global ceiling is 64800 s, or an average of `635.2941176470588 s/case` if all
102 are attempted. Because execution is serial, the instantaneous memory
ceiling is 4 GiB, not 102 times that value. Reaching either resource ceiling
is recorded and yields `INCONCLUSIVE_P1`; it does not authorize retry or
retuning.

## Single-attempt runner contract

`run_p1_6_campaign` consumes the locked manifest in order 1..102 and accepts
an injected case executor. Before a call it atomically records `started`; after
the call it atomically records `completed` or `interrupted`. Each case has
`attempt_count` in `{0,1}`. On resume, a stale `started` case is permanently
converted to `interrupted`; only `never_started` cases may run. Local failure
does not stop later cases. A closed campaign, existing outputs, changed lock,
duplicate attempt or output overwrite is rejected.

`execute_model_e_case` constructs the frozen centered dimer, invokes
`solve_model_be_nodal` once and collects each Model-E order once. For `N=2`,
the final E interaction force is taken from that same pair ledger and reused
as `B_E`; no second dimer solve is made. Model A is evaluated analytically
afterward. `execute_model_e_case_with_limits` reuses the P1.5-tested wall and
address-space guard. P1.6A tests these paths only with injected fake solvers.

Every checkpoint retains attempted/evaluated orders, the four channel vectors
and convergence fields, per-order numerical diagnostics, order/case timing,
peak RSS, Model A, B_E and E vectors, convergence, eligibility, stage and
reason. The pure module `p1_campaign_artifacts.py` imports neither NumPy,
SciPy nor any `acoustic_ms` solver module. It deterministically produces raw,
derived, plot, failure and performance CSV bytes. Two regenerations must be
identical before publication, and publication refuses any existing target.

## Gate G1 frozen before response

For normalized vectors `F/(a^2 E0)`, the identity and rotation residual use

```text
max_i ||F_observed,i - F_reference,i||_2
-------------------------------------------------
max(1, max_i ||F_observed,i||_2, max_i ||F_reference,i||_2)
```

with budget `1e-12`. G1 is evaluated with these immutable rules:

1. all 102 cases have exactly one attempted case-level call;
2. every one of the 12 `(ka, material_id)` strata has at least one eligible
   primary;
3. all six rotational audits and their frozen twins are eligible;
4. every eligible dimer satisfies `B_E=E` within `1e-12`;
5. all six audit–twin rotations satisfy covariance within `1e-12`;
6. only eligible primaries enter scientific/plot metrics; every exclusion is
   retained in the failure and performance tables;
7. `|B_E-A|` is a scientific response, never a G1 acceptance threshold.

Contract, identity or covariance failure has precedence and yields
`FAIL_G1 / NO_GO_P2`. Missing attempts, stratum/audit coverage or resource
exhaustion yields `INCONCLUSIVE_G1 / INCONCLUSIVE_P1`. Only all valid gates
yield `PASS_G1 / GO_P2`.

## P1.6A decision

Manifest recomputation, mutation guards, checkpoint/resume, one-attempt
semantics, local/global limits, schemas, deterministic hashes and two
no-solver regenerations pass with fake responses. The focused suite reports
**51 passed in 2.02 s**; the complete suite with warnings as errors reports
**591 passed, 1 xfailed in 667.88 s**. No confirmatory checkpoint
or response file exists in this commit, and `results/` and `papers/` are
unchanged. Decision: **`GO_P1.6B_EXECUTE`**, pending audit of this draft PR.
