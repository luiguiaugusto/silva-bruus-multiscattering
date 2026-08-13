# Decisions

- Python 3.11+ is the official implementation language.
- Scientific routines live in the importable `src/acoustic_ms` package.
- Notebooks, when introduced, are demonstrations only and will never contain
  the sole version of a scientific routine.
- T01 depends only on NumPy at runtime and pytest for development.
- Matplotlib is an optional `plot` dependency, used only by the reproducible Figure 2 script; the scientific package does not import it.
- The T02 corrected formula is a published two-particle benchmark, not a general multiple-scattering solution.
- T03 uses dense NumPy linear algebra at `Lmax=1` only; SciPy and SymPy are runtime dependencies for special functions and Gaunt coefficients.
- `Lmax=1` truncates multipolar order but not the number of rescattering events; no radiation-force API is added in T03.
- T04 implements Model C at Rayleigh level using Eq. (22)/(27) cross terms only, with no scattered--scattered products.
- The T03 production solver remains at `Lmax=1`; T04 uses local evaluation through `ell=2` only and reports no three-particle force results.
- T05 is restricted to canonical N=3 trimers. Model B and C deliberately share the T04 solver and observable, so C-B isolates multibody rescattering at Lmax=1.
- The scalar nodal-plane oracle is test-only. No zero-total-force constraint is imposed on the global scalene interaction observable, and T05 does not measure multipolar correction or introduce Model D.

- T05.1 defines numerical nullity relative to the global configuration scale, \(128\,\epsilon_{\mathrm{mach}}F_{\mathrm{scale}}\), without an absolute floor. Correction amplitudes are RMS vector magnitudes per particle, not component RMS; only derived metrics and corresponding artifacts changed, while A, B, C and their equations remained unchanged.
- T05.1a, T05.1b, and T05.1c are exclusively documentary. Binary determinism is assessed in the same numerical environment.
- T06 reports only planar \(N=4\) results at Rayleigh \(L_{\max}=1\); Model D and higher scattered multipoles remain out of scope.
- For \(N=4\), \(C-B\) is not an irreducible four-body contribution: it equals the embedded three-body sum plus \(\boldsymbol{\Phi}^{(4)}\).
- Every connected term is built exclusively from subsets solved by the same Model C. Model A is a comparison baseline and does not define \(\boldsymbol{\Phi}^{(3)}\) or \(\boldsymbol{\Phi}^{(4)}\).
- The decomposition is vectorial. No zero-sum constraint is imposed on the approved irregular-quartet observable.

- T06.1 is post-processing of the already versioned T05/T06 CSVs. It performs no new trimer or quartet force sweep; the only additional Model C evaluation is one centered dimer for the \(N=2,3,4\) comparison.
- \(\Lambda_{\max}\) is an exploratory geometric diagnostic. Its grouped improvement is reported descriptively and no universal validity threshold is defined.
- Within each fixed-shape dilation family, \(\Lambda_{\max}=C_g\eta\), so exponent and log-space fit quality are necessarily unchanged from the \(\eta\) fit.
- T06.1 changes no force model, solver, connected-body definition, or protected T03--T06 artifact. T07, Model D, higher multipoles, and new cluster families remain outside its scope.

## T07 decisions

- Model D is a new API; all A--C APIs retain their fixed \(L_{\max}=1\) meaning.
- Every positive multipole uses its leading Rayleigh coefficient. This extends multipole order but is not an exact Mie T-matrix.
- The planar active basis uses \(\ell+m\) odd and is checked against the complete basis. It is not reduced indiscriminately to odd \(\ell\).
- The balanced matrix is used for solution and conditioning; convergence is judged by forces and connected terms, not by the ill-scaled raw condition number.
- The strict odd-order dimer branch used to derive Eq. (30) is retained only as an independent benchmark. The general planar Model D also preserves symmetry-allowed even-\(\ell\), odd-\(m\) channels.
- Connected terms are recomputed by subset inclusion--exclusion at a common \(L\). Model A does not enter that definition.
- No 1,920-case T05 or T06 sweep was rerun. No universal claim that \(L=5\) is sufficient is made.
- T08 leaves every established force solver unchanged. Its \(B_L\) baseline is a diagnostic sum of isolated Model-D dimers at the same \(L\), and does not rename or redefine historical Model B.
- T08 calibrates exclusively on \(N\leq4\). The \(N=6,10\) cases are a locked external holdout and are never used for predictor selection, fits, or empirical thresholds.
- Predictor selection uses leave-\((N,\mathrm{family})\)-out log-RMSE for \(\varepsilon_A\). The selected predictor is \(\rho_1\), the spectral radius of the balanced dipolar rescattering operator.
- The 1%, 5%, and 10% cutoffs are conservative empirical nodal-plane thresholds within the sampled domain. Their successful holdout diagnosis is not a universal criterion or an analytic bound.
- Nonconverged-at-limit cases are labeled unconfirmed, not divergent, and are excluded from fits and thresholds. Raw conditioning is retained only as a scaling diagnostic; convergence decisions use successive force changes.
- T08 closes and freezes the computational datasets intended for the article. No T05/T06 sweep was rerun, and later work must treat the committed T08 raw tables as immutable provenance.

## T09 decisions

- T09 is the newly approved analytical foundation of \(\rho_1\), superseding the unrelated preliminary T09 label in the first project roadmap.
- Python/SymPy remains the official and sufficient implementation; Wolfram Mathematica is not a project dependency.
- The closed \(N\times N\) operator is an independent analytical reconstruction of the existing balanced \(L=1\) matrix. It does not change the T03, T07, or T08 solvers.
- The exact retarded entry is used for \(\rho_1\); the inverse-cube operator is only its near-field limit.
- \(\rho_1<1\) is used as the necessary and sufficient convergence condition for the finite-dimensional matrix Neumann series, not as a universal force-error threshold.
- The approximately linear force-error law is an asymptotic order statement. Its prefactor and the T08 fitted exponent remain empirical.
- Non-normality is reported explicitly. The spectral radius controls asymptotic convergence, while induced norms delimit possible finite-order amplification.
- T09 performs no new force sweep and does not rewrite any T01--T08 artifact.

## T10 decisions

- Exact fluid-sphere Mie coefficients live in a separate single-particle API;
  Models A--D and the multipolar solver retain their approved meanings.
- The exact point \(f_1=1\) uses the analytic rigid boundary condition. No
  finite surrogate density, clipping, or near-one tolerance is introduced.
- Validation uses an independent boundary-condition system, Rayleigh
  asymptotics, the rigid limit, and lossless unitarity.
- The production audit is restricted to \(10^{-3}\le ka\le0.1\) and
  \(\ell\le5\). Acceptance of larger \(ka\) by the isolated API is not a
  convergence guarantee for a caller-selected \(L_{\max}\).
- Exact isolated-sphere coefficients do not imply a complete collective-force
  theory. Integration with global multiple scattering and
  `scattered--scattered` force terms is deferred to T11.

## T11 decisions

- Model E is a new, isolated API combining the exact T10 T-matrix, global
  multiple scattering, and the complete multipolar force. Models A--D and all
  their established artifacts remain unchanged.
- T11.1 changes only the numerical solution path. Production solves the
  principal-square-root-balanced system
  \((I-D^{1/2}UD^{1/2})q=D^{1/2}a\), reconstructs \(d=D^{1/2}q\), and then
  reconstructs \(b=a+Ud\) without division by \(D^{1/2}\).
- The effective-incident and scattered systems remain explicit diagnostics.
  Their legacy public attributes retain their original meanings; neither raw
  condition number is interpreted as physical multipole divergence.
- The complete force includes both the quadratic recoil in \(\Gamma_n\) and
  the distinct scattered--scattered incident-field channel \(\mathcal F[c]\).
- An independent stress-tensor surface integral is the normalization and sign
  oracle. No empirical prefactor or sign adjustment is permitted.
- Model-E convergence is assessed separately for total, interaction,
  external--scattered, and scattered--scattered forces. Cases without two
  successive changes below \(10^{-5}\) by \(L_{\max}=9\) are retained and
  explicitly marked unconfirmed.
- Model E is not a recalibration of \(\rho_1\). T12, its sentinel campaign,
  and the T13--T14 holdout remain unopened.

## T12 decisions

- T12 is a preregistered audit of exactly 28 T08 calibration cases with
  \(N\leq4\). It neither recalibrates the frozen \(\rho_1\) law nor changes its
  1%, 5%, and 10% thresholds.
- The reference is the complete Model-E interaction force, not total force.
  All comparisons are three-dimensional and retain any computed \(F_z\).
- Frozen A and D vectors are independently reproduced through their public
  APIs before any E solve. No T08 sweep or earlier artifact is regenerated.
- E convergence is channel-specific, uses two applicable consecutive changes
  below \(10^{-5}\), reaches at least \(L_{\max}=5\), and is capped at 13.
  Cases reaching the cap without confirmation remain `unconfirmed`.
- The A–D–E mechanism decomposition is vectorial. RMS mechanism amplitudes
  cannot be treated as additive force fractions.
- The T12 gate is an internal decision about whether a separately authorized
  T13 may start. It is not a universal validation of \(\rho_1\), and T13–T14
  remain outside the scope of this implementation.

## T12.1 decisions

- T12.1 extends exactly ten preregistered T12 calibration sentinels. It copies
  the immutable \(L=2,\ldots,13\) records and computes only \(L=14,\ldots,21\)
  as necessary.
- A channel is confirmed only after two successive applicable relative changes
  no larger than \(10^{-5}\). Reaching \(L=21\) without this evidence is
  `unconfirmed_at_21`, never “divergent”.
- Models A--E, the definition of \(\rho_1\), the T08 fit, its thresholds, and
  all T01--T12 artifacts remain unchanged.
- Mechanism analysis is performed on signed vector fields. Cosines,
  projections, and RMS-amplitude ratios are not additive scalar force
  fractions.
- Candidate comparison uses deterministic leave-\((N,\mathrm{family})\)-out
  validation. No points are discarded and no predictor is augmented after
  inspecting residuals.
- The resulting recommendation is
  `READY_T12_2_RHO1_RECALIBRATION_STUDY`: it authorizes only a separately
  specified recalibration study on the 28 sentinels. It is not `GO_T13`, does
  not open the \(N=6,10\) holdout, and is not a universal validity claim.

## T12.2 decisions

- The only confirmatory candidate is the same unweighted power law in
  \(\rho_1\), fitted in log space. No alternative feature or model is tested.
- Generalization is assessed exclusively by seven deterministic
  leave-\((N,\mathrm{family})\)-out folds. The fit on all 28 cases is
  descriptive and cannot support the gate.
- Fold-specific thresholds are inverted only from training coefficients. No
  residual margin, clipping, outlier removal, or post-hoc safety factor is
  introduced.
- The 10% threshold retains one false-safe dimer. Consequently the exact gate
  is `NO_GO_T13_RHO1_NOT_QUANTITATIVE`, even though all other criteria pass.
- The final fit is retained as a documented candidate, not as an approved
  autonomous criterion. T13, T14, and the \(N=6,10\) holdout remain unopened.

## T12.3 decisions

- T12.3 is statistical post-processing of the 28 frozen, confirmed
  \(N\leq4\) sentinels. It performs no new acoustic solve and does not read the
  external \(N=6,10\) holdout.
- Exactly two candidates are admitted. M1 is an unweighted log-space power law
  in the geometric sum \(\Lambda_{\max}\). M2 adds \(\rho_1\) and can enter
  the gate only if M1 fails. P0 and P3 remain frozen baselines.
- Every decision metric is outer LOGO. A nested LOGO inside the six training
  groups calibrates one conservative residual factor; the external group
  enters neither the fit nor this factor. Safety uses the strict comparison
  \(\widehat\varepsilon_{\mathrm{safe}}<\tau\).
- The operational M2 identifiability rule is fixed before its gate: singular
  design, standardized condition above \(10^3\), or simultaneous positive
  mechanistic exponents in fewer than 80% of outer folds yields
  `UNSTABLE_COLLINEARITY`. No regularization is introduced.
- M1 satisfies all preregistered criteria, including zero conservative false
  safe and the 3/8/12 antivacuity counts. The exact result is
  `GO_T13_VALIDATE_LAMBDA_MAX`.
- This result does not validate \(\Lambda_{\max}\) externally. It authorizes
  only a separately specified T13 using the still unopened \(N=6,10\)
  holdout. T13 and T14 were not started.

## T13 decisions

- The exact 24-case selection, M1/P3 coefficients, margins, thresholds, gate,
  and three blind CSVs were committed and pushed before any new \(N=6,10\)
  Model-E solve.
- Selection uses only T08 holdout metadata, \(\Lambda_{\max}\), \(\rho_1\),
  and confirmation of the preexisting D reference. No A/D/E force or error
  participates.
- The complete Model-E interaction force is the sole external reference. A
  and D are independently audited diagnostics.
- Sufficiency precedes the scientific gate. An insufficient campaign is
  `INCONCLUSIVE`, not a scientific failure. P3 cannot rescue M1.
- The measured result is `PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX`; this
  validates the frozen criterion only on the sampled planar \(N=6,10\)
  domain. No coefficient was refitted and no universal claim is made.
- The campaign is executed once. Determinism is tested by two analysis-only
  passes, while eight stratified final-order solves provide an independent
  audit without repeating all 24 cases.
- T14 is only recommended by the frozen gate. It is not implemented or
  started by T13.

## T14 decisions

- T14 uses two chronological commits. The response-blind geometry,
  predictions, protocol, code, and tests were pushed before any of the 24 new
  Model-E solves. The four blind CSVs remained byte-identical afterward.
- The M1 coefficients, P3 comparator, safety factors, tolerances, target
  levels, case IDs, convergence rule, and gate were frozen. No fitting,
  outlier removal, or post-response selection was performed.
- The sole reference is the complete three-dimensional Model-E interaction
  force. P3 is transparent diagnostic evidence and cannot rescue or reject M1.
- Sufficiency is evaluated before science. All 24 cases were eligible and
  every literal scientific criterion passed, yielding
  `PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX` and
  `GO_T15_SYNTHESIS_AND_MANUSCRIPT`.
- This result supports the frozen rule only for the sampled planar families,
  positive fixed contrast, \(ka=0.1\), and \(N\leq28\). It does not establish
  a universal error bound or authorize recalibration within T14.

## T14.1 decisions

- T14.1 used two chronological commits. The complete response-blind code,
  24-case manifest, local-coupling vectors, M1/P3 predictions, protocol, and
  95 prior-artifact hashes were pushed in commit
  `538142b638dd59768d26bd16809b1def83bfdf8c` before the first T14.1 Model-E
  solve.
- The sample contains exactly 12 cases at \(N=45\) and 12 at \(N=105\), with
  eight per family and six per target level. No case, coefficient, margin,
  tolerance, convergence rule, or gate was changed after response revelation.
- The Model-E campaign used one worker and one BLAS thread. All 24 cases were
  eligible, all numerical gates passed, and the independent six-case
  post-revelation audit reproduced the official final-order results.
- M1 passed every literal criterion, yielding
  `PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX` and
  `GO_T15_SYNTHESIS_AND_MANUSCRIPT`. P3 remains reported as a comparator only.
- The matched-size diagnostic is `NO_SYSTEMATIC_DETERIORATION`; it is
  descriptive and does not add an \(N\)-dependent term to M1.
- The conclusion is limited to the prescribed deterministic planar families,
  \(N=45,105\), \(ka=0.1\), \(f_0=0\), \(f_1=0.8\), identical spheres, and
  the approved complete Model-E interaction force. It is not a universal
  theorem, does not authorize extrapolation or recalibration, and does not
  start T15.

## Paper redesign P0 decisions

- P0 inventories the actual checkout at post-T14.1 even though its original
  specification says post-T13. No later evidence is erased or relabeled.
- T01–T14.1 artifacts remain immutable. Their paper-protocol labels are
  `exploratory`, `development`, or `legacy_validation`; only a future
  response-blind P1–P6 campaign may be `confirmatory_new`.
- `docs/PAPER_CONFIRMATORY_PROTOCOL.md` is deliberately
  `DRAFT — NOT FROZEN`. P0 introduces no scientific parameter, fit,
  tolerance, force result or campaign.
- The editorial name \(B_E\) is reserved for a future sum of independently
  converged isolated Model-E dimer interaction forces. Historical B and
  \(B_L\) retain their approved meanings.
- P1 must be the canonical dimer benchmark. Model-E connected terms and any
  \(\Phi_E^{(5)}\) implementation remain deferred until P1 passes and P2 is
  separately frozen.
- New data have three immutable levels: raw solver output, derived analysis
  data and plot-ready data. Legacy CSVs are adapted into new files and are
  never rewritten.
- Manifest examples use JSON syntax with a `.yaml` suffix because JSON is a
  YAML subset. The P0 validator uses only the standard library and the
  versioned lightweight schema subset.
- `diagnostic-v1` is a reversible Matplotlib `rc_context`, not a final journal
  style. It supports one-/two-column physical sizes and PDF/SVG/PNG output.
- The P0 gate is `GO_P1_WITH_CONDITIONS`: scientific TBDs, B_E validation and
  a timed dimer pilot must precede any P1 campaign.

## P1.1 frozen decisions

- P0 was finalized by GitHub PR #1 at merge commit
  `926e639fe2d327eacd09a2542208500891399687`; its 17-file scope changed no
  file under `results/` or `papers/`.
- P1.1 freezes \(ka=\{0.05,0.1\}\), six ordered material cases, eight
  separations, 96 primary zero-angle dimers and six \(\pi/4\) covariance
  audits balanced across both \(ka\) values and the separation extremes.
  `campaigns/p1/campaign_manifest.yaml` contains the exact 102 IDs and
  immutable order.
- The rigid boundary uses `material_model=rigid`, `f1=1` and
  `f0_applicable=false`. Its stored finite `f0=0` is only an API sentinel, not
  a physical contrast.
- Campaign schema `1.1.0` moves `ka`, `k_rad_m`, material, contrasts,
  separation and angle to each case while retaining campaign constants and
  policies globally. Schema `1.0.0` remains accepted unchanged.
- Orders 2--21, minimum stop 5, all-applicable-channel two-step convergence,
  established numerical gates and the no-imputation failure policy are
  approved. The common-order \(B_E\) audit is assigned to P1.3.
- `campaigns/p1/pilot_manifest.yaml` is a separate disabled `development`
  manifest for the rigid \(ka=0.1,d/a=2.1,\theta=0\) resource pilot. Its force
  is prohibited from P1.6 scientific tables.
- Resource ceilings remain explicitly provisional: one worker/thread, 4 GiB
  per case, 30 minutes per case and 12 hours total.
- P1.1 adds no `B_E` code, solver runner, force output or campaign response.
  All 103 cases across the two manifests remain disabled, hashes remain `TBD`
  until P1.4, and the gate is `HOLD_P1 — P1.2, P1.3 AND P1.4 REQUIRED`.


## P1.2 decisions

- P1.1 was finalized from PR #2 only after its head was confirmed at
  `24ec933f366cb4950ad4050d83ce804d89d4eb43`. It was marked ready and merged
  with merge commit `4a5b58408dc40302568758b2bdea54701beb4747`.
- \(B_E\) is implemented by `solve_model_be_nodal` as the sum of the final
  Model-E `interaction_forces_xyz` from every isolated unordered pair
  \(i<j\), preserving original particle orientation and index association.
- Pair order is lexicographic and deterministic. Each pair converges
  independently over orders 2--21, cannot stop before 5, and requires two
  consecutive passes for every applicable force channel. Numerically null
  channels retain explicit non-applicability.
- The established Model-E finite, conditioning, residual, mode-dimension and
  planar gates are centralized in
  `evaluate_model_e_numerical_diagnostics` and reused; their scientific
  thresholds are unchanged.
- The ledger records individual dimer forces, attempted/evaluated/final/failed
  orders, full channel histories, diagnostics and explicit failure reasons.
  Later pairs are audited after a local failure, but the global force is
  unavailable unless every pair is eligible. No partial force is imputed.
- Historical Model B and \(B_L\) are unchanged. Manifests remain disabled and
  no pilot, campaign, production calculation, `results/` or `papers/` output
  is part of P1.2.
- Unit tests use injected fake solvers. Full physical identities,
  common-order sensitivity, rotation/reflection, action--reaction and
  asymptotic limits remain assigned to P1.3. The implementation decision is
  `GO_P1.3`; this does not authorize P1.4 or any execution.
