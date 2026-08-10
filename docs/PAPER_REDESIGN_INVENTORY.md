# Paper redesign inventory

Status: **P0 snapshot — no scientific result generated**

## Repository snapshot

- Initial commit: `e98080da520ec4f1f36b41ece2b76bc281df3a92`
  (`feat: validate frozen lambda max at large N`, 2026-08-01).
- The checkout is not merely post-T13: it includes completed T14 and T14.1.
  P0 therefore inventories the actual post-T14.1 state.
- Initial branch: `main`, tracking `origin/main`.
- Initial baseline: 464 passed in 87.93 s with
  `.venv/bin/python -m pytest`. The bare `.venv/bin/pytest` launcher failed
  collection because it omitted the repository root from `sys.path`.
- The initial worktree was mixed. Twenty-seven tracked prompt/task records were
  already deleted from the root; equivalent/additional records were untracked
  in `tarefas/`. `APOSTILA_SILVA_BRUUS_VOLUME_II_T09_A_T14_1.pdf` and
  `base_completa_silva_bruus_multiscattering.pdf` were also untracked.
  These user-owned changes are outside P0.
- The historical-result baseline is the 95-path path/size/SHA-256 manifest
  `results/data/t14_1_prior_artifact_hashes.csv`, enforced by
  `tests/test_t14_1_preregistration.py::test_preregistration_is_blind_idempotent_and_preserves_prior_hashes`.

## Scientific modules and responsibilities

Every code statement names its defining path and symbol.

| Module | Implemented responsibility |
|---|---|
| `contrasts.py` | Material contrast conversions \(f_0,f_1\) (`src/acoustic_ms/contrasts.py::monopole_contrast`, `dipole_contrast`). |
| `silva_bruus.py` | Model A: nodal Silva–Bruus signed pair force and cluster sum (`src/acoustic_ms/silva_bruus.py::nodal_pair_force_magnitude`, `nodal_pair_force_on_probe`, `nodal_pair_forces`). |
| `corrected_pair.py` | Fifth-order isolated-pair correction, not a global solver (`src/acoustic_ms/corrected_pair.py::corrected_pair_coefficients`, `corrected_nodal_pair_forces`). |
| `multipoles.py` | Complete multipole indexing (`src/acoustic_ms/multipoles.py::mode_index`, `mode_from_index`, `modes`). |
| `special.py` | Spherical functions and coordinates (`src/acoustic_ms/special.py::spherical_hankel1`, `spherical_harmonic`, `cartesian_to_spherical`). |
| `gaunt.py` | Cached 3j-based Gaunt coefficient (`src/acoustic_ms/gaunt.py::gaunt_coefficient`). |
| `translation.py` | Target-from-source re-expansion (`src/acoustic_ms/translation.py::separation_coefficient`, `translation_matrix`). |
| `incident.py` | Nodal standing-wave BSCs (`src/acoustic_ms/incident.py::nodal_standing_wave_coefficients`). |
| `scattering.py` | Monopole/dipole leading-Rayleigh coefficients (`src/acoustic_ms/scattering.py::rayleigh_scattering_coefficients`). |
| `solver.py` | Coupled Rayleigh \(L=1\) coefficient solve (`src/acoustic_ms/solver.py::solve_rayleigh_nodal`, `RayleighNodalSolution`). |
| `force.py` | Model C \(L=1\) nodal external–scattered force (`src/acoustic_ms/force.py::solve_rayleigh_nodal_interaction_forces`, `RayleighNodalInteractionResult`). |
| `geometries.py` | Canonical trimer/quartet coordinates (`src/acoustic_ms/geometries.py::linear_trimer`, `equilateral_trimer`, `scalene_trimer`, `linear_quartet`, `square_quartet`, `irregular_quartet`). |
| `comparison.py` | Historical A/B/C; B sums isolated Model-C dimers (`src/acoustic_ms/comparison.py::compare_nodal_force_models`, `NodalForceModelComparison`). |
| `cluster_expansion.py` | Model-C connected expansion at \(N=4,L=1\) (`src/acoustic_ms/cluster_expansion.py::decompose_nodal_quartet`, `NodalQuartetBodyExpansion`). |
| `metrics.py` | 2-D symmetric, RMS and angular metrics (`src/acoustic_ms/metrics.py::symmetric_particle_errors`, `rms_relative_error`, `rms_vector_magnitude`, `angular_errors_degrees`). |
| `scaling.py` | \(\eta\), \(\Lambda_{\max}\), log-space fits (`src/acoustic_ms/scaling.py::coupling_eta`, `maximum_geometric_coupling`, `fit_power_law`). |
| `multipolar_scattering.py` | Leading-Rayleigh coefficients for D (`src/acoustic_ms/multipolar_scattering.py::rayleigh_multipolar_scattering_coefficients`). |
| `multipolar_solver.py` | Balanced global Model-D solve (`src/acoustic_ms/multipolar_solver.py::solve_multipolar_nodal`, `MultipolarNodalSolution`). |
| `model_d.py` | Model-D interaction force and comparison (`src/acoustic_ms/model_d.py::solve_multipolar_nodal_interaction_forces`, `compare_nodal_model_d`). |
| `multipolar_expansion.py` | Common-\(L\) D inclusion–exclusion for \(2\le N\le4\) (`src/acoustic_ms/multipolar_expansion.py::decompose_multipolar_cluster`). |
| `cluster_families.py` | Transferability families through \(N=10\) (`src/acoustic_ms/cluster_families.py::cluster_family`, `enumerate_transferability_configurations`). |
| `transferability.py` | \(B_L,\rho_1\), convergence, fits and thresholds (`src/acoustic_ms/transferability.py::matched_multipolar_pairwise_baseline`, `spectral_radius_l1`, `two_step_converged`, `select_predictor_by_group_cv`, `conservative_threshold`). |
| `rho_foundation.py` | Closed \(L=1\) operators and Neumann audit (`src/acoustic_ms/rho_foundation.py::dipolar_coupling_entry`, `dipolar_balanced_coupling_matrix`, `dipolar_coupling_diagnostics`, `neumann_partial_solutions`). |
| `mie_scattering.py` | Exact isolated-sphere Mie/rigid coefficients (`src/acoustic_ms/mie_scattering.py::fluid_sphere_mie_scattering_coefficients`, `rigid_sphere_scattering_coefficients`, `mie_scattering_coefficients_from_contrasts`). |
| `mie_multiparticle.py` | Balanced exact-Mie global solve (`src/acoustic_ms/mie_multiparticle.py::solve_mie_multiparticle_nodal`, `MieMultiparticleSolution`). |
| `complete_force.py` | Complete 3-D radiation-force functional (`src/acoustic_ms/complete_force.py::complete_radiation_force_from_bsc`). |
| `model_e.py` | Model E and five force channels (`src/acoustic_ms/model_e.py::solve_model_e_nodal`, `ModelENodalResult`). |
| `model_e_comparison.py` | 3-D E errors and mechanism identity (`src/acoustic_ms/model_e_comparison.py::normalized_rms_error_xyz`, `compare_model_e_forces`). |
| `rho1_model_e_diagnostics.py` | T12.1 mechanisms, LOGO and tails (`src/acoustic_ms/rho1_model_e_diagnostics.py::mechanism_diagnostics`, `out_of_fold_metrics`, `convergence_tail_diagnostics`). |
| `rho1_model_e_recalibration.py` | T12.2 LOGO recalibration/gate (`src/acoustic_ms/rho1_model_e_recalibration.py::logo_power_law_predictions`, `classify_logo_safety`, `evaluate_recalibration_gate`). |
| `mechanistic_validity.py` | T12.3 M1/M2 nested-LOGO gate (`src/acoustic_ms/mechanistic_validity.py::nested_logo_predictions`, `audit_safety_thresholds`, `evaluate_mechanistic_gate`). |
| `external_validation.py` | T13 blind selection/prediction/gate (`src/acoustic_ms/external_validation.py::select_external_validation_cases`, `frozen_external_predictions`, `external_eligibility_mask`, `evaluate_external_validation_gate`). |
| `scale_out_validation.py` | T14 \(N=15,28\) scaling/gate (`src/acoustic_ms/scale_out_validation.py::build_scale_out_cases`, `geometric_coupling_sum`, `evaluate_scale_out_gate`). |
| `large_n_validation.py` | T14.1 \(N=45,105\) coupling/gate (`src/acoustic_ms/large_n_validation.py::build_large_n_cases`, `local_geometric_coupling`, `evaluate_large_n_gate`). |

P0 adds non-scientific infrastructure only:
`src/acoustic_ms/paper_pipeline.py::validate_manifest_file` and
`src/acoustic_ms/plot_style.py::diagnostic_rc_context`.

## Scripts by task

| Task | Scripts present |
|---|---|
| T01 | No script; importable implementation and tests only. |
| T02–T05 | `scripts/reproduce_figure_2.py`, `validate_t03_solver.py`, `validate_t04_force.py`, `validate_t05_trimers.py`. |
| T06/T06.1 | `scripts/validate_t06_quartets.py`, `analyze_t06_scaling.py`. |
| T07 | `scripts/validate_t07_multipolar.py`. |
| T08 | `scripts/run_t08_transferability.py`, `analyze_t08_transferability.py`. |
| T09 | `scripts/analyze_t09_rho_foundation.py`. |
| T10 | `scripts/analyze_t10_mie_rayleigh.py`. |
| T11/T11.1 | `scripts/analyze_t11_model_e.py`, `t11_stress_oracle.py`, `analyze_t11_1_model_e_stability.py`. |
| T12–T12.3 | `scripts/analyze_t12_model_e_sentinels.py`, `analyze_t12_1_rho1_failure.py`, `analyze_t12_2_rho1_recalibration.py`, `analyze_t12_3_mechanistic_validity.py`. |
| T13 | `scripts/preregister_t13_external_validation.py`, `run_t13_external_validation.py`, `analyze_t13_external_validation.py`. |
| T14/T14.1 | `scripts/preregister_t14_scale_out.py`, `run_t14_scale_out.py`, `analyze_t14_scale_out.py`, `preregister_t14_1_large_n.py`, `run_t14_1_large_n.py`, `analyze_t14_1_large_n.py`. |

## Existing datasets, tables, manifests, and figures

All paths below existed initially. CSVs are immutable legacy tables.

- T02–T04 data: `figure_2_relative_error.csv`,
  `t03_solver_validation.csv`, `t04_pair_force_validation.csv`;
  figure: `figure_2_relative_error.png`.
- T05 data: `t05_trimer_regression.csv`, `t05_trimer_sweep.csv`;
  figure: `t05_trimer_model_errors.png`.
- T06/T06.1 data: `t06_quartet_regression.csv`,
  `t06_quartet_sweep.csv`, `t06_1_body_order_summary.csv`,
  `t06_1_collapse_summary.csv`, `t06_1_scaling_fits.csv`;
  figures: `t06_quartet_body_decomposition.png`,
  `t06_quartet_model_errors.png`, `t06_1_eta_scaling.png`,
  `t06_1_lambda_scaling.png`.
- T07 data: `t07_pair_analytic_validation.csv`,
  `t07_dimer_convergence.csv`, `t07_cluster_convergence.csv`;
  figures: `t07_dimer_convergence.png`, `t07_cluster_convergence.png`.
- T08 data: `t08_cases.csv`, `t08_forces.csv`,
  `t08_convergence.csv`, `t08_predictor_fits.csv`,
  `t08_validity_thresholds.csv`;
  figures: `t08_predictor_comparison.png`, `t08_transferability.png`.
- T09 data: `t09_analytic_summary.csv`, `t09_neumann_convergence.csv`,
  `t09_operator_audit.csv`; figure: `t09_rho_foundation.png`.
- T10 data: `t10_mie_rayleigh_validation.csv`,
  `t10_mie_rayleigh_summary.csv`; figure: `t10_mie_rayleigh_error.png`.
- T11/T11.1 data: `t11_model_e_convergence.csv`,
  `t11_force_oracle.csv`, `t11_force_decomposition.csv`,
  `t11_1_solver_stability.csv`, `t11_1_high_precision_oracle.csv`;
  figures: `t11_model_e_validation.png`, `t11_1_model_e_stability.png`.
- T12 data/manifest: `t12_sentinel_manifest.csv`,
  `t12_model_e_convergence.csv`, `t12_model_comparison.csv`,
  `t12_threshold_audit.csv`; figure: `t12_model_e_sentinel_audit.png`.
- T12.1 data: `t12_1_extended_convergence.csv`,
  `t12_1_convergence_summary.csv`, `t12_1_resolved_comparison.csv`,
  `t12_1_mechanism_diagnostics.csv`,
  `t12_1_predictor_diagnostics.csv`,
  `t12_1_out_of_fold_predictions.csv`;
  figure: `t12_1_rho1_failure_diagnostics.png`.
- T12.2 data: `t12_2_logo_predictions.csv`, `t12_2_logo_fits.csv`,
  `t12_2_metrics.csv`, `t12_2_safety_audit.csv`,
  `t12_2_final_calibration.csv`, `t12_2_gate.csv`;
  figure: `t12_2_rho1_recalibration.png`.
- T12.3 data: `t12_3_case_influence.csv`, `t12_3_gate.csv`,
  `t12_3_group_bootstrap.csv`, `t12_3_logo_coefficients.csv`,
  `t12_3_metrics.csv`, `t12_3_nested_safety_factors.csv`,
  `t12_3_oof_predictions.csv`, `t12_3_threshold_audit.csv`;
  figure: `t12_3_mechanistic_validity.png`.
- T13 blind/manifest: `t13_holdout_manifest.csv`,
  `t13_frozen_predictions.csv`, `t13_frozen_protocol.csv`;
  revealed: `t13_model_e_convergence.csv`, `t13_forces.csv`,
  `t13_case_summary.csv`, `t13_external_predictions.csv`,
  `t13_metrics.csv`, `t13_threshold_audit.csv`, `t13_gate.csv`;
  figure: `t13_external_validation.png`.
- T14 blind/manifest: `t14_scale_manifest.csv`,
  `t14_frozen_predictions.csv`, `t14_frozen_protocol.csv`,
  `t14_prior_artifact_hashes.csv`; revealed:
  `t14_model_e_convergence.csv`, `t14_forces.csv`,
  `t14_case_summary.csv`, `t14_scale_predictions.csv`,
  `t14_metrics.csv`, `t14_threshold_audit.csv`,
  `t14_matched_scale_pairs.csv`, `t14_performance.csv`,
  `t14_gate.csv`; figure: `t14_scale_out_validation.png`.
- T14.1 blind/manifest: `t14_1_large_n_manifest.csv`,
  `t14_1_local_coupling.csv`, `t14_1_frozen_predictions.csv`,
  `t14_1_frozen_protocol.csv`, `t14_1_prior_artifact_hashes.csv`;
  revealed: `t14_1_model_e_convergence.csv`, `t14_1_forces.csv`,
  `t14_1_case_summary.csv`, `t14_1_large_n_predictions.csv`,
  `t14_1_metrics.csv`, `t14_1_threshold_audit.csv`,
  `t14_1_matched_large_n_pairs.csv`,
  `t14_1_combined_scale_sequence.csv`, `t14_1_performance.csv`,
  `t14_1_gate.csv`; figure: `t14_1_large_n_validation.png`.

All data paths above are under `results/data/`; all figures are under
`results/figures/`. References are
`papers/2014-silva-bruus(2).pdf` and
`papers/Acoustic_Interaction_Force (1)(2).pdf`.
No paper-ready PDF/SVG export, `fit_parameters.csv`, or figure manifest
existed before P0.

## Relevant legacy-column map

| Meaning | Existing columns/tables |
|---|---|
| Case key | `case_id` in T08 and T12–T14.1; T02–T07 use composite geometry/contrast/distance keys. |
| Geometry | `distance_ratio`/`d_min_over_a`; serialized `coordinates_xyz`; long `position_x/y/z` in T13/T14/T14.1. |
| A/B/C/D | `a_x/y`, `b_x/y`, `c_x/y`, `d_x/y`; newer tables use `model_a_*`, `model_d_*`. |
| E channels | Serialized `*_forces_xyz` in convergence tables; long `model_e_total_*`, `model_e_interaction_*`, `model_e_external_scattered_*`, `model_e_scattered_scattered_*` in force tables. |
| Predictors | `eta`, `lambda_max`, `rho_l1`; T14.1 adds `lambda_i` and distribution summaries. |
| Errors | `rms_a_vs_c`, `rms_b_vs_c`, `epsilon_a`, `epsilon_b`, `epsilon_a_e`, and mechanism amplitudes `x_*`. |
| Connected amplitudes | `relative_three_body_sum_amplitude`, `relative_four_body_amplitude`. |
| Convergence | `lmax`, channel `*_successive_change`, `*_change_applicable`, `*_minimum_confirmed_lmax`, `*_confirmed`, `stop_reason`. |
| Quality | `balanced_condition_number`, `balanced_backward_error`, closure/decomposition residuals, `finite`, symmetry/dimension checks. |
| Provenance | `coordinate_sha256`, `local_coupling_sha256`, `source_commit`, prior-artifact hash tables; most legacy files lack schema version and UTC time. |

Adapters must create derived stable keys without rewriting source bytes. New
canonical names are in `docs/PAPER_DATA_CONTRACT.md`.

## Models A–E and editorial aliases

- **A / SB** sums `nodal_pair_force_on_probe` over neighbors
  (`src/acoustic_ms/silva_bruus.py::nodal_pair_force_on_probe`;
  aggregation in `src/acoustic_ms/comparison.py::compare_nodal_force_models`).
- **Historical B** sums isolated Model-C Rayleigh dimers at \(L=1\), not the
  analytical correction or E
  (`src/acoustic_ms/comparison.py::compare_nodal_force_models`).
- **B_E / complete-dimer sum** is an editorial target only; no implementation
  loops over isolated pairs with Model E and sums interaction forces.
- **C / dipolar MS** is global \(L=1\) Rayleigh MS with the T04 observable
  (`src/acoustic_ms/force.py::solve_rayleigh_nodal_interaction_forces`).
- **D / multipolar bridge** uses leading-Rayleigh coefficients
  (`src/acoustic_ms/model_d.py::solve_multipolar_nodal_interaction_forces`;
  `src/acoustic_ms/multipolar_scattering.py::rayleigh_multipolar_scattering_coefficients`).
- **E / complete MS** combines exact Mie, global coupling and complete force
  (`src/acoustic_ms/model_e.py::solve_model_e_nodal`;
  `src/acoustic_ms/mie_multiparticle.py::solve_mie_multiparticle_nodal`;
  `src/acoustic_ms/complete_force.py::complete_radiation_force_from_bsc`).

No legacy name/column is redefined. **B_E** may appear only in new layers.

## Exact E interaction force and implemented metrics

`src/acoustic_ms/model_e.py::solve_model_e_nodal` computes

\[
\mathbf F^E_{\rm int}=\mathcal F[\mathbf b]-\mathcal F[\mathbf a]
=\mathbf F^E_{\rm ext-sc}+\mathbf F^E_{\rm ss},
\]

with API fields `ModelENodalResult.interaction_forces_xyz`,
`external_scattered_forces_xyz`, and
`scattered_scattered_forces_xyz`. The paper reference is interaction, not
total force.

- \(\mathcal R(F)=[N^{-1}\sum_i\lVert F_i\rVert_2^2]^{1/2}\):
  `src/acoustic_ms/metrics.py::rms_vector_magnitude` (2-D) and
  `src/acoustic_ms/model_e_comparison.py::rms_vector_magnitude_xyz` (3-D).
- \(\varepsilon_A^E=\mathcal R(F^A-F^E_{\rm int})/
  \mathcal R(F^E_{\rm int})\), with applicability:
  `src/acoustic_ms/model_e_comparison.py::normalized_rms_error_xyz`.
- \(\rho_1\): `src/acoustic_ms/transferability.py::spectral_radius_l1`,
  independently reconstructed by
  `src/acoustic_ms/rho_foundation.py::dipolar_balanced_coupling_matrix`.
- \(\Lambda_{\max}=|f_1|\max_i\sum_{j\ne i}(a/r_{ij})^3\):
  `src/acoustic_ms/scaling.py::maximum_geometric_coupling`; vector form:
  `src/acoustic_ms/large_n_validation.py::local_geometric_coupling`.
- \(Y_3=\mathcal R(\Phi_\Sigma^{(3)})/\mathcal R(F^C)\) and
  \(Y_4=\mathcal R(\Phi^{(4)})/\mathcal R(F^C)\):
  `scripts/analyze_t06_scaling.py::RESPONSES` and `_augment_t06`.
  These are Model-C amplitudes, not Model-E connected terms.

## Evidence classification under the new protocol

| Work | Classification | Reason |
|---|---|---|
| T01–T05 | `exploratory` | Model construction and first comparisons. |
| T06–T11.1 | `development` | Connected-order/predictor exploration, bridge and E construction. |
| T12/T12.1 | `legacy_validation` | Preregistered/follow-up checks predating P0. |
| T12.2/T12.3 | `development` | Recalibration and M1/M2 selection affected coefficients/margins. |
| T13/T14/T14.1 | `legacy_validation` | Historically blind, but completed before P0; never `confirmatory_new`. |
| P1–P6 | `confirmatory_new` only after preregistration | No case exists yet. |

Predictor/fit/threshold choices were influenced by T06.1, T08, T12.1,
T12.2 and T12.3. T13–T14.1 tested frozen choices without refitting.

## Gaps, compatibility risks, and debt

1. **B_E is absent.** E solves a dimer
   (`src/acoustic_ms/model_e.py::solve_model_e_nodal`) but there is no
   validated complete-dimer sum, cache, convergence propagation, or A/B_E/E
   identity API.
2. **\(\Phi_E^{(3..5)}\) is absent.** Inclusion–exclusion is only C at
   \(N=4,L=1\)
   (`src/acoustic_ms/cluster_expansion.py::decompose_nodal_quartet`) or D
   for \(N\le4\)
   (`src/acoustic_ms/multipolar_expansion.py::decompose_multipolar_cluster`).
   There is no common-confirmed-\(L\) E subset engine or fifth-order term.
3. T13–T14.1 are `legacy_validation`; P1–P6 samples remain TBD.
4. No adapter yet materializes canonical raw/derived/plot-ready data.
5. `docs/CONVENTIONS.md` has a pre-existing interrupted displayed
   definition of \(r_b\). P0 does not correct conventions.
6. Legacy CSVs mix 2-D/3-D names, serialized/long vectors and distance names.
7. Existing scripts have local serializers/styles; P0 does not convert them.
8. Naïve E inclusion–exclusion is exponential; D's current implementation is
   explicitly \(N\le4\)
   (`src/acoustic_ms/multipolar_expansion.py::decompose_multipolar_cluster`).
9. E requires \(L_{\max}\ge2\)
   (`src/acoustic_ms/model_e.py::solve_model_e_nodal`); B_E must propagate
   channel-specific convergence rather than silently share an order.
10. The task-record relocation is unresolved user state and is not staged by
    P0.
