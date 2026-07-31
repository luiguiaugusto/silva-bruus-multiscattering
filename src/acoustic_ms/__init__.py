"""Acoustic multiple-scattering research tools.

T01 provides the original Silva--Bruus pair force; T02 provides the corrected
two-particle analytical benchmark; T03 provides the coupled Rayleigh solver at
Lmax=1; T04 provides the Model C nodal interaction force with scattering
Lmax=1 and local evaluation through ell=2; T05 compares A/B/C trimers; T06 adds the connected N=4 body expansion; T07 adds multipolar Model D; T08--T09 audit rho_1 transferability and its operator; T10 adds exact isolated-sphere Mie coefficients; T11 adds globally coupled exact-Mie Model E and its complete multipolar radiation force without changing Models A--D.
T11.1 stabilizes Model E with a square-root-balanced linear solve; T12 adds
three-dimensional sentinel comparison metrics without recalibrating rho_1.
"""

from .contrasts import dipole_contrast, monopole_contrast
from .comparison import NodalForceModelComparison, compare_nodal_force_models
from .cluster_expansion import NodalQuartetBodyExpansion, decompose_nodal_quartet
from .cluster_families import TransferabilityConfiguration, cluster_family, compact_cluster, enumerate_transferability_configurations, irregular_cluster, linear_cluster
from .geometries import (
    equilateral_trimer, irregular_quartet, linear_quartet, linear_trimer,
    scalene_trimer, square_quartet,
)
from .metrics import angular_errors_degrees, rms_relative_error, rms_vector_magnitude, symmetric_particle_errors
from .scaling import PowerLawFit, coupling_eta, fit_power_law, maximum_geometric_coupling
from .incident import nodal_standing_wave_coefficients
from .force import RayleighNodalInteractionResult, solve_rayleigh_nodal_interaction_forces
from .multipoles import mode_count, mode_from_index, mode_index, modes
from .multipolar_scattering import rayleigh_multipolar_scattering_coefficients
from .mie_scattering import (
    fluid_sphere_mie_scattering_coefficients,
    material_ratios_from_contrasts,
    mie_scattering_coefficients_from_contrasts,
    rigid_sphere_scattering_coefficients,
)
from .mie_multiparticle import MieMultiparticleSolution, solve_mie_multiparticle_nodal
from .complete_force import complete_radiation_force_from_bsc
from .model_e import ModelENodalResult, solve_model_e_nodal
from .model_e_comparison import (
    ModelEForceComparison, compare_model_e_forces, normalized_rms_error_xyz,
    rms_vector_magnitude_xyz, symmetric_rms_error_xyz,
)
from .rho1_model_e_diagnostics import (
    ApplicableScalar, ConvergenceTailDiagnostics, LogLinearFit,
    MechanismDiagnostics, OutOfFoldMetrics, convergence_tail_diagnostics,
    fit_log_linear, leave_group_out_folds, mechanism_diagnostics,
    out_of_fold_metrics, spearman_correlation, vector_field_amplitude_ratio,
    vector_field_cosine, vector_field_inner_product, vector_field_projection,
)
from .rho1_model_e_recalibration import (
    BootstrapCalibration, ConfirmatoryMetrics, GateCriterion, LogoFoldFit,
    LogoPrediction, SafetyAudit, SafetyClassification,
    classify_logo_safety, confirmatory_metrics, evaluate_recalibration_gate,
    grouped_bootstrap_calibration, logo_power_law_predictions,
    power_law_threshold,
)
from .mechanistic_validity import (
    MechanisticGateCriterion, MechanisticOofPrediction,
    MechanisticPowerLawFit, MultiplicativeMetrics, NestedLogoFold,
    ThresholdAudit, audit_safety_thresholds, evaluate_mechanistic_gate,
    fit_mechanistic_power_law, fixed_baseline_nested_predictions,
    multiplicative_metrics, nested_logo_predictions,
    predict_mechanistic_power_law,
)
from .multipolar_solver import MultipolarNodalSolution, solve_multipolar_nodal
from .model_d import MultipolarNodalInteractionResult, NodalModelDComparison, compare_nodal_model_d, solve_multipolar_nodal_interaction_forces
from .multipolar_expansion import MultipolarClusterExpansion, MultipolarConnectedTerm, decompose_multipolar_cluster
from .rho_foundation import DipolarCouplingDiagnostics, dipolar_balanced_coupling_matrix, dipolar_coupling_diagnostics, dipolar_coupling_entry, near_field_dipolar_coupling_entry, near_field_dipolar_coupling_matrix, neumann_partial_solutions
from .scattering import rayleigh_scattering_coefficients
from .special import spherical_hankel1
from .solver import RayleighNodalSolution, solve_rayleigh_nodal
from .translation import separation_coefficient, translation_matrix
from .transferability import MatchedPairwiseBaseline, TransferabilityFit, conservative_threshold, fit_transferability_power_law, matched_multipolar_pairwise_baseline, normalized_rms_difference, select_predictor_by_group_cv, spectral_radius_l1, two_step_converged
from .corrected_pair import (
    corrected_nodal_pair_force_magnitude,
    corrected_nodal_pair_force_on_probe,
    corrected_nodal_pair_forces,
    corrected_pair_coefficients,
)
from .silva_bruus import (
    nodal_pair_force_magnitude,
    nodal_pair_force_on_probe,
    nodal_pair_forces,
)

__all__ = [
    "ApplicableScalar",
    "ConvergenceTailDiagnostics",
    "LogLinearFit",
    "MechanismDiagnostics",
    "OutOfFoldMetrics",
    "dipole_contrast",
    "NodalForceModelComparison",
    "NodalQuartetBodyExpansion",
    "NodalModelDComparison",
    "MatchedPairwiseBaseline",
    "TransferabilityConfiguration",
    "TransferabilityFit",
    "DipolarCouplingDiagnostics",
    "MultipolarClusterExpansion",
    "MultipolarConnectedTerm",
    "MultipolarNodalInteractionResult",
    "MultipolarNodalSolution",
    "MieMultiparticleSolution",
    "ModelENodalResult",
    "ModelEForceComparison",
    "PowerLawFit",
    "angular_errors_degrees",
    "compare_nodal_force_models",
    "compare_nodal_model_d",
    "compare_model_e_forces",
    "cluster_family",
    "compact_cluster",
    "conservative_threshold",
    "coupling_eta",
    "fit_power_law",
    "fit_transferability_power_law",
    "equilateral_trimer",
    "irregular_quartet",
    "irregular_cluster",
    "linear_quartet",
    "linear_trimer",
    "linear_cluster",
    "maximum_geometric_coupling",
    "material_ratios_from_contrasts",
    "matched_multipolar_pairwise_baseline",
    "normalized_rms_difference",
    "rms_relative_error",
    "rms_vector_magnitude",
    "rms_vector_magnitude_xyz",
    "scalene_trimer",
    "square_quartet",
    "symmetric_particle_errors",
    "symmetric_rms_error_xyz",
    "decompose_nodal_quartet",
    "decompose_multipolar_cluster",
    "dipolar_balanced_coupling_matrix",
    "dipolar_coupling_diagnostics",
    "dipolar_coupling_entry",
    "enumerate_transferability_configurations",
    "corrected_nodal_pair_force_magnitude",
    "corrected_nodal_pair_force_on_probe",
    "corrected_nodal_pair_forces",
    "corrected_pair_coefficients",
    "monopole_contrast",
    "near_field_dipolar_coupling_entry",
    "near_field_dipolar_coupling_matrix",
    "neumann_partial_solutions",
    "nodal_pair_force_magnitude",
    "nodal_pair_force_on_probe",
    "nodal_pair_forces",
    "normalized_rms_error_xyz",
    "RayleighNodalSolution",
    "RayleighNodalInteractionResult",
    "solve_rayleigh_nodal_interaction_forces",
    "mode_count",
    "mode_from_index",
    "mode_index",
    "modes",
    "nodal_standing_wave_coefficients",
    "rayleigh_scattering_coefficients",
    "rayleigh_multipolar_scattering_coefficients",
    "fluid_sphere_mie_scattering_coefficients",
    "mie_scattering_coefficients_from_contrasts",
    "complete_radiation_force_from_bsc",
    "rigid_sphere_scattering_coefficients",
    "separation_coefficient",
    "select_predictor_by_group_cv",
    "solve_rayleigh_nodal",
    "solve_multipolar_nodal",
    "solve_multipolar_nodal_interaction_forces",
    "solve_mie_multiparticle_nodal",
    "solve_model_e_nodal",
    "spherical_hankel1",
    "spectral_radius_l1",
    "two_step_converged",
    "translation_matrix",
    "convergence_tail_diagnostics",
    "fit_log_linear",
    "leave_group_out_folds",
    "mechanism_diagnostics",
    "out_of_fold_metrics",
    "spearman_correlation",
    "vector_field_amplitude_ratio",
    "vector_field_cosine",
    "vector_field_inner_product",
    "vector_field_projection",
    "BootstrapCalibration",
    "ConfirmatoryMetrics",
    "GateCriterion",
    "LogoFoldFit",
    "LogoPrediction",
    "SafetyAudit",
    "SafetyClassification",
    "classify_logo_safety",
    "confirmatory_metrics",
    "evaluate_recalibration_gate",
    "grouped_bootstrap_calibration",
    "logo_power_law_predictions",
    "power_law_threshold",
    "MechanisticGateCriterion",
    "MechanisticOofPrediction",
    "MechanisticPowerLawFit",
    "MultiplicativeMetrics",
    "NestedLogoFold",
    "ThresholdAudit",
    "audit_safety_thresholds",
    "evaluate_mechanistic_gate",
    "fit_mechanistic_power_law",
    "fixed_baseline_nested_predictions",
    "multiplicative_metrics",
    "nested_logo_predictions",
    "predict_mechanistic_power_law",
]
