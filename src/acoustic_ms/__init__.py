"""Acoustic multiple-scattering research tools.

T01 provides the original Silva--Bruus pair force; T02 provides the corrected
two-particle analytical benchmark; T03 provides the coupled Rayleigh solver at
Lmax=1; T04 provides the Model C nodal interaction force with scattering
Lmax=1 and local evaluation through ell=2; T05 compares A/B/C trimers; T06 adds the connected N=4 body expansion through four-body order at fixed Lmax=1; T06.1 adds post-processing predictors and log-space scaling diagnostics; T07 adds balanced multipolar Model D and convergence diagnostics.
"""

from .contrasts import dipole_contrast, monopole_contrast
from .comparison import NodalForceModelComparison, compare_nodal_force_models
from .cluster_expansion import NodalQuartetBodyExpansion, decompose_nodal_quartet
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
from .multipolar_solver import MultipolarNodalSolution, solve_multipolar_nodal
from .model_d import MultipolarNodalInteractionResult, NodalModelDComparison, compare_nodal_model_d, solve_multipolar_nodal_interaction_forces
from .multipolar_expansion import MultipolarClusterExpansion, MultipolarConnectedTerm, decompose_multipolar_cluster
from .scattering import rayleigh_scattering_coefficients
from .special import spherical_hankel1
from .solver import RayleighNodalSolution, solve_rayleigh_nodal
from .translation import separation_coefficient, translation_matrix
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
    "dipole_contrast",
    "NodalForceModelComparison",
    "NodalQuartetBodyExpansion",
    "NodalModelDComparison",
    "MultipolarClusterExpansion",
    "MultipolarConnectedTerm",
    "MultipolarNodalInteractionResult",
    "MultipolarNodalSolution",
    "PowerLawFit",
    "angular_errors_degrees",
    "compare_nodal_force_models",
    "compare_nodal_model_d",
    "coupling_eta",
    "fit_power_law",
    "equilateral_trimer",
    "irregular_quartet",
    "linear_quartet",
    "linear_trimer",
    "maximum_geometric_coupling",
    "rms_relative_error",
    "rms_vector_magnitude",
    "scalene_trimer",
    "square_quartet",
    "symmetric_particle_errors",
    "decompose_nodal_quartet",
    "decompose_multipolar_cluster",
    "corrected_nodal_pair_force_magnitude",
    "corrected_nodal_pair_force_on_probe",
    "corrected_nodal_pair_forces",
    "corrected_pair_coefficients",
    "monopole_contrast",
    "nodal_pair_force_magnitude",
    "nodal_pair_force_on_probe",
    "nodal_pair_forces",
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
    "separation_coefficient",
    "solve_rayleigh_nodal",
    "solve_multipolar_nodal",
    "solve_multipolar_nodal_interaction_forces",
    "spherical_hankel1",
    "translation_matrix",
]
