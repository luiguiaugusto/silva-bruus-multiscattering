"""Acoustic multiple-scattering research tools.

T01 provides the original nodal Silva--Bruus pair force; T02 provides the
corrected two-particle benchmark; T03 provides an Lmax=1 coupled Rayleigh
field-coefficient solver.
"""

from .contrasts import dipole_contrast, monopole_contrast
from .incident import nodal_standing_wave_coefficients
from .multipoles import mode_count, mode_from_index, mode_index, modes
from .scattering import rayleigh_scattering_coefficients
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
    "corrected_nodal_pair_force_magnitude",
    "corrected_nodal_pair_force_on_probe",
    "corrected_nodal_pair_forces",
    "corrected_pair_coefficients",
    "monopole_contrast",
    "nodal_pair_force_magnitude",
    "nodal_pair_force_on_probe",
    "nodal_pair_forces",
    "RayleighNodalSolution",
    "mode_count",
    "mode_from_index",
    "mode_index",
    "modes",
    "nodal_standing_wave_coefficients",
    "rayleigh_scattering_coefficients",
    "separation_coefficient",
    "solve_rayleigh_nodal",
    "translation_matrix",
]
