"""Analytical dipolar coupling operator underlying the T08 ``rho_1`` metric."""

from dataclasses import dataclass

import numpy as np

from .solver import _validate_positions


@dataclass(frozen=True)
class DipolarCouplingDiagnostics:
    """Spectral and norm diagnostics for the exact and near-field operators."""

    spectral_radius: float
    spectral_norm: float
    infinity_norm: float
    near_field_spectral_radius: float
    normalized_commutator: float


def _validated_coupling_scalars(
    distance: object, k: object, radius: object, f1: object
) -> tuple[float, float, float, float]:
    try:
        distance, k, radius, f1 = (
            float(value) for value in (distance, k, radius, f1)
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("distance, k, radius, and f1 must be real scalars") from exc
    if not np.all(np.isfinite((distance, k, radius, f1))):
        raise ValueError("distance, k, radius, and f1 must be finite")
    if distance <= 0.0 or k <= 0.0 or radius <= 0.0:
        raise ValueError("distance, k, and radius must be positive")
    if k * radius > 0.1:
        raise ValueError("the Rayleigh regime requires ka <= 0.1")
    if distance < 2.0 * radius:
        raise ValueError("distance must satisfy separation >= 2 * radius")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")
    return distance, k, radius, f1


def dipolar_coupling_entry(
    distance: object, k: object, radius: object, f1: object
) -> complex:
    r"""Return one off-diagonal entry of the balanced nodal ``L=1`` operator.

    For identical spheres in the nodal plane,

    .. math::

       K_{ij} = \frac{f_1}{2}\left(\frac{a}{r_{ij}}\right)^3
                e^{ikr_{ij}}(1-ikr_{ij}).

    This is exact within the leading-Rayleigh dipole model; it is not a
    near-field expansion and it is not an exact finite-frequency T-matrix.
    """
    distance, k, radius, f1 = _validated_coupling_scalars(
        distance, k, radius, f1
    )
    phase = k * distance
    return complex(
        0.5
        * f1
        * (radius / distance) ** 3
        * np.exp(1j * phase)
        * (1.0 - 1j * phase)
    )


def near_field_dipolar_coupling_entry(
    distance: object, k: object, radius: object, f1: object
) -> float:
    r"""Return the quasistatic term ``f1/2 * (a/r)**3`` of the dipole coupling."""
    distance, _, radius, f1 = _validated_coupling_scalars(
        distance, k, radius, f1
    )
    return float(0.5 * f1 * (radius / distance) ** 3)


def _validated_matrix_inputs(
    positions_xyz: object, k: object, radius: object, f1: object
) -> tuple[np.ndarray, float, float, float]:
    try:
        k, radius, f1 = (float(value) for value in (k, radius, f1))
    except (TypeError, ValueError) as exc:
        raise ValueError("k, radius, and f1 must be real scalars") from exc
    if not np.all(np.isfinite((k, radius, f1))) or k <= 0.0 or radius <= 0.0:
        raise ValueError("require finite k > 0, radius > 0, and finite f1")
    if k * radius > 0.1:
        raise ValueError("the Rayleigh regime requires ka <= 0.1")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")
    return _validate_positions(positions_xyz, radius), k, radius, f1


def dipolar_balanced_coupling_matrix(
    positions_xyz: object, k: object, radius: object, f1: object
) -> np.ndarray:
    """Return the exact balanced nodal ``L=1`` rescattering matrix."""
    positions, k, radius, f1 = _validated_matrix_inputs(
        positions_xyz, k, radius, f1
    )
    coupling = np.zeros((len(positions), len(positions)), dtype=complex)
    for target in range(len(positions)):
        for source in range(target + 1, len(positions)):
            distance = float(np.linalg.norm(positions[source] - positions[target]))
            entry = dipolar_coupling_entry(distance, k, radius, f1)
            coupling[target, source] = entry
            coupling[source, target] = entry
    return coupling


def near_field_dipolar_coupling_matrix(
    positions_xyz: object, k: object, radius: object, f1: object
) -> np.ndarray:
    """Return the real symmetric quasistatic approximation to the ``L=1`` matrix."""
    positions, k, radius, f1 = _validated_matrix_inputs(
        positions_xyz, k, radius, f1
    )
    coupling = np.zeros((len(positions), len(positions)), dtype=float)
    for target in range(len(positions)):
        for source in range(target + 1, len(positions)):
            distance = float(np.linalg.norm(positions[source] - positions[target]))
            entry = near_field_dipolar_coupling_entry(distance, k, radius, f1)
            coupling[target, source] = entry
            coupling[source, target] = entry
    return coupling


def dipolar_coupling_diagnostics(
    positions_xyz: object, k: object, radius: object, f1: object
) -> DipolarCouplingDiagnostics:
    """Return diagnostics that delimit what the spectral radius can certify."""
    coupling = dipolar_balanced_coupling_matrix(positions_xyz, k, radius, f1)
    near_field = near_field_dipolar_coupling_matrix(
        positions_xyz, k, radius, f1
    )
    eigenvalues = np.linalg.eigvals(coupling)
    near_field_eigenvalues = np.linalg.eigvalsh(near_field)
    frobenius_squared = float(np.linalg.norm(coupling, "fro") ** 2)
    if frobenius_squared == 0.0:
        nonnormality = 0.0
    else:
        commutator = (
            coupling.conj().T @ coupling
            - coupling @ coupling.conj().T
        )
        nonnormality = float(
            np.linalg.norm(commutator, "fro") / frobenius_squared
        )
    return DipolarCouplingDiagnostics(
        spectral_radius=float(np.max(np.abs(eigenvalues))),
        spectral_norm=float(np.linalg.norm(coupling, 2)),
        infinity_norm=float(np.linalg.norm(coupling, np.inf)),
        near_field_spectral_radius=float(
            np.max(np.abs(near_field_eigenvalues))
        ),
        normalized_commutator=nonnormality,
    )


def neumann_partial_solutions(
    coupling_matrix: object, source_vector: object, maximum_order: int
) -> np.ndarray:
    r"""Return ``sum_{p=0}^P K**p b`` for every ``0 <= P <= maximum_order``."""
    coupling = np.asarray(coupling_matrix, dtype=complex)
    source = np.asarray(source_vector, dtype=complex)
    if coupling.ndim != 2 or coupling.shape[0] != coupling.shape[1]:
        raise ValueError("coupling_matrix must be square")
    if source.shape != (coupling.shape[0],):
        raise ValueError("source_vector must match the coupling dimension")
    if not np.all(np.isfinite(coupling)) or not np.all(np.isfinite(source)):
        raise ValueError("coupling_matrix and source_vector must be finite")
    if (
        not isinstance(maximum_order, int)
        or isinstance(maximum_order, bool)
        or maximum_order < 0
    ):
        raise ValueError("maximum_order must be a non-negative integer")
    term = source.copy()
    partial = source.copy()
    solutions = [partial.copy()]
    for _ in range(maximum_order):
        term = coupling @ term
        partial = partial + term
        solutions.append(partial.copy())
    return np.asarray(solutions)
