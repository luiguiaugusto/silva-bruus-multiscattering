"""Matched pairwise baselines and leakage-safe T08 analysis utilities."""

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy.stats import rankdata

from .metrics import rms_vector_magnitude
from .model_d import MultipolarNodalInteractionResult, solve_multipolar_nodal_interaction_forces
from .scaling import PowerLawFit, fit_power_law


@dataclass(frozen=True)
class MatchedPairwiseBaseline:
    """Sum of isolated Model-D dimers at the same multipolar order."""

    forces_xy: np.ndarray
    lmax: int
    pair_count: int
    maximum_residual: float
    maximum_balanced_condition: float
    maximum_raw_condition: float


@dataclass(frozen=True)
class TransferabilityFit:
    """Power-law diagnostics augmented by Spearman rank correlation."""

    power_law: PowerLawFit
    spearman: float


@lru_cache(maxsize=None)
def _cached_dimer(
    distance: float, k: float, radius: float, energy_density: float,
    f0: float, f1: float, lmax: int,
) -> tuple[float, float, float, float]:
    positions = np.array([[-distance / 2, 0.0, 0.0], [distance / 2, 0.0, 0.0]])
    result = solve_multipolar_nodal_interaction_forces(
        positions, k, radius, energy_density, f0, f1, lmax
    )
    return (
        float(result.forces_xy[0, 0]), result.solution.residual_relative,
        result.solution.condition_number, result.solution.physical_condition_number,
    )


def matched_multipolar_pairwise_baseline(
    positions_xyz: object, k: object, radius: object, energy_density: object,
    f0: object, f1: object, lmax: int,
) -> MatchedPairwiseBaseline:
    """Return B_L, the sum of isolated Model-D pair solutions at order L."""
    positions = np.asarray(positions_xyz, dtype=float)
    if positions.ndim != 2 or positions.shape[1:] != (3,) or len(positions) < 2:
        raise ValueError("positions_xyz must have shape (N, 3), N >= 2")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions_xyz must be finite")
    forces = np.zeros((len(positions), 2), dtype=float)
    residuals, balanced, raw = [], [], []
    for first, second in combinations(range(len(positions)), 2):
        displacement = positions[second, :2] - positions[first, :2]
        distance = float(np.linalg.norm(positions[second] - positions[first]))
        if distance == 0.0:
            raise ValueError("particle centers must not coincide")
        radial, residual, condition, raw_condition = _cached_dimer(
            distance, float(k), float(radius), float(energy_density),
            float(f0), float(f1), lmax,
        )
        direction = displacement / distance
        forces[first] += radial * direction
        forces[second] -= radial * direction
        residuals.append(residual); balanced.append(condition); raw.append(raw_condition)
    return MatchedPairwiseBaseline(
        forces, lmax, len(residuals), max(residuals), max(balanced), max(raw)
    )


def spectral_radius_l1(result: MultipolarNodalInteractionResult) -> float:
    """Return the spectral radius of K_b = I - A_b for an L=1 solution."""
    solution = result.solution
    if solution.lmax != 1:
        raise ValueError("spectral_radius_l1 requires an L=1 solution")
    matrix = solution.system_matrix
    if matrix.size == 0:
        return 0.0
    coupling = np.eye(len(matrix), dtype=complex) - matrix
    return float(np.max(np.abs(np.linalg.eigvals(coupling))))


def two_step_converged(successive_differences: object, tolerance: float = 1e-3) -> bool:
    """Require the last two successive differences to satisfy the tolerance."""
    values = np.asarray(successive_differences, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        return False
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and non-negative")
    return bool(np.all(values[-2:] <= tolerance))


def normalized_rms_difference(first: object, second: object, reference: object) -> tuple[float, bool]:
    """Return an RMS difference normalized by a scale-aware reference."""
    first = np.asarray(first, dtype=float); second = np.asarray(second, dtype=float)
    reference = np.asarray(reference, dtype=float)
    denominator = rms_vector_magnitude(reference)
    scale = max(rms_vector_magnitude(first), rms_vector_magnitude(second), denominator)
    tolerance = 128 * np.finfo(float).eps * scale
    if denominator <= tolerance:
        return 0.0, False
    return rms_vector_magnitude(first - second) / denominator, True


def fit_transferability_power_law(x: object, y: object) -> TransferabilityFit:
    """Fit a power law and compute Spearman correlation without filtering."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    fit = fit_power_law(x, y)
    rx, ry = rankdata(x), rankdata(y)
    spearman = float(np.corrcoef(rx, ry)[0, 1])
    return TransferabilityFit(fit, spearman)


def select_predictor_by_group_cv(
    calibration_rows: list[dict], predictors: tuple[str, ...] = ("eta", "lambda_max", "rho_l1"),
) -> tuple[str, dict[str, float]]:
    """Select using calibration-only leave-(N,family)-out log RMSE."""
    if any(int(row["particle_count"]) > 4 for row in calibration_rows):
        raise ValueError("cross-validation rows must contain only N <= 4")
    groups = sorted({(int(row["particle_count"]), row["family"]) for row in calibration_rows})
    scores = {}
    for predictor in predictors:
        residuals = []
        for group in groups:
            train = [r for r in calibration_rows if (int(r["particle_count"]), r["family"]) != group]
            test = [r for r in calibration_rows if (int(r["particle_count"]), r["family"]) == group]
            fit = fit_power_law(
                [float(r[predictor]) for r in train], [float(r["epsilon_a"]) for r in train]
            )
            for row in test:
                predicted = fit.prefactor * float(row[predictor]) ** fit.exponent
                residuals.append(np.log(float(row["epsilon_a"]) / predicted))
        scores[predictor] = float(np.sqrt(np.mean(np.square(residuals))))
    priority = {"lambda_max": 0, "rho_l1": 1, "eta": 2}
    best_score = min(scores.values())
    tied = [name for name, score in scores.items() if score <= best_score + 1e-12]
    return min(tied, key=lambda name: priority[name]), scores


def conservative_threshold(
    calibration_rows: list[dict], predictor: str, tolerance: float, minimum_count: int = 8,
) -> tuple[float, bool, int]:
    """Return the largest sampled predictor whose full prefix is safe."""
    if any(int(row["particle_count"]) > 4 for row in calibration_rows):
        raise ValueError("threshold calibration rows must contain only N <= 4")
    ordered = sorted(calibration_rows, key=lambda row: float(row[predictor]))
    available = []
    for threshold in sorted({float(row[predictor]) for row in ordered}):
        prefix = [row for row in ordered if float(row[predictor]) <= threshold]
        if not all(float(row["epsilon_a"]) <= tolerance for row in prefix):
            break
        if len(prefix) >= minimum_count:
            available.append((threshold, len(prefix)))
    if not available:
        return 0.0, False, 0
    threshold, count = available[-1]
    return threshold, True, count
