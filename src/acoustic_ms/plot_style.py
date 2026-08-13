"""Central, reversible diagnostic plotting style for paper preparation."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal


MM_PER_INCH = 25.4
ONE_COLUMN_MM = 89.0
TWO_COLUMN_MM = 183.0

# Okabe--Ito-derived ordering.  Marker shape must also encode series identity.
ACCESSIBLE_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#56B4E9",
    "#CC79A7",
    "#000000",
)
REDUNDANT_MARKERS = ("o", "s", "^", "D", "v", "P", "X")


def millimetres_to_inches(value: float) -> float:
    """Convert a positive physical length in millimetres to inches."""

    result = float(value)
    if result <= 0.0:
        raise ValueError("millimetre length must be positive")
    return result / MM_PER_INCH


def figure_size(
    columns: Literal["one", "two"] = "one",
    *,
    height_mm: float | None = None,
    aspect_ratio: float = 0.68,
) -> tuple[float, float]:
    """Return final-size figure dimensions in inches."""

    if columns not in {"one", "two"}:
        raise ValueError("columns must be 'one' or 'two'")
    width_mm = ONE_COLUMN_MM if columns == "one" else TWO_COLUMN_MM
    if height_mm is None:
        if aspect_ratio <= 0.0:
            raise ValueError("aspect_ratio must be positive")
        height_mm = width_mm * float(aspect_ratio)
    return millimetres_to_inches(width_mm), millimetres_to_inches(height_mm)


def diagnostic_style() -> dict[str, object]:
    """Return the diagnostic-v1 rc mapping without changing global state."""

    from cycler import cycler

    return {
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9.0,
        "axes.titlesize": 9.0,
        "axes.labelsize": 9.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.0,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "axes.prop_cycle": cycler(color=ACCESSIBLE_COLORS)
        + cycler(marker=REDUNDANT_MARKERS),
        "svg.hashsalt": "acoustic-ms-diagnostic-v1",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }


def diagnostic_rc_context():
    """Return a ``matplotlib.rc_context`` for the reversible style mapping."""

    import matplotlib as mpl

    return mpl.rc_context(diagnostic_style())


def save_diagnostic_figure(
    figure: object,
    output_stem: str | Path,
    *,
    formats: Iterable[Literal["pdf", "svg", "png"]] = ("pdf", "svg", "png"),
    dpi: int = 300,
) -> tuple[Path, ...]:
    """Export a figure in the three supported editorial interchange formats."""

    stem = Path(output_stem)
    if stem.suffix:
        raise ValueError("output_stem must not include a file extension")
    requested = tuple(formats)
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("formats must be a non-empty sequence without duplicates")
    unsupported = set(requested) - {"pdf", "svg", "png"}
    if unsupported:
        raise ValueError(f"unsupported formats: {sorted(unsupported)}")
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    metadata = {
        "pdf": {"Creator": "acoustic-ms diagnostic-v1", "CreationDate": None, "ModDate": None},
        "svg": {"Creator": "acoustic-ms diagnostic-v1", "Date": None},
        "png": {"Software": "acoustic-ms diagnostic-v1"},
    }
    for file_format in requested:
        path = stem.with_suffix(f".{file_format}")
        figure.savefig(path, format=file_format, dpi=dpi, metadata=metadata[file_format])
        paths.append(path)
    return tuple(paths)
