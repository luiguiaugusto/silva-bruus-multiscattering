"""Multipole-mode ordering utilities."""


def mode_count(lmax: int) -> int:
    """Return the number of modes through angular order ``lmax``."""
    if not isinstance(lmax, int) or lmax < 0:
        raise ValueError("lmax must be a non-negative integer")
    return (lmax + 1) ** 2


def mode_index(ell: int, m: int) -> int:
    """Return ``ell**2 + ell + m`` in the project multipole ordering."""
    if not isinstance(ell, int) or not isinstance(m, int) or ell < 0 or abs(m) > ell:
        raise ValueError("ell must be non-negative and abs(m) must not exceed ell")
    return ell * ell + ell + m


def mode_from_index(index: int) -> tuple[int, int]:
    """Return ``(ell, m)`` from a non-negative mode index."""
    if not isinstance(index, int) or index < 0:
        raise ValueError("index must be a non-negative integer")
    ell = int(index**0.5)
    return ell, index - ell * ell - ell


def modes(lmax: int) -> tuple[tuple[int, int], ...]:
    """Return modes ordered by increasing ``ell`` then ``m=-ell,...,ell``."""
    return tuple((ell, m) for ell in range(lmax + 1) for m in range(-ell, ell + 1))
