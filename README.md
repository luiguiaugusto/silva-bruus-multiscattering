# Silva--Bruus multiscattering

Python tools for studying acoustic interaction forces between identical spheres
in a pressure nodal plane. T02 includes the Silva--Bruus pairwise force and the
published corrected fifth-order analytical two-particle benchmark; multiple
scattering remains out of scope.

## Setup

```bash
python -m venv .venv
.venv/bin/python -m pip install -e ".[dev,plot]"
.venv/bin/python -m pytest
.venv/bin/python scripts/reproduce_figure_2.py
```

See [the conventions](docs/CONVENTIONS.md) for the time and energy
normalizations used by the code.

The Figure 2 reproduction uses `ka = 0.1`, `radius = 1 m`, and `energy_density = 1 J m^-3`; the reported relative error is independent of this shared SI scale. Its contact-limit sample at `kd = 0.2 = 2ka` is included as the non-overlap boundary, not as a positive surface gap.
