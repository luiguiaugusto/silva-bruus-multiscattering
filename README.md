# Silva--Bruus multiscattering

Initial Python infrastructure for studying acoustic interaction forces between
spheres.  T01 implements only the Silva--Bruus pairwise force for identical
spheres in a pressure nodal plane; multiple scattering is deliberately out of
scope.

## Setup

```bash
python -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/python -m pytest
```

See [the conventions](docs/CONVENTIONS.md) for the time and energy
normalizations used by the code.
