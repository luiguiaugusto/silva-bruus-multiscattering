# Silva--Bruus multiscattering

Python tools for studying acoustic interaction forces between identical spheres
in a pressure nodal plane. T05 reaches canonical N=3 trimers: Model A is Silva--Bruus pairwise, Model B sums isolated T04 pair solves, and Model C is global Rayleigh multiple scattering at Lmax=1. Model D, higher multipole orders, and N>3 remain out of scope.

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

## T03 solver

The Rayleigh `Lmax=1` field solver is available through `solve_rayleigh_nodal`. It computes four multipole coefficients per particle and resums rescattering in a dense linear solve; it does not compute a radiation force. Generate validation data with:

```bash
.venv/bin/python scripts/validate_t03_solver.py
```

## T04 interaction force

T04 adds `solve_rayleigh_nodal_interaction_forces`, which returns the coupled T03 field solution, local scattered coefficients through `ell=2`, and nodal-plane interaction forces in newtons for one or two particles. It is the Rayleigh Model C cross-term force; it does not add off-nodal total force or multibody force results.

```bash
.venv/bin/python scripts/validate_t04_force.py
```
