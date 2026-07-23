# T01 handoff

## Files created

- Project configuration: `pyproject.toml`, `.gitignore`, `README.md`, `AGENTS.md`, and `TASKS.md`.
- Package: `src/acoustic_ms/__init__.py`, `contrasts.py`, and `silva_bruus.py`.
- Documentation: `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and this file.
- Tests: `tests/test_contrasts.py` and `tests/test_silva_bruus.py`.

## Implemented equations

`monopole_contrast` implements \(f_0=1-\kappa_p/\kappa_0\), and `dipole_contrast` implements \(f_1=2(\rho_p/\rho_0-1)/(2\rho_p/\rho_0+1)\). The nodal pair-force API implements the specified Silva--Bruus expression with \(-\tfrac32[\cos(kd)+kd\sin(kd)]\); the `kd sin(kd)` term therefore has the audited negative sign after multiplication.

## Adopted conventions

The temporal convention is \(e^{-i\omega t}\) and \(E_0=\rho_0|v_0|^2/4\). Positions are 2D SI coordinates in the pressure nodal plane, and \(\widehat{\mathbf d}_{ij}\) points from source to probe. Negative radial force is attractive. The public API validates the non-overlap domain \(d\ge2a\). Full details are in `CONVENTIONS.md`.

## Commands and verification

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/python -m pytest
```

The final run completed successfully: **16 passed in 0.07 s**. It includes an independent centred numerical derivative of the specified potential, asymptotic limits, planar symmetries, action--reaction, zeros, and invalid inputs.

## Limitations and scope

T01 intentionally contains no special functions, T-matrix, translation operator, solver, corrected 2026 two-body formula, notebooks, or plots. `nodal_pair_force_magnitude` is signed along source-to-probe rather than an unsigned norm, preserving physical attraction/repulsion.

## Git diff summary

`git diff --stat` não lista arquivos ainda não rastreados; nesta inicialização, todos os arquivos da T01 são novos e não rastreados. O estado final inclui 12 arquivos de infraestrutura, código, documentação e testes, além de `.gitignore`.

# T02 handoff

## Implemented model and equations

`corrected_pair.py` implements Eqs. (30a)--(30d) of the 2026 reference: the coefficients \(A_0\), \(A_2\), and \(D_0\), followed by the corrected signed nodal two-particle force. The direction, energy normalization, and non-overlap validation are unchanged from T01. It is a fifth-order analytical two-particle benchmark, not a multiple-scattering solver.

## Figure 2 reproduction

`scripts/reproduce_figure_2.py` uses public package force functions with `ka=0.1`, 501 samples on `0.2 <= kd <= 0.3`, and `f1 = 0.1, 0.4, 0.8, 1.0`. The error is `100 * abs(corrected - SB) / abs(corrected)`. The contact-limit point `kd=0.2=2ka` is retained as the permitted non-overlap boundary. The script writes the CSV and PNG below `results/`.

## Commands and verification

```bash
.venv/bin/python -m pip install -e ".[dev,plot]"
.venv/bin/python -m pytest -q
.venv/bin/python scripts/reproduce_figure_2.py
.venv/bin/python -m pytest -q -W error
```

The regression suite covers coefficient ratios, all published Figure 2 contact errors, the `kd=0.3` check, limits, physical validation, planar symmetries, and monotonic error curves. No T03 or `N >= 3` implementation was added.

The final verification completed with **32 passed** under `pytest -q -W error`. The reproduced contact-limit errors are 1.252519728707%, 5.160511340274%, 10.798343941865%, and 13.848266387733% for `f1 = 0.1, 0.4, 0.8, 1.0`.
