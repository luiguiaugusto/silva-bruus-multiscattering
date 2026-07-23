# Silva--Bruus multiscattering handoff

## T01

### Files created

- Project configuration: `pyproject.toml`, `.gitignore`, `README.md`, `AGENTS.md`, and `TASKS.md`.
- Package: `src/acoustic_ms/__init__.py`, `contrasts.py`, and `silva_bruus.py`.
- Documentation: `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and this file.
- Tests: `tests/test_contrasts.py` and `tests/test_silva_bruus.py`.

### Implemented equations

`monopole_contrast` implements \(f_0=1-\kappa_p/\kappa_0\), and `dipole_contrast` implements \(f_1=2(\rho_p/\rho_0-1)/(2\rho_p/\rho_0+1)\). The nodal pair-force API implements the specified Silva--Bruus expression with \(-\tfrac32[\cos(kd)+kd\sin(kd)]\); the `kd sin(kd)` term therefore has the audited negative sign after multiplication.

### Adopted conventions

The temporal convention is \(e^{-i\omega t}\) and \(E_0=\rho_0|v_0|^2/4\). Positions are 2D SI coordinates in the pressure nodal plane, and \(\widehat{\mathbf d}_{ij}\) points from source to probe. Negative radial force is attractive. The public API validates the non-overlap domain \(d\ge2a\). Full details are in `CONVENTIONS.md`.

### Commands and verification

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/python -m pytest
```

The final run completed successfully: **16 passed in 0.07 s**. It includes an independent centred numerical derivative of the specified potential, asymptotic limits, planar symmetries, action--reaction, zeros, and invalid inputs.

### Limitations and scope

T01 intentionally contains no special functions, T-matrix, translation operator, solver, corrected 2026 two-body formula, notebooks, or plots. `nodal_pair_force_magnitude` is signed along source-to-probe rather than an unsigned norm, preserving physical attraction/repulsion.

### Git diff summary

`git diff --stat` não lista arquivos ainda não rastreados; nesta inicialização, todos os arquivos da T01 são novos e não rastreados. O estado final inclui 12 arquivos de infraestrutura, código, documentação e testes, além de `.gitignore`.

## T02

### Implemented model and equations

`corrected_pair.py` implements Eqs. (30a)--(30d) of the 2026 reference: the coefficients \(A_0\), \(A_2\), and \(D_0\), followed by the corrected signed nodal two-particle force. The direction, energy normalization, and non-overlap validation are unchanged from T01. It is a fifth-order analytical two-particle benchmark, not a multiple-scattering solver.

### Figure 2 reproduction

`scripts/reproduce_figure_2.py` uses public package force functions with `ka=0.1`, 501 samples on `0.2 <= kd <= 0.3`, and `f1 = 0.1, 0.4, 0.8, 1.0`. The error is `100 * abs(corrected - SB) / abs(corrected)`. The contact-limit point `kd=0.2=2ka` is retained as the permitted non-overlap boundary. The script writes the CSV and PNG below `results/`.

### Commands and verification

```bash
.venv/bin/python -m pip install -e ".[dev,plot]"
.venv/bin/python -m pytest -q
.venv/bin/python scripts/reproduce_figure_2.py
.venv/bin/python -m pytest -q -W error
```

The regression suite covers coefficient ratios, all published Figure 2 contact errors, the `kd=0.3` check, limits, physical validation, planar symmetries, and monotonic error curves. No T03 or `N >= 3` implementation was added.

The final verification completed with **32 passed** under `pytest -q -W error`. The reproduced contact-limit errors are 1.252519728707%, 5.160511340274%, 10.798343941865%, and 13.848266387733% for `f1 = 0.1, 0.4, 0.8, 1.0`.

## T03

### Implementation

Added multipole indexing, special functions, cached Gaunt coefficients, target<-source translation, Rayleigh coefficients, nodal incident coefficients, and the dense coupled `Lmax=1` solver. The system is `(I - D_g U)s = D_g a_ext` and returns the residual and 2-norm condition number.

### Verification

The T03 suite includes the direct 3D reexpansion theorem test (including source mode `(1,1)`) with relative error below `1e-9`, analytic one- and two-particle benchmarks, and a structural three-particle test. The deterministic validation script writes `results/data/t03_solver_validation.csv`. No multibody radiation force was implemented.

### T03.1 coverage closure

Added independent public-API Hankel derivative recurrence checks (`ell=1,2,4`, tolerance `2e-13`), explicit odd-`m` harmonic conjugation, Gauss--Legendre/periodic-azimuth harmonic orthonormality, and numerical Gaunt quadratures for three nonzero couplings (including negative azimuthal index and `q=2`, tolerance `3e-12`).

The three-dimensional reexpansion test now verifies strict convergence for `L_test=2,4,6,8,10`; the terminal error is below `1e-9`. The Rayleigh small-`kd` test checks the coupling limit with relative tolerance `3e-6` and the pair/single coefficient ratio with `2e-7`. These are test-only quadratures; the production solver, equations, conventions, and validated CSV remain unchanged. T03 continues to compute scattered-field coefficients only, never multibody force.
Measured T03.1 results: the `L_test=10` reexpansion relative error is `3.6882199807785393e-13` (sequence: `6.008796997e-4`, `3.453822361e-6`, `1.756874298e-8`, `8.288116268e-11`, `3.688219981e-13`); the small-`kd` coupling-limit relative error is `1.999999778e-6`, and the pair/single ratio relative error is `1.052631574e-7`. Final verification: `58 passed` with `-W error`.

## T04: nodal Rayleigh interaction force

### Files and equations

Added `force.py`, `test_force.py`, `validate_t04_force.py`, and `t04_pair_force_validation.csv`. The public API calls the T03 `Lmax=1` solver once, reexpands other particles locally through `ell=2`, excludes self fields, and applies Eq. (27) via `b_{2,-1}` and `b_{2,1}`. It implements interaction force only: no primary off-nodal force and no scattered--scattered terms.

### Validation

The independent Cartesian oracle agrees with production for aligned and oblique pairs; the maximum scalar-reference error in the validation CSV is `4.276e-16`. Contact radial regressions (`a=E0=1`, `ka=0.1`, `d/a=2`) are `-0.011936371917121`, `-0.194729303800953`, `-0.799842697325624`, and `-1.26676999261163` for `f1=0.1,0.4,0.8,1.0`. The suite has 70 tests with warnings promoted to errors.

### Remaining limits

The solver is still `Lmax=1`; `ell=2` is only local force evaluation. T04 intentionally produces no interpreted or published `N>=3` force results, no dynamics, no torque, and no viscosity, streaming, wall, or scattered--scattered terms.
