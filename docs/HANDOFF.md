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

The independent Cartesian oracle agrees with production for aligned and oblique pairs; the maximum scalar-reference error in the validation CSV is `4.276e-16`. Contact radial regressions (`a=E0=1`, `ka=0.1`, `d/a=2`) are `-0.011936371917121`, `-0.194729303800953`, `-0.799842697325624`, and `-1.26676999261163` for `f1=0.1,0.4,0.8,1.0`. The T04.1 final suite has 85 tests with warnings promoted to errors.

### Remaining limits

The solver is still `Lmax=1`; `ell=2` is only local force evaluation. T04 intentionally produces no interpreted or published `N>=3` force results, no dynamics, no torque, and no viscosity, streaming, wall, or scattered--scattered terms.

### T04.1 coverage closure

Added explicit particle-order permutation checks for an oblique pair, particle by particle for forces, local coefficients, and T03 coefficients (`rtol=2e-12`, `atol=2e-14`). For the planar aligned pair, the six forbidden local modes `(0,0)`, `(1,-1)`, `(1,1)`, `(2,-2)`, `(2,0)`, `(2,2)` are bounded below `2e-13`, while the allowed `b_2,-1=-b_2,1` relation remains tested.

The public force API now has direct tests rejecting `energy_density=NaN,+inf,-inf`, `f1=-2.01`, `lmax=0,2`, and invalid position shapes. Existing negative-energy coverage is preserved. Final commands: editable install, both validation scripts, two T04 regenerations, `pytest -q`, `pytest -q -W error`, `git diff --check`, and the no-scientific-module diff check.

Measured final metrics: Cartesian-oracle relative error `6.522504517e-16`; scalar-pair-reference error `4.276e-16`; weak-contrast Silva--Bruus recovery error `6.373752577e-9`; large-separation recovery error `3.890732182e-8`; maximum validation residual `1.498342041e-16`; maximum condition number `1.285874345`; and measured pair action--reaction violation `4.163336342e-17`. The contact regressions remain `-0.011936371917121`, `-0.194729303800953`, `-0.799842697325624`, and `-1.26676999261163`.

Final count: **85 passed** with warnings promoted to errors. Both validation CSVs were regenerated byte-identically: T03 SHA-256 `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`; T04 SHA-256 `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`. T04.1 changed no scientific module, force equation, solver, or CSV bytes. The solver remains `Lmax=1`; `ell=2` remains local evaluation only; no `N>=3` force result was produced.

## T05: first N=3 multibody comparison

### Implementation

Added canonical linear, equilateral, and fixed-shape scalene trimer geometries; Model A/B/C comparison; exact vector decomposition `C-A=(B-A)+(C-B)`; and vector error metrics. B sums T04 solves of each isolated unordered pair, while C makes one global 12-coefficient T03/T04 solve.

### Regressions and artifacts

With `a=E0=1`, `ka=0.1`, `f0=0`, `f1=0.8`, the C RMS errors `(A,C)/(B,C)` are chain `0.0831973786371517/0.0430867438601629`, equilateral `0.0882674934567006/0.0461714463617981`, and scalene `0.0578621870496697/0.0270655878714466`. Generated artifacts are `t05_trimer_regression.csv`, `t05_trimer_sweep.csv`, and `t05_trimer_model_errors.png`.

### Scope

T05 publishes only N=3 canonical comparisons, keeps scattering `Lmax=1` and local evaluation `ell=2`, and adds no Model D, higher multipole order, or N>3 result.
Final validation: **89 passed** with warnings promoted to errors. Sweep maximum residual is `4.577510638855876e-16` and maximum condition number is `1.402967767286417`. Deterministic artifact SHA-256 hashes: regression `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`, sweep `a63e38ab625c518cf40c248209abd5371744f85d280cb0ce603099c6805eaca8`, figure `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`. Approved T03/T04 hashes remained `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a` and `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`.


## T05 / T05.1 / T05.1a / T05.1b / T05.1c audit record

T05.1 files, from the audited Git diff: `README.md`, `TAREFA_T05_1_FECHAMENTO_METRICAS_TESTES_DOCUMENTACAO.md`, `TASKS.md`, `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, `docs/HANDOFF.md`, `results/data/t05_trimer_sweep.csv`, `scripts/validate_t05_trimers.py`, `src/acoustic_ms/__init__.py`, `src/acoustic_ms/metrics.py`, and `tests/test_multibody.py`.

The canonical geometry parameters are \(a=1\), \(ka=0.1\), \(E_0=1\), \(f_0=0\), and \(f_1\in\{0.1,0.4,0.8,1.0\}\). The linear chain is \((-d,0,0),(0,0,0),(d,0,0)\); the equilateral geometry is \((-d/2,-\sqrt{3}d/6,0),(d/2,-\sqrt{3}d/6,0),(0,\sqrt{3}d/3,0)\); the scalene geometry is the centered form of \((0,0,0),(d,0,0),(3d/11,12d/11,0)\). Here \(d_{\min}/a\) is the geometry parameter. The sweep is 3 geometries × 4 values of \(f_1\) × 160 equally spaced values in \([2.1,10.0]\): 1,920 configurations. The regression CSV has nine rows.

| Geometry | Particle | A `(x,y)` | B `(x,y)` | C `(x,y)` | \(\Delta\mathbf F^{(2)}\) `(x,y)` | \(\Delta\mathbf F^{(3)}\) `(x,y)` |
|---|---:|---|---|---|---|---|
| linear_chain | 0 | (0.6648697224811182,0) | (0.6939581500131806,0) | (0.725204866335105,0) | (0.029088427532062422,0) | (0.031246716321924328,0) |
| linear_chain | 1 | (0,0) | (0,8.006173491782899e-17) | (-9.197387918200184e-17,8.069046260754971e-17) | (0,8.006173491782899e-17) | (-9.197387918200184e-17,6.287276897207196e-19) |
| linear_chain | 2 | (-0.6648697224811182,0) | (-0.6939581500131806,8.498349386392397e-17) | (-0.725204866335105,8.880778617305108e-17) | (-0.029088427532062422,8.498349386392397e-17) | (-0.031246716321924328,3.824292309127114e-18) |
| equilateral | 0 | (0.937368264537491,0.5411898198605332) | (0.9806479527421401,0.5661773594959305) | (1.0281176307855753,0.5935839908926523) | (0.043279688204649114,0.02498753963539735) | (0.04746967804343516,0.027406631396721837) |
| equilateral | 1 | (-0.937368264537491,0.5411898198605332) | (-0.9806479527421403,0.5661773594959304) | (-1.0281176307855755,0.5935839908926523) | (-0.043279688204649336,0.024987539635397238) | (-0.04746967804343516,0.027406631396721948) |
| equilateral | 2 | (0,-1.0823796397210663) | (-2.7755575615628914e-16,-1.1323547189918606) | (-2.299346979550046e-16,-1.1871679817853051) | (-2.7755575615628914e-16,-0.049975079270794254) | (4.7621058201284555e-17,-0.05481326279344456) |
| scalene | 0 | (0.5981114028472363,0.31570942814456443) | (0.6210882513068188,0.3245422033950039) | (0.6326450200727198,0.3306980965625752) | (0.022976848459582477,0.008832775250439462) | (0.01155676876590106,0.0061558931675713136) |
| scalene | 1 | (-0.6172222693777238,0.147057335349943) | (-0.6397215819060202,0.14965332217192848) | (-0.6578317179569251,0.15403905057183756) | (-0.0224993125282964,0.0025959868219854743) | (-0.018110136050904835,0.0043857283999090835) |
| scalene | 2 | (0.019110866530487547,-0.4627667634945074) | (0.018633330599201123,-0.47419552556693245) | (0.01945836875257903,-0.49335070738432696) | (-0.00047753593128642413,-0.011428762072425047) | (0.0008250381533779078,-0.019155181817394507) |

Verification commands: `python -m pytest -q`, `python -m pytest -q -W error`, `git diff --check`, `git status --short`, `git diff --stat`, `git diff --name-only`, and `sha256sum ...`. The audit measured 92 tests without warnings, force-oracle error \(\sim 2.72\times10^{-16}\), \(s_{10}\) error \(\sim 6.97\times10^{-17}\), forbidden mode \(\sim 1.04\times10^{-20}\), sweep residual \(\sim 4.58\times10^{-16}\), and condition number \(\sim 1.403\). The chain central angle is `NaN`; the maximum A--C/B--C angle over 640 chains is \(0^\circ\); 428 contaminated B--C angular entries were corrected. The only sweep changes are the two RMS columns, with \(F_{\mathrm{RMS,new}}=\sqrt{2}F_{\mathrm{RMS,old}}\) and maximum relative deviation \(4.44\times10^{-16}\). The figure was inspected and approved.

T05.1a environment: Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, Matplotlib 3.11.1. Two executions in the same environment are byte-identical; different numerical versions may change final floating-point digits or PNG bytes without scientific randomness. The scientific limits remain Rayleigh, \(N=3\), \(L_{\max}=1\), Models A/B/C only, no Model D, higher multipoles, or \(N>3\). A nonzero scalene sum is not automatically a total cluster force.
