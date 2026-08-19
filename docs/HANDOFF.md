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

## T06: connected \(N=4\) body expansion

### Implementation and equations

Created `src/acoustic_ms/cluster_expansion.py`, `tests/test_four_body.py`, `scripts/validate_t06_quartets.py`, the four T06 artifacts, and the T06 task record. Updated public exports, quartet geometries, README, task list, conventions, decisions, and this handoff. The expansion solves one quartet, the four lexicographic triplets `(0, 1, 2)`, `(0, 1, 3)`, `(0, 2, 3)`, `(1, 2, 3)`, and the isolated pairs supplied by the approved T05 comparison API.

For every particle,

\[
\mathbf F^C=\mathbf F^B+\boldsymbol{\Phi}_{\Sigma}^{(3)}+\boldsymbol{\Phi}^{(4)},
\qquad
\mathbf F^{(\le3)}=\mathbf F^B+\boldsymbol{\Phi}_{\Sigma}^{(3)}.
\]

The recursive and closed forms of \(\boldsymbol{\Phi}^{(4)}\) agree. Model A is used only as a comparison baseline.

### Geometries, regressions, and sweep

Production fixes \(a=E_0=1\), \(ka=0.1\), \(f_0=0\), uses \(f_1\in\{0.1,0.4,0.8,1.0\}\) and 160 equally spaced values of \(d_{\min}/a\in[2.1,10]\). The centered canonical geometries are the four-particle chain, square, and fixed irregular quadrilateral. The sweep has 1.920 configurations and the canonical CSV has 12 particle rows.

| Geometry | Particle | \(F_x^C\) | \(F_y^C\) | \(\Phi_{\Sigma,x}^{(3)}\) | \(\Phi_{\Sigma,y}^{(3)}\) | \(\Phi_x^{(4)}\) | \(\Phi_y^{(4)}\) |
|---|---:|---:|---:|---:|---:|---:|---:|
| linear_chain | 0 | 0.7410016700382743 | 0 | 0.03727662459009518 | 0 | 0.0015397402318367837 | 0 |
| linear_chain | 1 | 0.06998202349087587 | 8.088132637050065e-17 | 0.02998278511626904 | 7.97136892324476e-19 | -0.00019360981048134607 | 2.2454560347183963e-20 |
| linear_chain | 2 | -0.06998202349087597 | 8.948142298346614e-17 | -0.029982785116269336 | 4.460037048956324e-18 | 0.00019360981048154036 | 3.789207058585291e-20 |
| linear_chain | 3 | -0.7410016700382742 | 9.073887836290757e-17 | -0.03727662459009496 | 4.561082570518581e-18 | -0.0015397402318368947 | 1.866535328858388e-19 |
| square | 0 | 0.8218481000644532 | 0.8218481000644532 | 0.05275143070549826 | 0.05275143070549798 | 0.002227631154520071 | 0.002227631154520293 |
| square | 1 | -0.8218481000644532 | 0.8218481000644532 | -0.05275143070549801 | 0.052751430705497995 | -0.002227631154520293 | 0.002227631154520182 |
| square | 2 | -0.8218481000644534 | -0.8218481000644531 | -0.05275143070549815 | -0.05275143070549816 | -0.002227631154520293 | -0.002227631154520182 |
| square | 3 | 0.8218481000644529 | -0.8218481000644532 | 0.05275143070549834 | -0.052751430705498245 | 0.002227631154519738 | -0.002227631154520071 |
| irregular | 0 | 0.8082946943452195 | 0.3226000682788358 | 0.0366319939180071 | 0.016467994763673034 | 0.0012074413637551684 | 0.0004630733328369052 |
| irregular | 1 | -0.6457669671255071 | 0.49471866747294174 | -0.02209008713208649 | 0.023128523747474955 | -0.0006876380032539986 | 0.0004513247992271041 |
| irregular | 2 | 0.4730223188838932 | -0.48229814003930405 | 0.019230383145639768 | -0.02760191017146632 | 0.0007203083619243267 | -0.0006174373380510856 |
| irregular | 3 | -0.6308499575806539 | -0.3449859470334404 | -0.02920333472506556 | -0.021760503592003034 | -0.0011089784059681307 | -0.0004964168626586107 |

Canonical RMS errors \((A,C)/(B,C)/(\le3,C)\) are chain `0.09977255014752409/0.06573916615217347/0.0020850021070293536`, square `0.10419668455946075/0.06689686555910586/0.0027105144543686835`, and irregular `0.0787732847363454/0.047516541713439955/0.0014053466386455625`. Corresponding RMS amplitudes \(\Phi_{\Sigma}^{(3)}/\Phi^{(4)}\) are `0.03382686908479156/0.0010973342107698145`, `0.07460178873829999/0.0031503461906872187`, and `0.035686221945134364/0.001086879416669727`.

### Validation and artifacts

Final verification completed with 103 tests passing with warnings promoted to errors. The independent scalar oracle measured maximum relative errors: force `2.9594212636716024e-16`, \(s_{10}\) `1.569414410681675e-16`, embedded/summed \(\Phi^{(3)}\) `8.915326696739542e-15`, and \(\Phi^{(4)}\) `1.4551295456898224e-13`; the largest forbidden mode was `1.0901074636463385e-20`. The sweep maxima are quartet residual `3.3503948717080633e-16`, quartet condition number `1.4817617280277087`, triplet residual `4.0575294196463723e-16`, and triplet condition number `1.3583134910278372`. Maximum reconstruction and closed-form discrepancies are `1.1102230246251565e-16` and `4.440892098500626e-16`.

Two executions in Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib 3.11.1 produced identical hashes:

- regression CSV: `8d05db59dc4a44ee76118537af40db76aa386c8098f3f87f4359830cb5f9dea0`;
- sweep CSV: `36f64ebd16ea1df52bf4074d42bd83356306bfc17c613267b0def746901689b5`;
- model-error figure: `547945fa0fd658565fb837416f8f4a5c65bd4963c0e2e50226510f82f7af17d0`;
- body-decomposition figure: `e4ce5f5a0c5f1d212d72ddd1bc6d35e1ecded8e86e437efeaf4a19bce4ab6d16`.

Protected hashes remain T03 `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`, T04 `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`, T05 regression `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`, T05 sweep `dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1`, and T05 figure `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`.

The figures were visually inspected: panels are populated, axes and legends are legible, and curves are distinguishable. Scientific limits remain planar \(N=4\), Rayleigh \(L_{\max}=1\), Models A/B/C, external--scattered interaction force only, no higher multipoles, no unrestricted total-force interpretation for the irregular quartet, and no Model D. The recommended next step is Model D in a separate task.

## T06.1: connected-body scaling analysis

### Scope, implementation, and data provenance

T06.1 created `src/acoustic_ms/scaling.py`, `tests/test_scaling_analysis.py`,
`scripts/analyze_t06_scaling.py`, `TAREFA_T06_1_ANALISE_ESCALA.md`, three CSVs,
and two figures. It updated `src/acoustic_ms/__init__.py`, `README.md`,
`TASKS.md`, `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and this handoff.
No protected solver, force model, earlier test, or T03--T06 artifact changed.

The analysis reads all 1,920 existing T06 sweep rows directly: three
geometries, four contrasts \(f_1\in\{0.1,0.4,0.8,1.0\}\), and 160 distances
per geometry--contrast pair. It performs no new trimer or quartet force sweep.
The sole new Model C evaluation is the required centered dimer at
\(d/a=2.1\), evaluated with the public pair APIs for the seven-row body-order
table.

The predictors and responses are

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3,
\qquad
\Lambda_{\max}=\max_i\left[|f_1|\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3\right],
\]

\[
Y_3=\frac{F_{\mathrm{RMS}}(\boldsymbol{\Phi}_{\Sigma}^{(3)})}
{F_{\mathrm{RMS}}(\mathbf F^C)},
\qquad
Y_4=\frac{F_{\mathrm{RMS}}(\boldsymbol{\Phi}^{(4)})}
{F_{\mathrm{RMS}}(\mathbf F^C)}.
\]

Every unweighted fit uses every positive point in
\(\ln y=\ln C+p\ln x\). The diagnostics are \(R^2_{\log}\), log-space RMSE,
and maximum absolute log residual. No data are discarded, binned, smoothed,
weighted, or treated statistically.

### Geometric factors and 16 fits

The measured factors \(C_g=\Lambda_{\max}/\eta\) are:

| Geometry | \(C_g\) |
|---|---:|
| linear chain | 2.125 |
| square | 2.353553390593274 |
| irregular | 1.996580257145743 |

| Group | Predictor | Response | Points | \(C\) | \(p\) | \(p/o\) | \(R^2_{\log}\) | RMSE\(_{\log}\) | max \(|r_{\log}|\) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| linear chain | \(\eta\) | \(Y_3\) | 640 | 0.6536500039179222 | 0.9477457629207996 | 0.9477457629207996 | 0.9984681402009536 | 0.05819892689323387 | 0.11783376232604859 |
| linear chain | \(\eta\) | \(Y_4\) | 640 | 0.22802259922396498 | 1.9209166221002902 | 0.9604583110501451 | 0.9991340421510149 | 0.08865959328065456 | 0.1821914636121864 |
| linear chain | \(\Lambda_{\max}\) | \(Y_3\) | 640 | 0.3199574486842558 | 0.947745762920799 | 0.947745762920799 | 0.9984681402009536 | 0.05819892689323382 | 0.11783376232604326 |
| linear chain | \(\Lambda_{\max}\) | \(Y_4\) | 640 | 0.053598010989882774 | 1.9209166221002894 | 0.9604583110501447 | 0.9991340421510149 | 0.08865959328065454 | 0.1821914636121953 |
| square | \(\eta\) | \(Y_3\) | 640 | 0.6406525752847945 | 0.9487770531356282 | 0.9487770531356282 | 0.9983499381833476 | 0.060471898555038685 | 0.12137246541572999 |
| square | \(\eta\) | \(Y_4\) | 640 | 0.27236508623737304 | 1.9198426051755004 | 0.9599213025877502 | 0.9988356442244823 | 0.10276427955401718 | 0.210058674038331 |
| square | \(\Lambda_{\max}\) | \(Y_3\) | 640 | 0.28440637242529476 | 0.9487770531356281 | 0.9487770531356281 | 0.9983499381833476 | 0.06047189855503868 | 0.12137246541572733 |
| square | \(\Lambda_{\max}\) | \(Y_4\) | 640 | 0.05266229214840381 | 1.919842605175501 | 0.9599213025877505 | 0.9988356442244823 | 0.10276427955401728 | 0.210058674038331 |
| irregular | \(\eta\) | \(Y_3\) | 640 | 0.4552015284021664 | 0.9427068521180625 | 0.9427068521180625 | 0.9979957125494455 | 0.0662327902550376 | 0.13170610505534341 |
| irregular | \(\eta\) | \(Y_4\) | 640 | 0.14364901647817394 | 1.9128217876764342 | 0.9564108938382171 | 0.9987739440349565 | 0.10506952695317552 | 0.22289695800020937 |
| irregular | \(\Lambda_{\max}\) | \(Y_3\) | 640 | 0.23720362061363207 | 0.9427068521180623 | 0.9427068521180623 | 0.9979957125494455 | 0.0662327902550376 | 0.1317061050553443 |
| irregular | \(\Lambda_{\max}\) | \(Y_4\) | 640 | 0.0382743281607928 | 1.9128217876764346 | 0.9564108938382173 | 0.9987739440349566 | 0.10506952695317551 | 0.2228969580002076 |
| grouped | \(\eta\) | \(Y_3\) | 1920 | 0.575515712845375 | 0.9464098893914968 | 0.9464098893914968 | 0.9881180078496571 | 0.16270474366899415 | 0.343368958567277 |
| grouped | \(\eta\) | \(Y_4\) | 1920 | 0.20740091186607565 | 1.917860338317407 | 0.9589301691587035 | 0.9920672952487812 | 0.2688671716039741 | 0.5606023914972305 |
| grouped | \(\Lambda_{\max}\) | \(Y_3\) | 1920 | 0.2803940168979942 | 0.9477200372149447 | 0.9477200372149447 | 0.9927091201707364 | 0.1274516770701454 | 0.2716098723293241 |
| grouped | \(\Lambda_{\max}\) | \(Y_4\) | 1920 | 0.04835618532017296 | 1.9207265477838105 | 0.9603632738919052 | 0.996896030515444 | 0.16818445600557255 | 0.4160067894461079 |

Here \(o=1\) for \(Y_3\) and \(o=2\) for \(Y_4\), so \(p/o\) is the
exponent relative to the transformed predictor \(x^o\). Within each fixed
geometry, \(\Lambda_{\max}=C_g\eta\), and the predictor replacement preserves
\(p\), \(R^2_{\log}\), RMSE, and maximum residual to rounding while changing
the prefactor.

For the grouped fits, changing from \(\eta\) to \(\Lambda_{\max}\) reduces
log-space RMSE by 0.21666895386017415 for \(Y_3\) and
0.37447009613617466 for \(Y_4\), or approximately 21.67% and 37.45%.
This is a descriptive improvement for these three families, not evidence of a
universal coupling criterion.

### Canonical \(N=2,3,4\) comparison

All rows use \(a=E_0=1\), \(ka=0.1\), \(f_0=0\), \(f_1=0.8\), and
\(d_{\min}/a=2.1\). Zeros in unavailable body orders are structural zeros.

| \(N\) | Geometry | A--C | B--C | \((\le3)\)--C | \(\Phi_3/F^C\) | \(\Phi_4/F^C\) | \(F_{\mathrm{RMS}}(F^C)/(a^2E_0)\) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 2 | pair | 0.044133766948299784 | 0 | 0 | 0 | 0 | 0.6537653018280931 |
| 3 | linear chain | 0.08319737863715175 | 0.04308674386016288 | 0 | 0.04308674386016288 | 0 | 0.5921272938347617 |
| 3 | equilateral | 0.08826749345670055 | 0.04617144636179811 | 0 | 0.04617144636179811 | 0 | 1.1871679817853047 |
| 3 | scalene | 0.06634541632888573 | 0.031143703608647644 | 0 | 0.031143703608647644 | 0 | 0.7711809952716625 |
| 4 | linear chain | 0.09977255014752409 | 0.06573916615217347 | 0.0020850021070293536 | 0.06427311991532472 | 0.0020850021070293536 | 0.5262988498046471 |
| 4 | square | 0.10419668455946075 | 0.06689686555910586 | 0.0027105144543686835 | 0.06418635110473715 | 0.0027105144543686835 | 1.1622687293217102 |
| 4 | irregular | 0.0787732847363454 | 0.047516541713439955 | 0.0014053466386455625 | 0.046142664298696304 | 0.0014053466386455622 | 0.7733888471226102 |

The T05 regression CSV stores the scalene canonical row at \(d/a=2.2\), not
at the table's required \(2.1\). To preserve the no-new-trimer-solve rule, the
scalene \(F^C\) RMS at \(2.1\) is recovered from the exact sweep identity
`rms_irreducible_multibody / rms_b_vs_c`; the CSV records this provenance
explicitly. All other trimer/quartet amplitudes come from their per-particle
regression CSVs.

### Verification, determinism, and hashes

The observed environment is Python 3.12.3, NumPy 2.5.1, and Matplotlib 3.11.1.
The final suite reports 124 tests passing with warnings treated as errors. The
script validates all four input CSVs, and two consecutive executions produce
byte-identical outputs with 16 fit rows, two collapse rows, and seven
body-order rows.

New artifact SHA-256 hashes:

- `t06_1_scaling_fits.csv`: `2e71a1b5cb5df238b3bd9ed94f96ef36cae36964af87a77cef89238c0b4d3367`;
- `t06_1_collapse_summary.csv`: `e8120962b77a8fdd09b1f9209f939191fe1aba9b597731aedba201f92b4eca8e`;
- `t06_1_body_order_summary.csv`: `4258146c21e0a183fe3baa3e3f92d4e8bd59e0782344c8687fe5600a6b785ae9`;
- `t06_1_eta_scaling.png`: `4aef5f9654533c5b2e694e0feaf0933a84266aef1c06c4d1862ac8289d32a064`;
- `t06_1_lambda_scaling.png`: `f817b7ce31a3ee6ed2b829dd98fa8d0534a18e5e6633f51f796b547ebd4720c8`.

Preserved SHA-256 hashes:

- T03 CSV: `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`;
- T04 CSV: `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`;
- T05 regression: `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`;
- T05 sweep: `dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1`;
- T05 figure: `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`;
- T06 regression: `8d05db59dc4a44ee76118537af40db76aa386c8098f3f87f4359830cb5f9dea0`;
- T06 sweep: `36f64ebd16ea1df52bf4074d42bd83356306bfc17c613267b0def746901689b5`;
- T06 model figure: `547945fa0fd658565fb837416f8f4a5c65bd4963c0e2e50226510f82f7af17d0`;
- T06 body-decomposition figure: `e4ce5f5a0c5f1d212d72ddd1bc6d35e1ecded8e86e437efeaf4a19bce4ab6d16`.

Verification commands include editable installation, `python -m pytest -q`,
`python -m pytest -q -W error`, two executions of
`python scripts/analyze_t06_scaling.py`, SHA-256 checks, CSV row/finiteness
checks, control-character inspection, `git diff --check`, and Git scope
inspection. Both 2400×1000 PNGs were visually inspected: axes and legends are
legible, all geometry/contrast series and grouped fits are present, the axes
are logarithmic, and no clipping, empty panels, `NaN`, or `inf` is visible.

The numerical evidence is limited to \(ka=0.1\), positive \(f_1\), three fixed
planar dilation families, Rayleigh \(L_{\max}=1\), and the external--scattered
force observable. The sweep varies \(kd\), whereas the exact coupling contains
Hankel functions rather than only \(r^{-3}\). The fitted exponents are not an
analytic proof. The amplitude of the embedded vector sum
\(\boldsymbol{\Phi}_{\Sigma}^{(3)}\) is not the sum of individual trimer
amplitudes, and multipolar corrections may exceed \(\boldsymbol{\Phi}^{(4)}\).
T07 and Model D were not started. The recommended next step is a separate T07
beginning with multipolar convergence for \(N=2\).

## T07 — Model D and multipolar convergence

### Implementation and protected scope

T07 added `multipolar_scattering.py`, `multipolar_solver.py`, `model_d.py`, and `multipolar_expansion.py`; exported their public APIs; added three T07 test modules and `validate_t07_multipolar.py`; generated three CSVs and two PNGs; and updated this documentation, conventions, decisions, task registry, README, and the self-contained T07 specification. No protected T01--T06.1 scientific module, test, script, or artifact was modified.

The positive-order Rayleigh coefficients are

\[
s_\ell=i\frac{3\ell f_1}{(2\ell-1)!!(2\ell+1)!![2(2\ell+1)-(\ell-1)f_1]}(ka)^{2\ell+1}.
\]

The solver uses the complete project ordering, a planar active basis with \(\ell+m\) odd, principal complex square roots in the balanced equation, and reports a residual evaluated in the original physical equation. Local fields are reexpanded with source order \(L_{\max}\) and target order 2; the T04 force expression is unchanged.

### Independent validation

Across an isolated particle, two dimers, the three canonical trimers, and the three canonical quartets, the maximum absolute differences between Model D at \(L=1\) and Model C were \(4.44\times10^{-16}\) for force components and \(4.41\times10^{-19}\) for scattering coefficients. The corrected one-particle symmetry reduction of the dimer agreed with the global solution at \(L=1,3,5\); the largest observed coefficient difference was \(1.13\times10^{-19}\).

The strict odd-order branch used in the derivation of Eq. (30) gave relative discrepancies \(1.7061326608\times10^{-3}\), \(4.3074064108\times10^{-4}\), and \(1.0794930439\times10^{-4}\) at \(ka=0.1,0.05,0.025\), respectively. Thus its asymptotic error decreases by approximately a factor four per halving of \(ka\). The general planar Model D retains additional symmetry-allowed even-\(\ell\), odd-\(m\) channels; its direct discrepancy from Eq. (30) is separately recorded in the analytic CSV and must not be conflated with the strict reduced benchmark.

### Dimer convergence

The base map contains exactly \(5\times4\times5=100\) rows. Nine targeted \(L=11\) rows were added only where needed to test or confirm convergence. Minimum confirmed orders were:

| \(d/a\) | \(f_1=0.1\) | \(f_1=0.4\) | \(f_1=0.8\) | \(f_1=1.0\) |
|---:|---:|---:|---:|---:|
| 2.00 | 7 | 9 | not confirmed through 11 | not confirmed through 11 |
| 2.05 | 7 | 9 | not confirmed through 11 | not confirmed through 11 |
| 2.10 | 7 | 9 | not confirmed through 11 | not confirmed through 11 |
| 2.50 | 5 | 7 | 7 | 7 |
| 3.00 | 5 | 5 | 5 | 7 |

For the stress case \(d/a=2,f_1=1\), the \(L=11\) successive force error is \(6.0377621229\times10^{-3}\); no convergence claim is made.

### Canonical clusters

All subsets were solved at the same order. The final available RMS force and successive errors are:

| geometry | final \(L\) | \(F_{\rm RMS}(\mathbf F^{D,L})\) | \(\epsilon_L^F\) | \(\epsilon_L^{(3)}\) | \(\epsilon_L^{(4)}\) |
|---|---:|---:|---:|---:|---:|
| trimer linear | 11 | 0.6518606954444777 | 1.45622e-4 | 2.35237e-3 | not applicable |
| trimer equilateral | 11 | 1.417220834396371 | 2.11264e-4 | 2.81969e-4 | not applicable |
| trimer scalene | 9 | 0.8611163977289301 | 8.09874e-4 | 6.69608e-4 | not applicable |
| quartet linear | 13 | 0.5750165079059795 | 3.01722e-5 | 3.40585e-4 | 2.07603e-4 |
| quartet square | 11 | 1.344690911288789 | 1.82564e-4 | 2.28772e-4 | 3.39363e-4 |
| quartet irregular | 9 | 0.8492744782490987 | 7.15549e-4 | 5.43663e-4 | 7.84447e-4 |

The linear trimer's total force is converged by the force criterion, but its three-body term is not confirmed through \(L=11\); this distinction is intentional. Four clusters required \(L=11\), and only the linear quartet required \(L=13\). The maximum physical residual was \(8.01\times10^{-16}\). The largest balanced condition number was 1.6753, while the raw matrix reached \(5.57\times10^{40}\); the latter is scaling pathology, not physical divergence.

### Artifacts and environment

The verification environment was Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib 3.11.1. Two consecutive runs in this environment were byte-identical. Official T07 hashes:

```text
c0e115f87a8d80a58f7a9188c69be4ce0f0d60518326bd3128726b369c93735e  results/data/t07_pair_analytic_validation.csv
fba9e51c5a93c0161a3a04c9ff98505ee2032d3e42f2f55490ca644cf4eb2dc6  results/data/t07_dimer_convergence.csv
a272d80770a19a015a5bcc6245a2d35eb2acb165a3b6fd79bae0947c988b1db0  results/data/t07_cluster_convergence.csv
c38297886a664bb82193a1a500812902b85daf315ff9ede196705837715d65c7  results/figures/t07_dimer_convergence.png
0b19fa7b7cffdbf06e1144ad72f881393b3fc658b43586fd969f9ee59e6c45de  results/figures/t07_cluster_convergence.png
```

The pre-existing T03--T06.1 artifact hashes were rechecked byte for byte. Both figures were inspected for readable labels, visible tolerance lines, distinct geometry/contrast curves, and explicit adaptive orders.

### Commands and limitations

Verification completed with 160/160 tests passing and no warnings. It included editable installation, `pytest -q`, `pytest -q -W error`, two executions of the T07 validator, SHA-256 comparisons, `git diff --check`, protected-file diff checks, and an ASCII-control scan. No T05/T06 1,920-case sweep was rerun.

Scientific scope remains \(ka\leq0.1\), \(N\leq4\), identical fixed spheres in the nodal plane, leading Rayleigh coefficients at each multipole, and the external--scattered force. Exact Mie coefficients, scattered--scattered forces, viscosity, streaming, walls, torque, and dynamics remain outside scope. No universal multipolar cutoff is inferred from these canonical cases.

## T08 — transferability through \(N=10\) and frozen article data

### Scope, implementation, and provenance

T08 started from `bd3752ef2d4fa61ec10887445575f507d8c8cd6b` and leaves all
T01--T07 scientific modules and artifacts unchanged. It adds the deterministic
cluster families, matched multipolar pairwise baseline, transferability
analysis, three focused test files, raw and derived CSVs, two figures, and the
self-contained T08 specification. The exact changed-file set is:

```text
TAREFA_T08_TRANSFERIBILIDADE_CRITERIO_ACOPLAMENTO.md
README.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
src/acoustic_ms/__init__.py
src/acoustic_ms/cluster_families.py
src/acoustic_ms/transferability.py
scripts/run_t08_transferability.py
scripts/analyze_t08_transferability.py
tests/test_t08_cluster_families.py
tests/test_t08_transferability.py
tests/test_t08_analysis.py
results/data/t08_cases.csv
results/data/t08_forces.csv
results/data/t08_convergence.csv
results/data/t08_predictor_fits.csv
results/data/t08_validity_thresholds.csv
results/figures/t08_predictor_comparison.png
results/figures/t08_transferability.png
```

The physical grid is fixed at \(a=E_0=1\), \(k=ka=0.1\), \(f_0=0\),
\(f_1\in\{0.1,0.4,0.8,1.0\}\), and
\(d_{\min}/a\in\{2.1,2.5,3,4,6,10\}\). The thirteen deterministic
geometries are one pair and linear, compact, and irregular families for each
of \(N=3,4,6,10\). This produces exactly 312 configurations: 168 calibration
cases with \(N\leq4\) and 144 external-holdout cases with \(N=6,10\).
The long force table contains 1,704 particle rows and the adaptive convergence
table contains 1,270 evaluated orders.

The new matched baseline is

\[
\mathbf F_i^{B_L}=\sum_{j\ne i}\mathbf F_{ij}^{D,N=2,L}.
\]

Each unordered dimer is solved at the same \(L\) as the global cluster and is
cached deterministically by distance and physical parameters. Across the 13
families at a representative point, the maximum component difference in
\(B_1=B\) was \(4.440892098500626\times10^{-16}\). For all dimer rows,
\(B_L=D_L\) to a maximum recorded relative discrepancy of
\(2.703147394448881\times10^{-16}\). The vector identity

\[
\mathbf F^{D_L}-\mathbf F^A=
(\mathbf F^{B_L}-\mathbf F^A)+(\mathbf F^{D_L}-\mathbf F^{B_L})
\]

closed with maximum absolute component error exactly zero in the frozen CSV.

The analyzed errors and amplitude diagnostics are

\[
\varepsilon_A=\frac{F_{\mathrm{RMS}}(\mathbf F^A-\mathbf F^D)}
{F_{\mathrm{RMS}}(\mathbf F^D)},
\qquad
\varepsilon_B=\frac{F_{\mathrm{RMS}}(\mathbf F^{B_L}-\mathbf F^D)}
{F_{\mathrm{RMS}}(\mathbf F^D)},
\]

\[
Y_{\mathrm{2B}}=\frac{F_{\mathrm{RMS}}(\mathbf F^{B_L}-\mathbf F^A)}
{F_{\mathrm{RMS}}(\mathbf F^D)},\quad
Y_{\mathrm{coll}}=\frac{F_{\mathrm{RMS}}(\mathbf F^D-\mathbf F^{B_L})}
{F_{\mathrm{RMS}}(\mathbf F^D)},\quad
Y_{\mathrm{mp}}=\frac{F_{\mathrm{RMS}}(\mathbf F^D-\mathbf F^{D_1})}
{F_{\mathrm{RMS}}(\mathbf F^D)}.
\]

### Convergence and numerical audit

Every case was evaluated at \(L=1,3,5,7,9,11\), stopping early only after
both latest normalized differences satisfied \(10^{-3}\); \(L=13\) was
available only for \(N\leq4\). Total, matched-pairwise, and collective
residual convergence were assessed separately. Eleven calibration cases
required \(L=13\): 2 dimers, 5 trimers, and 4 quartets.
The final status by particle count and family is:

| \(N\) | family | cases | total confirmed | joint confirmed | residual resolved |
|---:|---|---:|---:|---:|---:|
| 2 | pair | 24 | 24 | 24 | 3 |
| 3 | linear | 24 | 24 | 24 | 23 |
| 3 | compact | 24 | 24 | 24 | 24 |
| 3 | irregular | 24 | 24 | 24 | 24 |
| 4 | linear | 24 | 24 | 24 | 23 |
| 4 | compact | 24 | 24 | 24 | 24 |
| 4 | irregular | 24 | 24 | 24 | 24 |
| 6 | linear | 24 | 23 | 23 | 22 |
| 6 | compact | 24 | 22 | 22 | 22 |
| 6 | irregular | 24 | 23 | 23 | 23 |
| 10 | linear | 24 | 23 | 23 | 22 |
| 10 | compact | 24 | 23 | 23 | 23 |
| 10 | irregular | 24 | 23 | 23 | 23 |

Seven close, strong holdout cases reached the \(L=11\) limit without two-step
confirmation and are explicitly `unconfirmed`, not divergent:

| case | penultimate \(\delta_D\) | last \(\delta_D\) |
|---|---:|---:|
| `n6_linear_f1.0_d2.1` | 0.001456726809597844 | 0.00032289030982289487 |
| `n6_compact_f0.8_d2.1` | 0.0010763361885936945 | 0.0001852716698896735 |
| `n6_compact_f1.0_d2.1` | 0.0020863427787391302 | 0.0004162852041807814 |
| `n6_irregular_f1.0_d2.1` | 0.0010913381614902689 | 0.00023661282579013175 |
| `n10_linear_f1.0_d2.1` | 0.0014522538245897031 | 0.00032191409609328694 |
| `n10_compact_f1.0_d2.1` | 0.0019299115954099244 | 0.00038728538373673338 |
| `n10_irregular_f1.0_d2.1` | 0.001095178367888512 | 0.00023650133997409569 |

Holdout total-convergence coverage is \(137/144=0.9513888888888888\).
The maximum physical residual is \(1.2851934906867192\times10^{-15}\), the
maximum balanced condition number is 2.236483305457277, and the maximum raw
condition number is \(2.2449893674832796\times10^{43}\). The raw number is a
basis-scaling diagnostic and was not used to judge convergence. All forces,
required metrics, residues, and conditions are finite.

The six canonical T07 clusters were compared at every overlapping stored
order: all 35 RMS-force comparisons had zero absolute and relative difference
in this environment. No T07 CSV was modified.

### Predictor calibration and leakage-safe selection

The three predictors are

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3,
\qquad
\Lambda_{\max}=|f_1|\max_i\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3,
\]

\[
\rho_1=\max_\nu|\lambda_\nu(\mathbf I-\mathbf A_b^{(1)})|.
\]

All fits use only eligible calibration cases. The \(\varepsilon_B\) fit also
excludes dimers, unresolved residuals, and cases without joint convergence.

| predictor | response | points | \(C\) | \(p\) | \(R^2_{\log}\) | RMSE\(_{\log}\) | max \(|r_{\log}|\) | Spearman |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| \(\eta\) | \(\varepsilon_A\) | 168 | 1.8559963067143632 | 1.0565911785455437 | 0.9811693912231740 | 0.26950576804208048 | 0.7018318973222302 | 0.9877944657410198 |
| \(\eta\) | \(\varepsilon_B\) | 142 | 0.47433468176396937 | 0.9499523012966420 | 0.9550720993519037 | 0.38026694230919372 | 0.9662089673176322 | 0.9727713531608475 |
| \(\Lambda_{\max}\) | \(\varepsilon_A\) | 168 | 0.9484449956388588 | 1.0487353446836789 | 0.9861594514823149 | 0.23105349130712635 | 0.8069362842422012 | 0.9916858600913990 |
| \(\Lambda_{\max}\) | \(\varepsilon_B\) | 142 | 0.24645304782633867 | 0.9530358680835566 | 0.9641571551868680 | 0.33965009913058924 | 0.9580978250544607 | 0.9784272925999035 |
| \(\rho_1\) | \(\varepsilon_A\) | 168 | 2.6353684041458636 | 1.1088518115798773 | 0.9872268644670024 | 0.22196507422252662 | 0.8018841610319449 | 0.9922707416949703 |
| \(\rho_1\) | \(\varepsilon_B\) | 142 | 0.6506524433812200 | 1.0144638088529034 | 0.9787397449907628 | 0.26158614905731969 | 0.7845782241353017 | 0.9871581819439177 |

Leave-\((N,\mathrm{family})\)-out validation on \(\varepsilon_A\) gave
log-RMSE 0.29572460861555683 for \(\eta\), 0.24630178255151594 for
\(\Lambda_{\max}\), and 0.23203390779877014 for \(\rho_1\). Therefore
\(\rho_1\) was selected without consulting the holdout.

### External holdout and empirical thresholds

For frozen \(\rho_1\) calibration predictions, the holdout results are:

| scope | eligible | RMSE\(_{\log}\) | median factor | p90 factor | max factor | within factor 2 | Spearman |
|---|---:|---:|---:|---:|---:|---:|---:|
| all holdout | 137 | 0.2320007775714556 | 1.2294242478993347 | 1.4345931747220635 | 1.5933123541981697 | 1.0 | 0.9970504228349543 |
| \(N=6\) | 68 | 0.1965479872155958 | 1.1853099998758223 | 1.3508031534732079 | 1.4452142044835998 | 1.0 | 0.9967553536664503 |
| \(N=10\) | 69 | 0.26229261958175215 | 1.264449698281531 | 1.4834377134937253 | 1.5933123541981697 | 1.0 | 0.9965655827548413 |

The conservative thresholds were fixed from calibration prefixes only:

| tolerance | \(\rho_{1,\tau}\) | calibration count | safe holdout | coverage | false safe | worst safe error |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.0053990295322641655 | 75 | 48 | 0.35036496350364965 | 0 | 0.007151263247680959 |
| 0.05 | 0.02000077753569526 | 117 | 86 | 0.6277372262773723 | 0 | 0.03756638733704074 |
| 0.10 | 0.03914887870730305 | 141 | 101 | 0.7372262773722628 | 0 | 0.05274196806688931 |

All five prespecified diagnostic conditions are satisfied, so
`criterion_supported = true`. This means only that \(\rho_1\) improved the
transferable description for these deterministic nodal clusters and sampled
parameters. It is not a universal validity criterion, an analytic proof, or a
licensed extrapolation beyond the domain.

### Determinism, artifacts, and environment

The complete expensive sweep was executed once. Its three raw files were then
audited by recomputing eight stratified cases covering \(N=2,4,10\), weak and
strong contrast, near and far separation, and all three \(N=10\) families.
The audit agreed numerically and did not rewrite the raw files. The analysis
was executed twice; both derived CSVs and both PNGs were byte-identical.
No 1,920-case T05 or T06 sweep was rerun.

The verification environment was Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0,
Matplotlib 3.11.1, on Linux 7.0.0-28-generic x86_64 with glibc 2.39. Binary
reproducibility is asserted only within this numerical environment.

Official T08 SHA-256 hashes:

```text
62e013b8429846cb085633cfe2fb0eb530e944dd8059484bd20055a02a6b94fb  results/data/t08_cases.csv
73f870eaa2774a47d2eb83d07618b29290d83ee16516f027a5245e6e6d4beed4  results/data/t08_forces.csv
87e0320b7ff023884213c562997912fb7492b5aec7bd5caa505a3039e037f3a3  results/data/t08_convergence.csv
707024d14c79b4ee6aebc5e20ae0037b33fd65790a082b9b561e60dee30776a1  results/data/t08_predictor_fits.csv
27ad27a0b5b2a7bc4521118d4b4e32c04e77a6ec19085141e56c9ac0e3114f2c  results/data/t08_validity_thresholds.csv
f67c0c41f0c590af6127a50f9041da5f22b082049ac568cc1d74d3f260269578  results/figures/t08_predictor_comparison.png
078d6543fa819033e01389d2c76bb81c4d8a0e2bdf34b4838dc25748bb48723e  results/figures/t08_transferability.png
```

The full pre-T08 results manifest was compared after generation. Every one of
its 21 hashes remained unchanged; the audit includes Figure 2 and all
T03--T07 data and figures. The hashes are:

```text
c19b025db20f97b04a19b8ec14afbdf4440760b980ee533a464528b29074de40  results/data/figure_2_relative_error.csv
7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a  results/data/t03_solver_validation.csv
15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5  results/data/t04_pair_force_validation.csv
e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1  results/data/t05_trimer_regression.csv
dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1  results/data/t05_trimer_sweep.csv
4258146c21e0a183fe3baa3e3f92d4e8bd59e0782344c8687fe5600a6b785ae9  results/data/t06_1_body_order_summary.csv
e8120962b77a8fdd09b1f9209f939191fe1aba9b597731aedba201f92b4eca8e  results/data/t06_1_collapse_summary.csv
2e71a1b5cb5df238b3bd9ed94f96ef36cae36964af87a77cef89238c0b4d3367  results/data/t06_1_scaling_fits.csv
8d05db59dc4a44ee76118537af40db76aa386c8098f3f87f4359830cb5f9dea0  results/data/t06_quartet_regression.csv
36f64ebd16ea1df52bf4074d42bd83356306bfc17c613267b0def746901689b5  results/data/t06_quartet_sweep.csv
a272d80770a19a015a5bcc6245a2d35eb2acb165a3b6fd79bae0947c988b1db0  results/data/t07_cluster_convergence.csv
fba9e51c5a93c0161a3a04c9ff98505ee2032d3e42f2f55490ca644cf4eb2dc6  results/data/t07_dimer_convergence.csv
c0e115f87a8d80a58f7a9188c69be4ce0f0d60518326bd3128726b369c93735e  results/data/t07_pair_analytic_validation.csv
678cb8f086d7dacfd9be9f7960556dc8ad6ccdac3a38335fe5a3827719490f05  results/figures/figure_2_relative_error.png
5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62  results/figures/t05_trimer_model_errors.png
4aef5f9654533c5b2e694e0feaf0933a84266aef1c06c4d1862ac8289d32a064  results/figures/t06_1_eta_scaling.png
f817b7ce31a3ee6ed2b829dd98fa8d0534a18e5e6633f51f796b547ebd4720c8  results/figures/t06_1_lambda_scaling.png
e4ce5f5a0c5f1d212d72ddd1bc6d35e1ecded8e86e437efeaf4a19bce4ab6d16  results/figures/t06_quartet_body_decomposition.png
547945fa0fd658565fb837416f8f4a5c65bd4963c0e2e50226510f82f7af17d0  results/figures/t06_quartet_model_errors.png
0b19fa7b7cffdbf06e1144ad72f881393b3fc658b43586fd969f9ee59e6c45de  results/figures/t07_cluster_convergence.png
c38297886a664bb82193a1a500812902b85daf315ff9ede196705837715d65c7  results/figures/t07_dimer_convergence.png
```

Both T08 figures were visually inspected. Their logarithmic axes, labels,
legends, calibration/holdout distinction, family and particle-count encodings,
calibration-only fit lines, open holdout and crossed unconfirmed markers, tolerance lines, and the
observed-versus-predicted factor-two band are legible and populated without
clipping or physical `NaN`/`inf` values.

### Commands and final limitations

Final verification reports 193/193 tests passing, including warnings promoted
to errors.

The workflow used editable installation; baseline and final `pytest -q` and
`pytest -q -W error`; one full `run_t08_transferability.py`; one
`run_t08_transferability.py --audit-existing`; two executions of
`analyze_t08_transferability.py`; SHA-256 manifests before and after; CSV
count, uniqueness, finiteness, and identity audits; ASCII-control scans;
visual PNG inspection; `git diff --check`; and protected-path inspection.

The evidence is limited to the nodal plane, \(ka=0.1\), positive \(f_1\),
identical fixed spheres, \(N\leq10\), and only three deterministic families at
\(N=6,10\). It uses the external--scattered force and leading Rayleigh
coefficient of every retained multipole. It excludes scattered--scattered
forces, exact Mie T-matrices, negative contrasts, random ensembles, viscosity,
streaming, walls, contact, torque, trajectories, and dynamics. The empirical
thresholds apply only to the sampled domain. T08 is the final computational
sweep planned for the article, and the committed datasets are frozen.

## T09 — analytical foundation of \(\rho_1\)

T09 leaves the frozen T08 tables untouched and analytically reconstructs the
balanced \(L=1\) rescattering operator. With \(f_0=0\), planar nodal symmetry
leaves one \((1,0)\) channel per particle. Combining

\[
s_1=i f_1(ka)^3/6
\]

with

\[
h_0^{(1)}(x)+h_2^{(1)}(x)
=-3e^{ix}(x+i)/x^3
\]

gives

\[
(K_b)_{ii}=0,\qquad
(K_b)_{ij}
=
\frac{f_1}{2}\left(\frac{a}{r_{ij}}\right)^3
e^{ikr_{ij}}(1-ikr_{ij}).
\]

Its near-field limit is \(f_1(a/r_{ij})^3/2\), which derives the scales
\(\eta\) and \(\Lambda_{\max}\). For a dimer,

\[
\rho_1
=
\frac{|f_1|}{2}\left(\frac{a}{d}\right)^3
\sqrt{1+(kd)^2}.
\]

The expansion

\[
(\mathbf I-\mathbf K_b)^{-1}
=
\sum_{p=0}^{\infty}\mathbf K_b^p
\]

converges as a matrix series exactly when \(\rho_1<1\). The \(p=0\) force is Model A in the
\(L=1\) observable and the omitted correction begins at \(O(\mathbf K_b)\),
explaining the fitted exponent near one. A connected \(n\)-body path first
requires \(p=n-2\), so \(\Phi^{(3)}=O(\rho_1)\) and
\(\Phi^{(4)}=O(\rho_1^2)\), consistent with T06.1.

The independent closed matrix reproduced all 312 frozen T08 radii with
maximum absolute difference \(4.163336342344337\times10^{-16}\). Across that
domain, the maximum radius was \(0.2544601331856266\), the maximum
\(\|\mathbf K_b\|_2/\rho_1\) was \(1.0174745602981063\), and the maximum
\(\|\mathbf K_b\|_\infty/\rho_1\) was \(1.4935853868714868\). Thus every
sampled operator is Neumann-convergent and close to normal in the
spectral-norm ratio, although non-normality remains an explicit theoretical
caveat.

SymPy verified the Hankel identity and

\[
e^{ix}(1-ix)
=1+x^2/2+ix^3/3-x^4/8+O(x^5).
\]

No Mathematica step was required. Three sentinel cases with
\(\rho_1=0.00599,\ 0.0672,\ 0.254\) showed monotone convergence of the
partial Neumann sums to the direct solution down to the floating-point floor.
T09 adds 24 focused tests. Final verification gives **217/217 tests passing**
with warnings promoted to errors.

Official T09 SHA-256 hashes:

```text
00278f05d92b4040bf2abd1572e1073bc69c9aec0048e096fba4d80dbdb30ff9  results/data/t09_analytic_summary.csv
be512cfe95c19fa6491e8f43a6f5f1645c56648c48881f43805a26d26c751f29  results/data/t09_neumann_convergence.csv
9ff6a120ceca132203ba83de88000d7e94ee79845e70d0a7f3ef81fb8db24e92  results/data/t09_operator_audit.csv
2abb82bdb19627b813ec45693483bb4d0f998598707ea1a750b62c8220f4f6ea  results/figures/t09_rho_foundation.png
```

The analytical result justifies the order and physical meaning of the
predictor. It does not derive the empirical prefactor \(2.635\), make the
exponent \(1.109\) exact, or turn the T08 thresholds into universal
constants. The complete derivation and artifact list are in
`TAREFA_T09_FUNDAMENTACAO_ANALITICA_RHO1.md`.

## T10 — exact isolated-sphere Mie coefficients

### Scientific scope and implementation

T10 adds the exact diagonal partial-wave response of a homogeneous, lossless
fluid sphere under the project's \(e^{-i\omega t}\) convention. It changes no
Model A--D equation and does not connect the new coefficients to the global
multipolar solver. The implemented definitions are

\[
x=ka,\qquad y=x\sqrt{\widetilde\rho\widetilde\kappa},\qquad
\beta=\sqrt{\widetilde\kappa/\widetilde\rho},
\]

\[
s_\ell^{\mathrm{Mie}}=-
\frac{\beta j_\ell(x)j_\ell'(y)-j_\ell(y)j_\ell'(x)}
{\beta h_\ell^{(1)}(x)j_\ell'(y)-j_\ell(y){h_\ell^{(1)}}'(x)},
\]

with

\[
\widetilde\kappa=1-f_0,\qquad
\widetilde\rho=\frac{2+f_1}{2(1-f_1)},\qquad
\frac{c_p}{c_0}=(\widetilde\rho\widetilde\kappa)^{-1/2}.
\]

Exactly \(f_1=1\) selects the analytic sound-hard limit

\[
s_\ell^{\mathrm{rigid}}=-\frac{j_\ell'(x)}{{h_\ell^{(1)}}'(x)}.
\]

There is no clipping or artificial large density. Exact material matching
returns bitwise zeros. The API accepts \(L_{\max}=0\), validates real finite
physical inputs, and leaves the truncation choice to the caller outside the
audited Rayleigh interval.

### Independent validation and campaign

An independent test-only \(2\times2\) boundary-condition solve covers two
materials, \(ka\in\{0.01,0.05,0.1\}\), and
\(\ell=0,\ldots,5\). The maximum coefficient difference was
\(9.359714502574357\times10^{-20}\) in absolute value. The largest relative
difference, \(3.1694037722230735\times10^{-12}\), occurs in the nearly
matched, poorly scaled monopole at \(ka=0.01\). The maximum boundary-condition
residual in the 2,424-row production campaign was
\(2.2204460492505655\times10^{-16}\), and the maximum lossless-unitarity
defect \(|\operatorname{Re}s_\ell+|s_\ell|^2|\) was
\(1.3234889800848443\times10^{-23}\).

The campaign uses \(f_0=0\),
\(f_1\in\{0.1,0.4,0.8,1\}\), 101 logarithmic values over
\(10^{-3}\le ka\le0.1\), and \(\ell=0,\ldots,5\). It yields 2,424 validation
rows and 24 summary rows. At \(ka=0.1\), the complex relative dipole errors
are:

| \(f_1\) | \(\varepsilon_{s_1}\) |
|---:|---:|
| 0.1 | 0.0018021200413863723 |
| 0.4 | 0.0012045872800514213 |
| 0.8 | 0.00042193608732433199 |
| 1.0, rigid | 0.0030243505842077237 |

Every positive-order asymptotic relative-error slope is numerically 2 to the
reported precision. The \(f_0=0\) Rayleigh monopole vanishes while its exact
dynamic correction need not; the CSV therefore retains absolute error and
applicability flags and documents that this channel is inactive in the
current nodal symmetry. Increasing finite densities
\(10^8,10^{10},10^{12}\) approached the direct rigid result with relative
vector errors \(1.6337683061174943\times10^{-3}\),
\(2.0490033945573912\times10^{-4}\), and
\(1.474460765039206\times10^{-5}\).

### Files, verification, and hashes

Created:

```text
src/acoustic_ms/mie_scattering.py
tests/test_mie_scattering.py
tests/test_t10_artifacts.py
scripts/analyze_t10_mie_rayleigh.py
results/data/t10_mie_rayleigh_validation.csv
results/data/t10_mie_rayleigh_summary.csv
results/figures/t10_mie_rayleigh_error.png
TAREFA_T10_COEFICIENTES_EXATOS_MIE.md
```

Updated:

```text
src/acoustic_ms/__init__.py
README.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
```

Verification uses Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib
3.11.1. The final suite reports 258 tests passing with warnings promoted to
errors. Two consecutive campaign executions are byte-identical. T10 hashes:

```text
a4b8c3c496da33b899516dab912adc2102b84c9ea070bb62da2c5d10b009ca29  results/data/t10_mie_rayleigh_validation.csv
6a240084e1e108a2cb31efb749abf8f7e10e2ef8cf9a419a00b1c0407057a3e3  results/data/t10_mie_rayleigh_summary.csv
96a23e87c3504d697f3d2a16f920f6a98df42c78bf0006098a57664446caacc9  results/figures/t10_mie_rayleigh_error.png
```

The figure was inspected for readable logarithmic axes, four distinguishable
contrasts, all five positive multipole orders, unclipped labels, and absence
of physical `NaN` or `inf` values.

### Interpretation and remaining limitations

\[
\boxed{
\text{exact Mie coefficients}
\ne
\text{complete collective force}
}
\]

T10 validates an isolated ideal-fluid-sphere T-matrix, not Silva--Bruus or
Model D as complete force theories. The campaign is limited to
\(10^{-3}\le ka\le0.1\), real lossless properties, and \(\ell\le5\). It
does not include absorption, viscosity, elastic solids, walls, nonspherical
particles, or `scattered--scattered` force terms. Global Mie integration and
the complete-force extension belong to T11.
## T11 — complete-reference Model E

### Scope and implementation

T11 starts from \`153403c59571ad081098248165b5b184c3721179\`
with 258 baseline tests and leaves every Model A--D module and artifact
unchanged. Model E is

\[
\boxed{\text{exact Mie}+\text{global multiple scattering}+\text{complete multipolar force}.}
\]

The effective incident and scattered coefficients satisfy

\[
(I-UD)b=a,\qquad d=Db.
\]

Production solves the first equation directly with \`numpy.linalg.solve\`; an
independent test solves \((I-DU)d=Da\). Planar symmetry retains \(n+m\) odd,
while all public arrays use the complete \(n^2+n+m\) ordering and exact zeros
for inactive modes.

The force coupling is

\[
\Gamma_n=s_n+s_{n+1}^*+2s_ns_{n+1}^*,
\]

with the complete transverse and longitudinal formulas documented in
\`docs/CONVENTIONS.md\`. Their prefactors include
\(E_{\mathrm{LAS}}=2E_0\). With \(c=b-a\),

\[
F_{\mathrm{int}}=\mathcal F[b]-\mathcal F[a]
=F_{\mathrm{ext-sc}}+F_{\mathrm{ss}},
\qquad
F_{\mathrm{ss}}=\mathcal F[c].
\]

The recoil term inside \(\Gamma_n\) and the field channel \(\mathcal F[c]\)
are distinct.

### Files

Created:

\`\`\`text
src/acoustic_ms/mie_multiparticle.py
src/acoustic_ms/complete_force.py
src/acoustic_ms/model_e.py
tests/test_mie_multiparticle.py
tests/test_complete_force.py
tests/test_model_e.py
tests/test_t11_artifacts.py
scripts/t11_stress_oracle.py
scripts/analyze_t11_model_e.py
results/data/t11_model_e_convergence.csv
results/data/t11_force_oracle.csv
results/data/t11_force_decomposition.csv
results/figures/t11_model_e_validation.png
TAREFA_T11_MODELO_E_REFERENCIA_COMPLETA.md
\`\`\`

Updated:

\`\`\`text
src/acoustic_ms/__init__.py
README.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
\`\`\`

The specification `PROMPT_T11_MODELO_E_REFERENCIA_COMPLETA.md` was included in
the T11 commit together with the executed task record; it is not untracked.

### Independent validation

The stress-tensor oracle independently reconstructs the local regular plus
outgoing field and integrates

\[
\overline S/E_0=-(|g|^2-|\psi|^2)I+2\operatorname{Re}(gg^\dagger).
\]

At \(L_{\max}=4\), two particles, two radii
\(R/a\in\{1.01,1.04\}\), and angular grids \(24\times48\) and
\(32\times64\), the maximum relative component error for resolved forces is
\(6.297712539249374\times10^{-15}\). The result is independent of radius and
resolution at the reported scale. The Rayleigh cross-channel test with
\(s_0=s_2=0\) and \(s_1=i f_1(ka)^3/6\) reproduces the approved Model-D
\(L=1\) force. Full and planar bases agree at approximately
\(4.3\times10^{-16}\) in the audited dimer, and the complete tests cover
trimers and quartet geometry as well.

### Compact campaign

All cases use \(a=E_0=1\) and \(f_0=0\). Coordinates are recorded at full
precision in \`t11_model_e_convergence.csv\`.

| case | \(N\) | \(ka\) | \(f_1\) | \(d_{\min}/a\) | final \(L\) | RMS total | RMS ext-sc | RMS ss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dimer axis | 2 | 0.1 | 0.8 | 2.5 | 9 | 0.35488509499395798 | 0.33265034508023195 | 0.022234749913726034 |
| dimer diagonal | 2 | 0.05 | 0.4 | 4.0 | 9 | 0.011957953437788547 | 0.011905867214909946 | 0.000052086222878601739 |
| dimer rigid | 2 | 0.1 | 1.0 | 3.0 | 9 | 0.25220422902022588 | 0.24342909911142702 | 0.008775129908798876 |
| trimer equilateral | 3 | 0.1 | 0.8 | 3.0 | 9 | 0.28851299734919822 | 0.27523599980056734 | 0.01327699754863086 |
| trimer scalene | 3 | 0.1 | 0.8 | 2.7 | 9 | 0.29442036247659009 | 0.27961658721656057 | 0.014961914394038276 |
| quartet irregular | 4 | 0.1 | 0.8 | 2.8 | 9 | 0.25007642652434453 | 0.23873909624244133 | 0.011452274463533961 |

The RMS values are normalized by \(a^2E_0\). Convergence requires two
successive applicable changes at or below \(10^{-5}\) for each channel:

| case | total | interaction | ext-sc | ss |
|---|---:|---:|---:|---:|
| dimer axis | unconfirmed | unconfirmed | unconfirmed | unconfirmed |
| dimer diagonal | 7 | 7 | 6 | 9 |
| dimer rigid | 9 | 9 | 8 | unconfirmed |
| trimer equilateral | 9 | 9 | 8 | unconfirmed |
| trimer scalene | unconfirmed | unconfirmed | 9 | unconfirmed |
| quartet irregular | 9 | 9 | 8 | unconfirmed |

Thus a converged total does not imply convergence of the smaller
scattered--scattered term. Unconfirmed cases are not labeled divergent. The
largest system residual is \(9.575091508106486\times10^{-5}\), in the axial
dimer at \(L=9\), and the largest condition number is
\(5.906478133855654\times10^{24}\). The raw diagnostics are retained without
regularization, clipping, or substitution by the alternative scattered-field
system.

### Verification, determinism, and environment

Commands executed include:

\`\`\`bash
.venv/bin/python -m pip install -e ".[dev,plot]"
.venv/bin/python scripts/analyze_t11_model_e.py
.venv/bin/python scripts/analyze_t11_model_e.py
.venv/bin/python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
git diff --name-only
sha256sum results/data/t11_*.csv results/figures/t11_model_e_validation.png
\`\`\`

The final suite reports **295 passed** with warnings treated as errors. The
campaign has 48 convergence rows, 72 oracle rows, and 16 per-particle
decomposition rows; no physical \`NaN\` or \`inf\` is present. Two executions
in the same environment are byte-identical. The environment is Python 3.12.3,
NumPy 2.5.1, SciPy 1.18.0, and Matplotlib 3.11.1.

Official T11 hashes:

\`\`\`text
f017993c893a6a1d8db5161007ba8361dc55770922849101e8ac193d45ccf893  results/data/t11_model_e_convergence.csv
d1a8e89c62a248ac339d5f8c1e51c35b30651034460ca5b2bed0419d05f585fb  results/data/t11_force_oracle.csv
6c2b83134b306ab6b113e058439348a9de21f5f735a468a891e0868ad6f53986  results/data/t11_force_decomposition.csv
b96f9e089a833b175bea0621fef6b708c66f73a47f523655300c24a67ac9301f  results/figures/t11_model_e_validation.png
\`\`\`

All 35 pre-T11 result artifacts match the baseline manifest. The figure was
visually inspected: all three panels are populated, labels and legends are
readable, logarithmic scales and the \(10^{-5}\) line are visible, and
unconfirmed interaction cases are marked at their final order.

### Remaining limitations

Model E is limited to identical, lossless fluid spheres fixed in the nodal
plane of an ideal unbounded fluid. It includes neither viscosity, streaming,
walls, elasticity, absorption, nonspherical particles, nor dynamics. The
compact campaign is not the T12 transferability study:

\[
\boxed{\text{internal Model-E convergence}\ne
\text{validation of the }\rho_1\text{ thresholds}.}
\]

No T12 sentinel was evaluated, no T08 threshold was recalibrated, and the
T13--T14 holdout remains unopened.

## T11.1 numerical stabilization of Model E

### Scope and derivation

T11.1 changes the numerical path only. The exact T10 coefficients, translation
operator, complete force, four force channels, planar symmetry, six physical
cases, and \(L_{\max}=2,\ldots,9\) campaign are unchanged. With
\(S=D^{1/2}\) taken from the principal NumPy complex square root, production
solves

\[
(I-SUS)q=Sa,
\]

then reconstructs

\[
d=Sq,
\qquad
b=a+Ud.
\]

No \(S^{-1}\), inverse, pseudoinverse, least-squares solve, or magnitude
threshold is used. The legacy system \(A_b=I-UD\) and scattered system
\(A_d=I-DU\) remain explicit diagnostics. The legacy names system_matrix,
right_hand_side, condition_number, and residual_relative still denote
\(A_b\), \(a\), \(\kappa(A_b)\), and the legacy residual, respectively.

### Files

Created:

    scripts/analyze_t11_1_model_e_stability.py
    tests/test_t11_1_stability.py
    tests/test_t11_1_artifacts.py
    results/data/t11_1_solver_stability.csv
    results/data/t11_1_high_precision_oracle.csv
    results/figures/t11_1_model_e_stability.png
    TAREFA_T11_1_ESTABILIZACAO_NUMERICA_MODELO_E.md

Updated:

    src/acoustic_ms/mie_multiparticle.py
    src/acoustic_ms/__init__.py
    scripts/analyze_t11_model_e.py
    results/data/t11_model_e_convergence.csv
    results/data/t11_force_oracle.csv
    results/data/t11_force_decomposition.csv
    results/figures/t11_model_e_validation.png
    pyproject.toml
    README.md
    TASKS.md
    docs/CONVENTIONS.md
    docs/DECISIONS.md
    docs/HANDOFF.md

The user input PROMPT_T11_1_ESTABILIZACAO_NUMERICA_MODELO_E.md remains
untracked, unchanged, and outside the T11.1 change set.

### Conditioning and physical closures

Across all 48 rows:

| diagnostic | maximum |
|---|---:|
| \(\kappa(A_b)\) | \(5.906478133855654\times10^{24}\) |
| \(\kappa(A_d)\) | \(6.561723110752275\times10^{24}\) |
| \(\kappa(A_q)\) | \(1.1329535742333885\) |
| legacy residual | \(6.305093914217009\times10^{-5}\) |
| balanced backward error | \(5.344776630489665\times10^{-17}\) |
| \(r_b\) | \(2.0138083755836\times10^{-17}\) |
| \(r_d\) | \(3.0898191582801966\times10^{-16}\) |

Thus the raw condition number diagnoses coefficient scaling, not physical
divergence. It is also distinct from multipole convergence and from empirical
Silva–Bruus validity.

The audit solves the scattered formulation directly and applies one explicit
residual-refinement step. Maximum resolved force discrepancies between that
audit and production are:

| channel | maximum relative discrepancy |
|---|---:|
| total | \(7.892084928414508\times10^{-16}\) |
| interaction | \(7.892084928414508\times10^{-16}\) |
| external–scattered | \(7.48946401245463\times10^{-16}\) |
| scattered–scattered | \(2.1188592886720627\times10^{-15}\) |

The unrefined legacy solve remains diagnostic and differs from production by at
most \(1.3560036588172588\times10^{-7}\) among the four resolved force
channels at high order.

### High-precision linear oracle

The dimer_axis and trimer_scalene sentinels at \(L_{\max}=9\) use mpmath 1.3.0
with 80 decimal digits. The official complex128 \(A_q\) and \(Sa\) are
converted element by element, so this is an oracle for the linear solve rather
than an independent arbitrary-precision Mie or translation implementation.

| quantity | maximum relative discrepancy |
|---|---:|
| \(q\) | \(2.510742075062353\times10^{-16}\) |
| \(d\) | \(2.3455518820419097\times10^{-16}\) |
| \(b\) | \(2.407681443774223\times10^{-16}\) |
| total force | \(3.128402140403269\times10^{-16}\) |
| interaction force | \(3.128402140403269\times10^{-16}\) |
| external–scattered force | \(3.3375083154927963\times10^{-16}\) |
| scattered–scattered force | \(2.4654006544895997\times10^{-16}\) |

The independent surface-stress oracle remains satisfied; its maximum resolved
relative component error after regeneration is
\(6.4551553527306055\times10^{-15}\).

### Effect on T11 results

Relative to the official pre-stabilization T11 artifacts, the maximum changes
of RMS force channels are:

| channel | maximum relative change |
|---|---:|
| total | \(2.7919410267986013\times10^{-8}\) |
| interaction | \(2.7919410267986013\times10^{-8}\) |
| external–scattered | \(2.5925511856521585\times10^{-8}\) |
| scattered–scattered | \(6.92534995800785\times10^{-8}\) |

No convergence classification changed. The differences are numerical effects
of stabilizing the same equations, not a physical model change. The old T11
hashes in the preceding historical section are superseded by the stabilized
official hashes below.

### Verification and environment

Commands executed include:

    .venv/bin/python -m pip install -e ".[dev,plot]"
    .venv/bin/python scripts/analyze_t11_model_e.py
    .venv/bin/python scripts/analyze_t11_1_model_e_stability.py
    .venv/bin/python scripts/analyze_t11_model_e.py
    .venv/bin/python scripts/analyze_t11_1_model_e_stability.py
    .venv/bin/python -m pytest -q -W error
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only
    sha256sum results/data/t11*.csv results/figures/t11*.png

The final suite reports **307 passed** with warnings treated as errors.
The environment is Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0,
Matplotlib 3.11.1, and mpmath 1.3.0.

Two consecutive executions in this environment produced identical hashes:

    99805905800b155c7ff7d278989ef55e42f424363b5d410765959e1e221b1e6e  results/data/t11_model_e_convergence.csv
    c0738004a2fa2de8ee081e633beb8eab725c71161e6eec1111c0f2b37fc6ea9e  results/data/t11_force_oracle.csv
    8c2b83134b306ab6b113e058439348a9de21f5f735a468a891e0868ad6f53986  results/data/t11_force_decomposition.csv
    afcee969e41cb1e0d58d1a2ded3bc091b7809c953e59299b2773d00118eb2217  results/figures/t11_model_e_validation.png
    51f34cf9f89612615b1d31e8a458c14590a93b394bfd072897d43a905f8e6693  results/data/t11_1_solver_stability.csv
    eaf8c7688ad4607d6f1119551dd869ddbc7e4b5dd994bd0a6708b2b31fcd2ff7  results/data/t11_1_high_precision_oracle.csv
    06ac914313d5454b3e455b2dd511de6544da9fdf3f01e988eb349d48adea7420  results/figures/t11_1_model_e_stability.png

All pre-existing T01–T10 artifacts remain byte-identical to the initial
manifest. Safe raster inspection confirmed RGB images of 3300 by 1056 and
3410 by 1056 pixels, three populated panels in each, nonempty color content,
clear outer margins, and no invalid pixel values. Direct image export to the
conversation was blocked by the sandbox privacy policy.

### Limitations and next-step gate

Model E remains restricted to identical lossless fluid spheres fixed in the
nodal plane of an ideal unbounded fluid. It excludes viscosity, streaming,
walls, elasticity, absorption, nonspherical particles, and dynamics. The
high-precision oracle validates the linear system built in complex128, not the
full Mie and translation pipeline at arbitrary precision.

No T12 sentinel was evaluated, no \(\rho_1\) threshold was recalibrated, and
the T13–T14 holdout remains unopened. With all numerical acceptance criteria
satisfied, the solver is technically ready for a separately authorized T12;
this task does not start it.

## T12 — preregistered Model-E sentinel audit of \(\rho_1\)

### Scope and files

T12 evaluates exactly 28 preregistered T08 calibration cases: four \(\rho_1\)
bands in each of seven \((N,\text{family})\) strata with \(N\leq4\). It does
not evaluate the \(N=6,10\) holdout, recalibrate the frozen law, alter a force
model, or regenerate an earlier artifact.

Created:

    src/acoustic_ms/model_e_comparison.py
    scripts/analyze_t12_model_e_sentinels.py
    tests/test_t12_sentinels.py
    tests/test_t12_artifacts.py
    results/data/t12_sentinel_manifest.csv
    results/data/t12_model_e_convergence.csv
    results/data/t12_model_comparison.csv
    results/data/t12_threshold_audit.csv
    results/figures/t12_model_e_sentinel_audit.png
    TAREFA_T12_SENTINELAS_MODELO_E_CRITERIO_RHO1.md

Updated:

    src/acoustic_ms/__init__.py
    README.md
    TASKS.md
    docs/CONVENTIONS.md
    docs/DECISIONS.md
    docs/HANDOFF.md

The user-provided prompt files remain untracked and outside the T12 change
set.

### Reference, convergence, and metrics

The principal reference is the complete three-dimensional interaction force,

\[
\mathbf F^E=\mathbf F^E_{\mathrm{int}}
=\mathbf F^E_{\mathrm{ext-sc}}+\mathbf F^E_{\mathrm{ss}}.
\]

Models A and D are padded with \(F_z=0\); the E component \(F_z\) is retained.
All force amplitudes use

\[
\mathcal R(\mathbf F)=
\left[\frac1N\sum_i\lVert\mathbf F_i\rVert_2^2\right]^{1/2}.
\]

The audit verifies the signed-vector identity

\[
\mathbf F^E-\mathbf F^A=
(\mathbf F^D-\mathbf F^A)
+(\mathbf F^E_{\mathrm{ext-sc}}-\mathbf F^D)
+\mathbf F^E_{\mathrm{ss}}.
\]

The 28 frozen A and D vectors were reproduced through public APIs with
`rtol=5e-12, atol=5e-14` before E was evaluated. Each E case ran from
\(L_{\max}=2\) through at least 5 and stopped only after two applicable
successive changes no larger than \(10^{-5}\) in all four force channels, or
at the cap \(L_{\max}=13\). The raw table contains 289 case-order rows.
The final suite reports **328 passed**, with warnings treated as errors.

Interaction converged for 22/28 cases. All four channels converged for 18/28.
The six interaction-unconfirmed cases, all retained at \(L_{\max}=13\), are:

- `n2_pair_f1.0_d2.1`;
- `n3_compact_f0.8_d2.1`;
- `n3_irregular_f1.0_d2.1`;
- `n3_linear_f1.0_d2.1`;
- `n4_irregular_f0.8_d2.1`;
- `n4_linear_f0.8_d2.1`.

Four additional cases confirm interaction but not scattered–scattered and are
therefore ineligible for the mechanism decomposition. They are classified as
`unconfirmed`, never divergent.

### Numerical diagnostics

Across all calculated orders:

| diagnostic | maximum |
|---|---:|
| \(\kappa_2(A_q)\) | 1.8619902606818648 |
| balanced backward error | \(1.0426962134786671\times10^{-16}\) |
| physical closure error | \(6.912079812246232\times10^{-16}\) |
| force-channel decomposition residual | \(1.322894305049479\times10^{-16}\) |

The maximum relative A–D–E identity residual is
\(3.175205390812495\times10^{-14}\), and the maximum interaction \(|F_z|\)
is exactly zero. Every solve used `balanced_sqrt`; all numerical diagnostics
satisfy their computational tolerances.

### Frozen prediction and threshold audit

The unchanged prediction is

\[
\widehat\varepsilon_A=2.6353684041458636\,
\rho_1^{1.1088518115798773}.
\]

Among the 22 threshold-eligible sentinels:

| tolerance | predicted safe | observed safe | false safe | false unsafe | worst predicted-safe error |
|---:|---:|---:|---:|---:|---:|
| 1% | 7 | 7 | 0 | 0 | 0.005814217887209692 |
| 5% | 14 | 14 | 0 | 0 | 0.049538760108543481 |
| 10% | 21 | 20 | 1 | 0 | 0.12057318984999543 |

The 10% false-safe case is `n2_pair_f0.8_d2.5`. There is no false safe at
5%. Stratified records by \(N\) and family are preserved in the audit CSV.

Frozen-law performance over 22 applicable predictions is:

| metric | value |
|---|---:|
| log-space RMSE | 0.7810049747126869 |
| median multiplicative factor | 1.4773270018337583 |
| 90th-percentile factor | 2.5519818159527858 |
| maximum factor | 15.157639966319506 |
| fraction within factor 2 | 0.7272727272727273 |
| Spearman | 0.9503105590062113 |

The largest applicable normalized amplitudes are 0.11598562790520174 for
\(X_{D-A}\), 0.005971980174515004 for \(X_{\mathrm{Mie/ext-sc}}\), and
0.11924421121442162 for \(X_{\mathrm{ss}}\). They are not additive fractions.

### Gate result

`t12_gate_supported=false` and the recommendation is `NO-GO_T13`. Three gate
conditions failed: interaction coverage is 78.57%, below 80%; log-space RMSE
exceeds \(\ln2\); and only 72.73% of applicable predictions are within a
factor 2. Numerical diagnostics pass, the 5% threshold has predicted-safe
coverage, and it has no false safe. The failed scientific gate does not negate
computational completion and no parameter was changed to force approval.

### Determinism, environment, and commands

The verification environment is Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0,
and Matplotlib 3.11.1. The campaign and two explicit `--analyze-only` runs use
no random state and produced byte-identical derived artifacts.

Direct image display was unavailable because the sandbox could not configure
its loopback device. Safe raster inspection found a finite 2750×2090 RGB
image, 1,953 distinct colors, populated content in all four panels, and clear
outer margins. No invalid pixels or empty panels were detected; labels and
layout are generated with `constrained_layout=True`.

Commands include:

    .venv/bin/python -m pip install -e ".[dev,plot]"
    .venv/bin/python scripts/analyze_t12_model_e_sentinels.py
    .venv/bin/python scripts/analyze_t12_model_e_sentinels.py --analyze-only
    .venv/bin/python scripts/analyze_t12_model_e_sentinels.py --analyze-only
    .venv/bin/python -m pytest -q
    .venv/bin/python -m pytest -q -W error
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only
    sha256sum results/data/t12_*.csv results/figures/t12_model_e_sentinel_audit.png

Official T12 artifact hashes in this environment are:

    a46cf99bee4802ce42e29d5a9970c9fb8da7ae63940fa772ee2e9dc5f77befe1  results/data/t12_sentinel_manifest.csv
    a1cf482541c5d95fc5145d8a01a0e69c39fae737070094e0e03071175aaf8524  results/data/t12_model_e_convergence.csv
    3fd672da68099a264c36497ca7b6ee548f5e6430a4d9ea0f493b1f06fdd5cf91  results/data/t12_model_comparison.csv
    a33544c8040693f6be7607c8ecce28a33c9bfd0ab58034dd98bad09b3e88a516  results/data/t12_threshold_audit.csv
    3ad7ab8569640e98e86bf405e3740630ccf043941aa3c83c7d66145ff1124059  results/figures/t12_model_e_sentinel_audit.png

The initial and final manifests of all earlier versioned result files are
identical.

### Limitations

The result applies only to identical fixed spheres in an ideal unbounded
fluid, the nodal plane, \(ka=0.1\), sampled positive contrasts, and seven fixed
calibration strata with \(N\leq4\). It does not validate the external
\(N=6,10\) holdout and is not a universal criterion:

\[
\boxed{
\text{approval in }N\leq4
\ne
\text{external validation in }N=6,10
\ne
\text{universal criterion}.
\]

T13 and T14 were not started.

## T12.1 — convergence diagnosis and frozen-rho1 failure analysis

T12.1 is a restricted follow-up to the 28-case T12 calibration audit. It did
not alter Models A--E, the exact-Mie response, the force observable, the
balanced solver, \(\rho_1\), the frozen T08 power law, or its thresholds. It
copied the versioned T12 records for \(L=2,\ldots,13\) and evaluated only ten
preregistered cases at \(L=14,\ldots,21\) as needed.

### Files

Created:

- `src/acoustic_ms/rho1_model_e_diagnostics.py`;
- `scripts/analyze_t12_1_rho1_failure.py`;
- `tests/test_t12_1_diagnostics.py`;
- `tests/test_t12_1_artifacts.py`;
- `results/data/t12_1_extended_convergence.csv`;
- `results/data/t12_1_convergence_summary.csv`;
- `results/data/t12_1_resolved_comparison.csv`;
- `results/data/t12_1_mechanism_diagnostics.csv`;
- `results/data/t12_1_predictor_diagnostics.csv`;
- `results/data/t12_1_out_of_fold_predictions.csv`;
- `results/figures/t12_1_rho1_failure_diagnostics.png`;
- `TAREFA_T12_1_DIAGNOSTICO_CONVERGENCIA_FALHA_RHO1.md`.

Updated: `src/acoustic_ms/__init__.py`, `README.md`, `TASKS.md`,
`docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and `docs/HANDOFF.md`.
The three local `PROMPT_*.md` inputs were preserved as untracked user files
and were not staged.

### Convergence extension

The raw extension has 174 deterministic case-order rows. Its first 12 records
per case are byte-derived from T12 and carry `source=t12`; newly evaluated
orders carry `source=t12_1`. The 40-row summary contains one record for every
case/channel pair. Confirmation requires two successive applicable changes no
larger than \(10^{-5}\).

All 28 interaction channels are directly confirmed. Twenty-six cases confirm
all four channels. The only unconfirmed results at the hard cap are the
scattered--scattered channels of `n2_pair_f1.0_d2.1` and
`n3_irregular_f1.0_d2.1`; both are classified `unconfirmed_at_21`, not
“divergent”. Final orders for the ten extended cases are:

| case | final \(L\) | interaction confirmation | all channels |
|---|---:|---:|:---:|
| `n2_pair_f1.0_d2.1` | 21 | 20 | no |
| `n3_compact_f0.8_d2.1` | 20 | 18 | yes |
| `n3_irregular_f1.0_d2.1` | 21 | 19 | no |
| `n3_linear_f1.0_d2.1` | 21 | 19 | yes |
| `n4_irregular_f0.8_d2.1` | 20 | 17 | yes |
| `n4_linear_f0.8_d2.1` | 21 | 17 | yes |
| `n3_compact_f0.1_d2.1` | 15 | 12 | yes |
| `n4_compact_f0.1_d2.1` | 15 | 12 | yes |
| `n4_irregular_f0.1_d2.1` | 15 | 11 | yes |
| `n4_linear_f0.1_d2.1` | 15 | 12 | yes |

The linear high-contrast trimer and quartet show oscillatory changes in the
last five evaluated orders; this is recorded descriptively and does not alter
the direct two-change rule.

Across the extension, the maximum balanced condition number is
1.8777295339624336. The maxima of balanced backward error, effective-incident
closure, scattering closure, and force decomposition residual are respectively
\(7.41903330032037\times10^{-17}\),
\(3.654410728345143\times10^{-17}\),
\(6.770687779265224\times10^{-16}\), and
\(1.322894305049479\times10^{-16}\). The maximum force-channel \(|F_z|\) is
exactly zero. Every new solve used `balanced_sqrt`, remained finite, and had
condition number below 10.

### Signed mechanism diagnostics

With

\[
C_D=F^D-F^A,\qquad C_M=F^E_{\mathrm{ext-sc}}-F^D,\qquad
C_S=F^E_{\mathrm{ss}},
\]

\[
C=F^E-F^A=C_D+C_M+C_S,
\]

T12.1 records all pairwise signed cosines, the projections
\(p_j=\langle C_j,C\rangle/\langle C,C\rangle\), amplitude ratios, and the
cancellation ratio. The projection identity \(p_D+p_M+p_S=1\) and the vector
closure are tested. Mechanism fields are inapplicable when any necessary
channel lacks direct convergence.

The three preregistered special cases show:

- `n2_pair_f1.0_d6.0`: \(\varepsilon_A^E=2.4653935235086007\times10^{-4}\),
  frozen factor 15.157639966319506, \(R(C_S)/R(C_D)=0.9734004311297458\),
  \(p_S=12.441623851320271\), and cancellation ratio 49.44646599067913;
- `n2_pair_f0.8_d2.5`: the 10% false-safe result remains
  \(\varepsilon_A^E=0.12057318984999543\), with factor 2.575360131418254,
  \(p_S=0.519679558324466\), and cancellation ratio 1.0069608311364993;
- `n2_pair_f1.0_d2.1`: after direct interaction confirmation at \(L=20\),
  the large error persists at \(\varepsilon_A^E=0.38790756344523031\) and
  frozen factor 3.6574113766652632. Its mechanism decomposition remains
  inapplicable because the scattered--scattered channel is unconfirmed at 21.

### Predictor diagnostics

P0 is the untouched frozen law

\[
\widehat\varepsilon_A=2.6353684041458636\rho_1^{1.1088518115798773}.
\]

P1--P4 are fitted inside each deterministic leave-\((N,\mathrm{family})\)-out
fold. P1 uses \(\eta\), P2 uses \(\Lambda_{\max}\), P3 recalibrates \(\rho_1\),
and P4 uses the reference-derived \(\varepsilon_A^D\). All 28 directly
confirmed interaction cases appear exactly once out of fold for each candidate.

| candidate | RMSE log | median factor | p90 factor | maximum factor | within factor 2 | Spearman |
|---|---:|---:|---:|---:|---:|---:|
| P0 frozen \(\rho_1\) | 0.8106033557995027 | 1.6049065112056535 | 2.637524883542658 | 15.157639966319513 | 0.6071428571428571 | 0.9709906951286261 |
| P1 \(\eta\) | 0.7950905911225619 | 1.311059322494296 | 2.378254419902941 | 35.404797373179854 | 0.8214285714285714 | 0.9299397920087575 |
| P2 \(\Lambda_{\max}\) | 0.6293890920247307 | 1.3769284227385432 | 1.890186928563853 | 14.554217960790472 | 0.8928571428571429 | 0.9469074986316366 |
| P3 recalibrated \(\rho_1\) | 0.6458489737104017 | 1.253161254310177 | 1.820045294720189 | 18.43960615363461 | 0.9285714285714286 | 0.9600437876299945 |
| P4 \(\varepsilon_A^D\) | 0.5610640793133741 | 1.236272533986774 | 1.382567067681843 | 15.37947584862918 | 0.9642857142857143 | 0.9945265462506840 |

P4 is not an independent predictor because it requires the Model-D reference.
The largest absolute descriptive Spearman correlations of the frozen residual
are with distance ratio (-0.8549809035587734), cancellation ratio
(-0.6772649572649572, reference-derived), \(R(C_M)/R(C_D)\)
(-0.6690598290598290, reference-derived), and \(\rho_1\)
(0.6382905982905983). These 28 sentinels do not establish a universal
mechanistic criterion.

All interaction channels and numerical diagnostics pass. P3 has RMSE below
\(\ln2\), 92.86% of predictions within factor 2, and lies within 0.05 RMSE of
the best among P1--P3. The preregistered recommendation is therefore
`READY_T12_2_RHO1_RECALIBRATION_STUDY`. Historical T12 remains `NO-GO_T13`;
T13 and T14 were not started.

### Artifacts, determinism, and environment

The six CSV row counts are 174, 40, 28, 28, 82, and 140 for extended
convergence, convergence summary, resolved comparison, mechanisms, predictor
diagnostics, and out-of-fold predictions. Two explicit `--analyze-only` runs
preserved the raw extension and reproduced all derived files byte for byte.
The T12.1 hashes are:

    d41a956e9c58e5d49ab06f94b7574f8c9f987610223f9692e2c1b67297019e23  results/data/t12_1_extended_convergence.csv
    89addfd05f6d1c33160bd9bc1b3cbb6ff05f93a99dc1ce80c4c7d96a2184cbf0  results/data/t12_1_convergence_summary.csv
    5097cd7014bac635e09179e5bd4f49a0308dc4f6f02eb7bd76d60f18c2e89f39  results/data/t12_1_resolved_comparison.csv
    d3a968f322feddb7da9fd9f3fb470564c691b4b4db49ddca65248627e30e4334  results/data/t12_1_mechanism_diagnostics.csv
    d63c16d216a6235a41d44346e8c1f981e2ffd1ba3977c728713c8c57bd71842f  results/data/t12_1_predictor_diagnostics.csv
    7686c091a0323d011e79e50ebfe9dc096dc10240ca21c60b0749c9af488ac4d6  results/data/t12_1_out_of_fold_predictions.csv
    40bd04d570282fd7b6734d943f6f602e3ac6fd4e80f2d4132a235b16f80865dc  results/figures/t12_1_rho1_failure_diagnostics.png

The verification environment is Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0,
and Matplotlib 3.11.1. The final suite reports **354 passed**, with warnings
treated as errors. The initial and final SHA-256 manifests of all 47 earlier
versioned result files are identical.

The sandbox image viewer could not configure its loopback device. Structural
inspection verified a finite 2860×2090 RGBA raster, full image bounds, no empty
canvas, and deterministic bytes. The figure uses `constrained_layout=True`,
populates all four panels, uses log scales only for positive data, and carries
the prescribed identity, factor-two, and \(10^{-5}\) reference marks.

Commands executed include:

    .venv/bin/python -m pip install -e ".[dev,plot]"
    .venv/bin/python -m pytest -q -W error
    .venv/bin/python scripts/analyze_t12_1_rho1_failure.py
    .venv/bin/python scripts/analyze_t12_1_rho1_failure.py --analyze-only
    .venv/bin/python scripts/analyze_t12_1_rho1_failure.py --analyze-only
    sha256sum results/data/t12_1_*.csv results/figures/t12_1_rho1_failure_diagnostics.png
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only

### Limitations

The diagnosis is limited to the 28 frozen calibration sentinels with
\(N\leq4\), \(ka=0.1\), positive \(f_1\), fixed planar geometry families,
identical spheres, and the complete Model-E interaction-force observable. It
performs no new \(N=6,10\) computation and does not open the external holdout.
Two scattered--scattered channels remain unconfirmed at \(L=21\), so their
mechanism diagnostics are intentionally absent. The result is descriptive,
does not establish a universal error bound, and authorizes only a separately
specified T12.2 study.

## T12.2 — controlled recalibration of rho1 against Model E

### Scope and provenance

T12.2 reads the 28 canonical rows of
`results/data/t12_1_resolved_comparison.csv` and channel status from the
versioned T12/T12.1 convergence tables. It performs no multipolar or Model-E
solve. The domain remains \(N\in\{2,3,4\}\), \(ka=0.1\), \(f_0=0\), positive
sampled \(f_1\), fixed planar families, identical spheres, and the confirmed
complete Model-E interaction force. The seven generalization units are
`n2_pair`, `n3_compact`, `n3_irregular`, `n3_linear`, `n4_compact`,
`n4_irregular`, and `n4_linear`.

All 28 interaction/total references are confirmed. The
scattered--scattered channels of `n2_pair_f1.0_d2.1` and
`n3_irregular_f1.0_d2.1` remain `unconfirmed_at_21` and are explicitly marked;
this does not remove their confirmed interaction force from the regression.
No \(N=6,10\) case was read as a calibration row or evaluated.

Created files are:

- `src/acoustic_ms/rho1_model_e_recalibration.py`;
- `scripts/analyze_t12_2_rho1_recalibration.py`;
- `tests/test_t12_2_recalibration.py`;
- `tests/test_t12_2_artifacts.py`;
- six `results/data/t12_2_*.csv` files;
- `results/figures/t12_2_rho1_recalibration.png`;
- `TAREFA_T12_2_RECALIBRACAO_CONTROLADA_RHO1.md`.

Updated files are `src/acoustic_ms/__init__.py`, `README.md`, `TASKS.md`,
`docs/DECISIONS.md`, and `docs/HANDOFF.md`. The four local prompt inputs remain
untracked and are not part of the commit.

### Confirmatory protocol

The single candidate is

\[
\log\widehat\varepsilon_A^E=\beta_0+\beta_1\log\rho_1,
\qquad
\widehat\varepsilon_A^E=C_E\rho_1^{\alpha_E}.
\]

Each LOGO fold fits only the other six groups. The 28 strictly OOF predictions
retain canonical sentinel order. The frozen P0 law is evaluated without
refitting. No epsilon floor is used because all 28 target errors are strictly
positive; the recorded value is zero.

| held-out group | train/test | \(C_E^{(-g)}\) | \(\alpha_E^{(-g)}\) | 1% threshold | 5% threshold | 10% threshold |
|---|---:|---:|---:|---:|---:|---:|
| n2_pair | 24/4 | 8.664479283491534 | 1.276944106033447 | 0.005004909122702858 | 0.0176511631671110 | 0.0303749680210992 |
| n3_compact | 24/4 | 15.08273461506406 | 1.435586621602670 | 0.006108525564888686 | 0.0187423683582234 | 0.0303749623820695 |
| n3_irregular | 24/4 | 16.14180680161521 | 1.447883002680308 | 0.006086735301654885 | 0.0184985425625791 | 0.0298571228848490 |
| n3_linear | 24/4 | 16.03273298265274 | 1.444265785492683 | 0.006037736222652936 | 0.0184007831001668 | 0.0297349674845416 |
| n4_compact | 24/4 | 16.86608146832317 | 1.453133997997103 | 0.006015514813023603 | 0.0182088047039304 | 0.0293386808154164 |
| n4_irregular | 24/4 | 16.49085211796474 | 1.453121927799739 | 0.006109118707505293 | 0.0184923113498588 | 0.0297955949711974 |
| n4_linear | 24/4 | 16.23827532420431 | 1.449085672013724 | 0.006087479079470081 | 0.0184837428561375 | 0.0298213846896951 |

### OOF performance and safety

| model | RMSE log | MAE log | median absolute log ratio | within factor 2 | within factor 1.5 | Spearman | maximum log underestimation |
|---|---:|---:|---:|---:|---:|---:|---:|
| P0 frozen | 0.8106033557995027 | 0.6330031366128833 | 0.4726269784443000 | 0.6071428571428571 | 0.4285714285714286 | 0.9709906951286261 | 1.2967556230417525 |
| recalibrated rho1 | 0.6458489737104012 | 0.3786559305780956 | 0.2256534353283870 | 0.9285714285714286 | 0.7142857142857143 | 0.9600437876299945 | 0.7513842151522394 |

| tolerance | predicted safe | groups covered | true safe | false safe | false unsafe | worst false-safe excess |
|---:|---:|---:|---:|---:|---:|---:|
| 1% | 7 | 7 | 7 | 0 | 0 | 0 |
| 5% | 14 | 7 | 14 | 0 | 0 | 0 |
| 10% | 20 | 7 | 19 | 1 | 1 | 0.020573189849995427 |

The 10% false-safe case is `n2_pair_f0.8_d2.5`. It is retained without
clipping, exclusion, margin, or post-hoc safety factor.

### Final descriptive calibration and uncertainty

Only after freezing OOF predictions, the same model fitted to all 28 cases is

\[
\widehat\varepsilon_A^E
=14.73950709797405\rho_1^{1.4226504975598322}.
\]

Its candidate thresholds are 0.005926947606709601,
0.01837157635582504, and 0.029905042165737895 for 1%, 5%, and 10%.
The fold ranges are 8.664479283491534--16.86608146832317 for \(C_E\) and
1.276944106033447--1.453133997997103 for \(\alpha_E\).

The whole-group bootstrap used seed 1202 and produced 10,000 valid samples in
10,000 attempts. Percentile 95% intervals are:

- \(C_E\): [7.579782806369310, 47.07370137527077];
- \(\alpha_E\): [1.249649355494450, 1.728753959643017];
- 1% threshold: [0.004875978405960654, 0.007605269077613946];
- 5% threshold: [0.01714149487572852, 0.01942975087950808];
- 10% threshold: [0.02830655600345593, 0.03155362429778696].

These intervals are descriptive and do not change the gate.

### Gate decision

Nine of ten criteria pass: all predictions and fold coefficients are positive;
RMSE is below \(\ln2\); 92.86% lie within factor 2; Spearman is 0.9600; both
primary metrics improve on P0; coverage is nonempty; and integrity checks pass.
The zero-false-safe criterion fails at 10%. The exact decision is therefore:

```text
NO_GO_T13_RHO1_NOT_QUANTITATIVE
```

The result supports \(\rho_1\) as an ordinal or mechanistic indicator but does
not validate the simple power law as an autonomous quantitative criterion. No
second candidate was tried. T13 and T14 remain unopened.

### Verification, determinism, and artifacts

The final suite reports **371 passed** with warnings treated as errors. The
verification environment is Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and
Matplotlib 3.11.1. Two consecutive script executions in this environment
produce identical bytes. Different numerical/plotting environments need not
produce identical binary representations.

Artifact hashes are:

    e0c4b24078ea64d142808f64e170e346621d8252155f58f407620b972eab4b48  results/data/t12_2_logo_predictions.csv
    f75b2d04fc342c1b211bd6171ad2b2fe7668141d7787b7baf63d0966efc55e4b  results/data/t12_2_logo_fits.csv
    4a4a013e4e83eacf32bd766851349e52a10100f2fdbdb555ebca67f7da898e56  results/data/t12_2_metrics.csv
    b262f8c4576501a93e4d1dbc5a69d344f225a318965db603d55742afd7607328  results/data/t12_2_safety_audit.csv
    ea43f118469b46c41acf03b03780750fe977c79d6816eab7c044383c39398da6  results/data/t12_2_final_calibration.csv
    0059c29c652af7a3ab84d1b1469c8824b1679196a2dbab4b4877cf43d9a6bbbf  results/data/t12_2_gate.csv
    bc9eee3afd9d717cb16aeebe8a2482aece9351930fb0568e5a71f7c841654104  results/figures/t12_2_rho1_recalibration.png

The initial and final manifests of all 54 earlier versioned artifacts are
identical. The sandbox image viewer again failed while configuring loopback;
structural raster inspection verified a finite 2904×2156 RGBA image with full
bounds, all channels spanning valid ranges, four populated panels,
`constrained_layout=True`, and deterministic bytes.

Commands include:

    .venv/bin/python -m pip install -e ".[dev,plot]"
    .venv/bin/python scripts/analyze_t12_2_rho1_recalibration.py
    .venv/bin/python scripts/analyze_t12_2_rho1_recalibration.py
    .venv/bin/python -m pytest -q -W error
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only
    sha256sum results/data/t12_2_*.csv results/figures/t12_2_rho1_recalibration.png

### Limitations

This is a small, internally calibrated \(N\leq4\) data set, not independent
validation. It covers only \(ka=0.1\), \(f_0=0\), sampled positive contrasts,
fixed planar geometry families, identical spheres, and the Model-E interaction
force. It provides no result for \(N=6,10\), no universal threshold, and no
new physical descriptor. In particular,

\[
\boxed{\text{candidate calibration}\ne\text{external validation}\ne
\text{universal criterion}.}
\]

## T12.3 — grouped mechanistic validity criterion

### Scope and implementation

T12.3 reads only the 28 canonical rows of
`results/data/t12_1_resolved_comparison.csv`. It performs no Model-A, Model-D
or Model-E solve and never reads or inspects the external \(N=6,10\) holdout.
The seven frozen groups are `n2_pair`, `n3_compact`, `n3_irregular`,
`n3_linear`, `n4_compact`, `n4_irregular`, and `n4_linear`.

Created files:

- `src/acoustic_ms/mechanistic_validity.py`;
- `scripts/analyze_t12_3_mechanistic_validity.py`;
- `tests/test_t12_3_mechanistic_validity.py`;
- `tests/test_t12_3_artifacts.py`;
- eight `results/data/t12_3_*.csv` tables;
- `results/figures/t12_3_mechanistic_validity.png`;
- `TAREFA_T12_3_CRITERIO_MECANISTICO_VALIDACAO_AGRUPADA.md`.

Updated files are `src/acoustic_ms/__init__.py`, `README.md`, `TASKS.md`,
`docs/DECISIONS.md`, and `docs/HANDOFF.md`. Local `PROMPT_*.md` inputs remain
untracked and are not part of the change.

The primary candidate is

\[
\widehat\varepsilon_{M1}=C_\Lambda\Lambda_{\max}^{\alpha_\Lambda},
\qquad
\Lambda_{\max}=|f_1|\max_i\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3.
\]

The second candidate is

\[
\widehat\varepsilon_{M2}=C_{\Lambda\rho}
\Lambda_{\max}^{\alpha_\Lambda}\rho_1^{\alpha_\rho}.
\]

M1 precedes M2 because it adds one independent geometric descriptor to the
failed isolated-\(\rho_1\) criterion without searching combinations. Both use
ordinary least squares in natural-log coordinates. P0 and P3 are evaluated
with their frozen coefficients and are never refitted.

Each outer LOGO train has 24 cases and each test has four. Inside every train,
a six-group LOGO produces 24 honest predictions and

\[
s=\exp\!\left[\max_j
(\log\varepsilon_j-\log\widehat\varepsilon_j^{\mathrm{inner}})
\right],
\qquad
\widehat\varepsilon_{\mathrm{safe}}=s\widehat\varepsilon_{\mathrm{OOF}}.
\]

Safety uses strict inequalities: both predicted and observed safe mean values
strictly below \(\tau\). A synthetic test changes the held-out group's errors
by a factor \(10^{12}\) and confirms that its fit, inner predictions, margin,
and point predictions are unchanged.

### Fits, OOF metrics, and collinearity

Only after the OOF decision, the descriptive complete fits are

\[
\widehat\varepsilon_{M1}
=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\]

\[
\widehat\varepsilon_{M2}
=5.0396777007270535\,
\Lambda_{\max}^{1.2602714475189609}\rho_1^{0.13234182295409233}.
\]

| model | RMSE log | MAE log | within factor 2 | within factor 1.5 | Spearman | worst factor |
|---|---:|---:|---:|---:|---:|---:|
| P0 | 0.810603355799503 | 0.633003136612883 | 0.607142857142857 | 0.428571428571429 | 0.970990695128626 | 15.1576399663195 |
| P3 | 0.577854414091078 | 0.336406995083200 | 0.928571428571429 | 0.750000000000000 | 0.970990695128626 | 13.2499440166714 |
| M1 | 0.629389092024730 | 0.400097164805058 | 0.892857142857143 | 0.714285714285714 | 0.946907498631637 | 14.5542179607905 |
| M2 | 0.662457567037336 | 0.426097483531842 | 0.928571428571429 | 0.678571428571429 | 0.943623426382047 | 17.3299357001287 |

For M2, the full log-predictor correlation is 0.995454787228431 and the
standardized design condition is 20.95288536113825. Four of seven outer folds
have both exponents positive; three flip the sign of \(\alpha_\rho\). The
group-bootstrap 95% intervals are [-0.0786779184889365,
3.92257728418912] for \(\alpha_\Lambda\) and [-2.44090704448650,
1.36092235099522] for \(\alpha_\rho\). M2 is therefore labeled
`UNSTABLE_COLLINEARITY` under the preregistered fold-identifiability rule.

### Conservative safety audit

| model | tolerance | predicted safe | observed safe | false safe | false unsafe | safe coverage |
|---|---:|---:|---:|---:|---:|---:|
| P3 | 1% | 7 | 7 | 0 | 0 | 1.000000000000000 |
| P3 | 5% | 13 | 14 | 0 | 1 | 0.928571428571429 |
| P3 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |
| M1 | 1% | 4 | 7 | 0 | 3 | 0.571428571428571 |
| M1 | 5% | 9 | 14 | 0 | 5 | 0.642857142857143 |
| M1 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |
| M2 | 1% | 4 | 7 | 0 | 3 | 0.571428571428571 |
| M2 | 5% | 9 | 14 | 0 | 5 | 0.642857142857143 |
| M2 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |

M1 meets the 3/8/12 minima and has zero false safe. The former 10% false-safe
case `n2_pair_f0.8_d2.5` has observed error 0.12057318984999543, point M1
prediction 0.07162592847780698, outer-fold factor 2.1299122614975046, and safe
prediction 0.15255694330602437. It is now correctly unsafe.

### Gate, artifacts, and determinism

All eight M1 criteria pass. The hierarchical decision is

```text
GO_T13_VALIDATE_LAMBDA_MAX
```

M2 is diagnostic because M1 already passes; independently, M2 fails its sign
stability and incremental-value items. A GO here is a candidate-selection
result only, not external validation.

The new artifact hashes are:

    a8c081bf93c1a0d46c8cb230b5415648bafb4e2c2d85df4cd03d1f8f8a83e63a  results/data/t12_3_case_influence.csv
    e7a38c8d4148d1e4e44521be793b1a97f9a4b0bcc81e634350794f1d9e351f0a  results/data/t12_3_gate.csv
    d429a621f9e1fa99ce3a1c745b9d865b28167866fb9216c890b7234a27f91341  results/data/t12_3_group_bootstrap.csv
    4f63024eebcdbe509d796fc54856d867ee212f95d463a278b45ce44501e068f5  results/data/t12_3_logo_coefficients.csv
    8ff189eb43c265dba37f8e823ecbcd09879c5102a0deab356c51407d94ab9ee3  results/data/t12_3_metrics.csv
    ff671d2251b89ef57507919357e1df9ad450b35389d38f3074a9dcc01deb9bc8  results/data/t12_3_nested_safety_factors.csv
    c558ac04458044ee95a3a36248a0fab1d46a387be1affc1e9951d5d895a43e04  results/data/t12_3_oof_predictions.csv
    effb98332c54b1382f7fd9deee9e1e0dbc518f88ac9a747cf951fbb54658d80c  results/data/t12_3_threshold_audit.csv
    279f8848fa6ab1173c6c8e0b2eceb404c13563db19d336f34df5d80246ac31c5  results/figures/t12_3_mechanistic_validity.png

Two complete executions produced identical bytes. The bootstrap used public
seed 1203 and 10,000 valid samples in 10,000 attempts. The environment is
Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib 3.11.1. All 61
versioned earlier artifacts — 54 predating T12.2 and seven from T12.2 — match
their initial SHA-256 manifest.

Commands include:

    .venv/bin/python -m pip install -e ".[dev,plot]"
    .venv/bin/python scripts/analyze_t12_3_mechanistic_validity.py
    .venv/bin/python scripts/analyze_t12_3_mechanistic_validity.py
    .venv/bin/python -m pytest -q
    .venv/bin/python -m pytest -q -W error
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only

### Limitations

This internal candidate-selection result is limited to 28 sentinels with
\(N\leq4\), \(ka=0.1\), \(f_0=0\), positive sampled contrasts, fixed planar
families, identical spheres, and the complete Model-E interaction force. The
sample is small, M2 is strongly correlated, and the bootstrap is descriptive.
No \(N=6,10\) case, new force solve, T13, or T14 was executed.

## T13 — external validation of the frozen Lambda-max criterion

### Chronology, scope, and files

The response-blind Fase A was committed as
`af29faf89cb8f8c6883cc8bea0d44073e7caf020` and pushed before the first new
Model-E solve for \(N=6,10\). It fixed 24 IDs, the M1/P3 formulas, conservative
margins, target levels, thresholds, metrics, convergence protocol, and gate.
The three blind files remained byte-identical throughout Fase B.
The pre-T13 baseline was 388 passing tests; the published Fase-A suite had 395.

Created or updated implementation and documentation files are:

- `src/acoustic_ms/external_validation.py` and `src/acoustic_ms/__init__.py`;
- `scripts/preregister_t13_external_validation.py`,
  `scripts/run_t13_external_validation.py`, and
  `scripts/analyze_t13_external_validation.py`;
- `tests/test_t13_preregistration.py` and
  `tests/test_t13_external_validation.py`, plus the future-proof artifact filter in
  `tests/test_t12_3_artifacts.py`;
- `TAREFA_T13_VALIDACAO_EXTERNA_LAMBDA_MAX.md`, `README.md`, `TASKS.md`,
  `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and `docs/HANDOFF.md`;
- three blind CSVs, seven revealed CSVs, and
  `results/figures/t13_external_validation.png`.

No solver, scattering coefficient, force equation, T08 raw table, or T01–T12.3
artifact was changed. The local `PROMPT_*.md` inputs remain untracked.

### Frozen sample and protocol

The sample contains four target levels for each of `n6_linear`, `n6_compact`,
`n6_irregular`, `n10_linear`, `n10_compact`, and `n10_irregular`. Thus it has
12 cases for each \(N\), eight per family, six per level, and 24 total. The
physical parameters are \(a=E_0=1\), \(ka=0.1\), \(f_0=0\), positive sampled
\(f_1\), identical planar spheres, and the complete Model-E interaction force.

The frozen confirmatory law is

\[
\widehat\varepsilon_{M1}
=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\qquad
\widehat\varepsilon_{M1,\mathrm{safe}}
=2.5699703122019222\,\widehat\varepsilon_{M1}.
\]

The blind conservative safe counts were 6, 12, and 18 at 1%, 5%, and 10%.
Convergence requires two successive applicable channel changes no larger than
\(10^{-5}\), with no early stop below \(L_{\max}=5\), standard cap 13, and
interaction-only extension cap 21.

### Campaign, convergence, and diagnostics

Exactly 24 cases were solved once, producing 205 order rows and 192 long-form
particle-force rows. All 24 interaction references are confirmed, eligible,
finite, and diagnostically approved. No case required extension beyond 13.
Final orders were: three cases at 6, one at 7, two at 8, seven at 9, four at
10, two at 11, two at 12, and three at 13.

The eight-case stratified `--audit-existing` sample reproduced both values of
\(N\), all three families, all four target levels, all force channels, and the
balanced condition numbers. Across all calculated orders, the maxima were:

| diagnostic | maximum |
|---|---:|
| balanced condition number | 1.073351390634244 |
| balanced backward error | \(8.327818760262406\times10^{-17}\) |
| effective-incident closure | \(4.458332851489474\times10^{-17}\) |
| scattering closure | \(5.913315086486320\times10^{-16}\) |
| force-channel decomposition residual | \(1.111617380757604\times10^{-16}\) |
| planar \(\max|F_z|\) | 0 |

Six holdout cases lie outside the development range of
\(\Lambda_{\max}\); they remain present and are explicitly flagged rather
than discarded. The signed mechanism identity has maximum absolute closure
error \(1.0408340855860843\times10^{-17}\). Mechanism and cancellation
diagnostics do not participate in the gate.

### External metrics and strict safety audit

| model/scope | RMSE log | MAE log | factor-2 fraction | Spearman | worst factor |
|---|---:|---:|---:|---:|---:|
| M1 all | 0.465095692587546 | 0.405022711598156 | 0.875000000000000 | 0.964347826086956 | 2.493764309928004 |
| M1 \(N=6\) | 0.454863744834485 | 0.369519569189719 | 0.833333333333333 | 0.965034965034965 | 2.493764309928004 |
| M1 \(N=10\) | 0.475107335411828 | 0.440525854006593 | 0.916666666666667 | 0.958041958041958 | 2.405678357508736 |
| P3 all | 0.428876277730838 | 0.342637422617774 | 0.875000000000000 | 0.964347826086956 | 2.568595139431090 |

M1 produced zero false-safe cases at every tolerance. At 1%, 5%, and 10%,
its predicted-safe counts were 6, 12, and 18; observed-safe counts were 11,
20, and 24; false-unsafe counts were 5, 8, and 6. The worst observed errors
inside its predicted-safe regions were 0.00337799714600463,
0.0145233336544098, and 0.0367211292312633. P3 is reported transparently but
cannot alter the M1 decision.

Every sufficiency and scientific item passed. The exact result is:

```text
PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX
GO_T14_SCALE_OUT_WITH_FROZEN_LAMBDA_MAX
```

This is external evidence for the frozen M1 criterion in the sampled domain,
not a universal error theorem and not an implementation of T14.

### Artifacts, determinism, and commands

The blind SHA-256 hashes are:

```text
25d79db59d9dd6d52c5674d0a64fe2fea351cf213a0cdcd92b45845a9ecc2b38  results/data/t13_holdout_manifest.csv
581a748dca2e5d161890284fca673ed20f2a4fbcbc0ff356d5d31db6ec8ac9c2  results/data/t13_frozen_predictions.csv
eb1878e3425ede7a2b599fd20f63550d2fdb23d177264a043d443694907dc650  results/data/t13_frozen_protocol.csv
```

The revealed SHA-256 hashes are:

```text
ac73b4b12d1ab937fd39d7b62e5446e021058e21b207295e57fd4c818ae0d95d  results/data/t13_model_e_convergence.csv
95719a996df27f0b11c1828bb2589403a4d9a8b8df2447f0c37e556532d8490b  results/data/t13_forces.csv
a07dfab5386a15be32894ca69638e81fe9389f5adefca01f1d8e301e1559ae97  results/data/t13_case_summary.csv
1936ccdf7dbda9eec6e8f68fe181157f18fd43f059f0002ce3fc560075910ae1  results/data/t13_external_predictions.csv
1e2fac933d0684fcd953a389255e0ac1ef83ea9c8500ffe380f82f0e0cdb2788  results/data/t13_metrics.csv
81d1c3d0dc84e232f24b1165303a16b0035af1a336baa76905dba918478dd6ff  results/data/t13_threshold_audit.csv
8d6f8484c44c8c26499862c9572b62ae7ec770efdb5a46a232119aa5b0933761  results/data/t13_gate.csv
76dd31bf849855330ae3c61bc438c62ee38ab9be7079631c01d06b2a5c243776  results/figures/t13_external_validation.png
```

Two `--analyze-only` executions in the same environment produced identical
bytes. The figure is a finite 2904×1804 RGBA PNG; visual inspection confirmed
six populated panels, legible text, visible identities/thresholds, distinct
\(N\), families and target levels, and no destructive overlap. The verification
environment is Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib 3.11.1.

Commands executed include:

```bash
.venv/bin/python scripts/preregister_t13_external_validation.py
.venv/bin/python -m pytest -q -W error
.venv/bin/python scripts/run_t13_external_validation.py --workers 2
.venv/bin/python scripts/run_t13_external_validation.py --audit-existing
.venv/bin/python scripts/run_t13_external_validation.py --analyze-only
git diff --check
git status --short
git diff --stat
git diff --name-only
sha256sum results/data/t13_*.csv results/figures/t13_external_validation.png
```

All 70 pre-T13 result artifacts retained their original hashes. The three
blind files remained identical to the first commit. The final suite reports
405 passed with warnings treated as errors.

### Limitations

The result covers only \(N=6,10\), three fixed planar family constructions,
the sampled positive contrasts, \(ka=0.1\), \(f_0=0\), identical fixed fluid
spheres, and the complete interaction force of the present Model E. It does
not establish universality, other frequencies, negative contrasts, arbitrary
geometries, dynamics, viscosity, streaming, or walls. T14 remains a separate
future task.

## T14 — scale-out validation of frozen Lambda-max

### Chronology and implementation

The response-blind phase A was committed as
`6520173359b29cffa3a3d6432cefafcf17310f69` and pushed to `origin/main`
before the first T14 Model-E solve. It created the task specification,
`scale_out_validation.py`, three scripts, preregistration tests, and four blind
CSVs. Phase B added the single official campaign, result tests, derived
tables, figure, and this documentation. A serialization-only cache fix was
needed after the first computed case; no official raw CSV had yet been
published and no scientific quantity or blind artifact changed.

The sample has 24 cases: 12 each for \(N=15\) and \(N=28\), eight per family
and six per target level. The physical values are \(a=E_0=1\), \(ka=0.1\),
\(f_0=0\), and \(f_1=0.8\). Linear, triangular compact, and deterministic
irregular templates are centered, planar, normalized to unit minimum distance,
and scaled analytically to the four frozen \(\Lambda_{\max}\) targets.

### Campaign and numerical diagnostics

The sequential campaign (`workers=1`) produced 162 order rows and 516
particle-force rows. All 24 cases confirmed total, interaction,
external–scattered, and scattered–scattered channels and were eligible. Final
orders were five cases at L=6, seven at L=7, six at L=8, three at L=9, one at
L=10, and two at L=11. No solve exceeded L=11.

Across all order rows, the maximum balanced condition number was
1.0685477306640647, balanced backward error
\(4.675035074483183\times10^{-17}\), incident closure
\(4.395906289090826\times10^{-17}\), scattering closure
\(4.33178657257686\times10^{-16}\), and force decomposition residual
\(1.232528360511498\times10^{-16}\). The maximum \(|F_z|\) was zero. The
largest final system dimension was 1848, recorded peak process memory was
1,032,216 KiB, and recorded cumulative solve time was 2963.8310701187 s.

### Frozen-prediction results

| model/scope | RMSE log | MAE log | factor-2 fraction | Spearman | worst factor |
|---|---:|---:|---:|---:|---:|
| M1 all | 0.343011370242051 | 0.298275154329041 | 1.0 | 0.923685071658958 | 1.770801654598976 |
| M1 \(N=15\) | 0.370986850641528 | 0.326241864212104 | 1.0 | 0.945531595595981 | 1.770801654598976 |
| M1 \(N=28\) | 0.312541768219449 | 0.270308444445977 | 1.0 | 0.936642944221318 | 1.765778382477152 |
| P3 all | 0.314409582355776 | 0.255714075434696 | 1.0 | 0.904347826086957 | 1.918738739670054 |
| P3 \(N=15\) | 0.315406152386347 | 0.268024171274447 | 1.0 | 0.895104895104895 | 1.847431421543507 |
| P3 \(N=28\) | 0.313409843481514 | 0.243403979594945 | 1.0 | 0.916083916083916 | 1.918738739670054 |

At 1%, 5%, and 10%, M1 retained blind predicted-safe counts 6, 12, and 18,
with observed-safe counts 8, 18, and 24. It produced zero false-safe cases;
the worst errors within predicted-safe regions were 0.002598850242074866,
0.01121330957997045, and 0.02564715331758862. All sufficiency and scientific
criteria passed. The literal decision is:

```text
PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX
GO_T15_SYNTHESIS_AND_MANUSCRIPT
```

At matched \(\Lambda_{\max}\), the 12 ratios
\(\varepsilon_{N=28}/\varepsilon_{N=15}\) ranged from
0.729329016341299 to 1.002844791946528. This is descriptive evidence inside
the fixed families, not a new fitted size correction.

### Artifacts, determinism, and commands

The phase-A hashes are:

```text
7ab7e9ee0965ce560dbac5d33a48dad7e3cde1ce83bee73a23a38879a2c60682  results/data/t14_scale_manifest.csv
cb12db65123c694f8e3936a03ca8466829a254e7bcf493e04b22b4dbce39d86c  results/data/t14_frozen_predictions.csv
1bfc1a499fc28ffc47d99c3c9566d8360a9ed3ede1b9a918011132f2a5bce2c6  results/data/t14_frozen_protocol.csv
9b85642f889e92b626da1e15cec39954e2ebbf07419801597c9d874184355200  results/data/t14_prior_artifact_hashes.csv
```

The phase-B hashes are recorded in the final T14 commit and were reproduced by
two analysis-only executions in the same environment. The verification
environment was Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0, and Matplotlib
3.11.1. Commands included:

```bash
.venv/bin/python scripts/preregister_t14_scale_out.py
.venv/bin/python scripts/run_t14_scale_out.py --workers 1
.venv/bin/python scripts/run_t14_scale_out.py --audit-existing
.venv/bin/python scripts/run_t14_scale_out.py --analyze-only
.venv/bin/python -m pytest -q -W error
git diff --check
sha256sum results/data/t14_*.csv results/figures/t14_scale_out_validation.png
```

The 81 prior result artifacts retained path, size, and SHA-256. The four blind
files stayed identical to phase A. Visual inspection confirmed six populated,
legible panels, distinct sizes and families, visible identity/tolerance guides,
and no NaN, infinity, clipped axes, or destructive legend overlap.

### Limitations

The conclusion is restricted to \(N=15,28\), the three prescribed uniformly
scaled planar families, \(ka=0.1\), \(f_0=0\), \(f_1=0.8\), identical fixed
fluid spheres, and the complete Model-E interaction force. It does not cover
negative contrast, other frequencies, arbitrary geometries, viscosity,
streaming, walls, or dynamics. `GO_T15_SYNTHESIS_AND_MANUSCRIPT` is the next
gate; T15 was not implemented by T14.

### Official T14 artifact hashes

```text
c8b70a9ad2f87aeb231c7deab9f5fe43d6308baa6a2013cda5b0a4b97a4fa164  results/data/t14_model_e_convergence.csv
301aaee1c464d2cb71928d6263a08610b5ae808254c95aa4aaaa2187c1d83186  results/data/t14_forces.csv
018402f57571ae8295f7218346403f9f5801ec51ca8053f1a8e47d6c25a2e9f2  results/data/t14_case_summary.csv
dd0334bdb00409bbfb229e01f3c31a2b2cba230e2b9bbe6a3445256e65f8b4dc  results/data/t14_scale_predictions.csv
690c0aad9b44d2ab619a56ea74042741b49575bfbc780605edffc8fa86cf2e04  results/data/t14_metrics.csv
e07273b7b973f03fec45df102b906390abb4ff6e0c6952abb679daccc5450ad8  results/data/t14_threshold_audit.csv
cdb800c2f99b0d085d168a3124c48771d1a8e878cb02284f1bd8e89fbd83ba1d  results/data/t14_matched_scale_pairs.csv
aad0a941997abd2898c538113137ae43c2fac3e3882e096c909eb99b48e26d7b  results/data/t14_performance.csv
715d51a1c3e82242f295de5b184f7d39d3f6297a09565f0dee83f7201221a315  results/data/t14_gate.csv
d54185efa1279d61a769e41f0d53cd165ba853da89ffc0a2e28ffb9e510123ec  results/figures/t14_scale_out_validation.png
```

### Files created or updated

- `src/acoustic_ms/scale_out_validation.py` and public exports in
  `src/acoustic_ms/__init__.py`;
- `scripts/preregister_t14_scale_out.py`, `scripts/run_t14_scale_out.py`, and
  `scripts/analyze_t14_scale_out.py`;
- `tests/test_t14_preregistration.py`, `tests/test_t14_scale_out_results.py`,
  and the future-proof T14 exclusion in `tests/test_t12_3_artifacts.py`;
- `TAREFA_T14_SCALE_OUT_LAMBDA_MAX.md`, `README.md`, `TASKS.md`,
  `docs/CONVENTIONS.md`, `docs/DECISIONS.md`, and `docs/HANDOFF.md`;
- four blind CSVs, the raw convergence CSV, eight derived CSVs, and
  `results/figures/t14_scale_out_validation.png`.

The phase-A suite reported 425 passing tests. The final suite reported 434
passing tests with warnings treated as errors. No solver, scattering
coefficient, force equation, prior result artifact, or local untracked prompt
was modified or committed.

## T14.1 — frozen Lambda-max confirmation at N=45 and N=105

### Chronology, scope, and files

The response-blind phase was committed and pushed as
`538142b638dd59768d26bd16809b1def83bfdf8c` at
2026-07-31T19:53:48-03:00 with message
`chore: preregister T14.1 large-N confirmation`. The official raw campaign was
published afterward, at 2026-08-01T02:23:05-03:00. The runner verified the
phase-A artifact and code hashes before every solve and before the independent
audit.

Phase A created `src/acoustic_ms/large_n_validation.py`, public exports,
three T14.1 scripts, `tests/test_t14_1_preregistration.py`, the tracked task
specification, and five blind CSVs. Phase B created
`tests/test_t14_1_large_n_results.py`, the raw convergence table, nine derived
CSVs, and the six-panel figure, and updated `README.md`, `TASKS.md`, and the
three project documents. No solver, force equation, scattering module, prior
result artifact, or untracked local prompt was modified.

The 24 IDs are ordered as \(N=45\) followed by \(N=105\); within each size the
families are `linear`, `compact`, and `irregular`, and each family uses levels
1--4. Thus there are 12 cases per size, eight per family, and six per level.
The fixed parameters are \(a=E_0=1\), \(ka=0.1\), \(f_0=0\), and \(f_1=0.8\).
The four targets are
\(0.0031111241226691642\), \(0.011108933664494051\),
\(0.025457132710914911\), and \(0.065350897425260762\).

### Eligibility and convergence

All 24 cases are eligible; no case was imputed or discarded. The raw CSV has
160 case-order rows and every final record has stop reason
`all_channels_confirmed`. Final orders were:

| \(N\) | family | levels 1, 2, 3, 4 |
|---:|---|---|
| 45 | linear | 7, 8, 9, 11 |
| 45 | compact | 6, 6, 7, 8 |
| 45 | irregular | 6, 7, 8, 9 |
| 105 | linear | 7, 8, 9, 11 |
| 105 | compact | 6, 6, 7, 8 |
| 105 | irregular | 6, 7, 8, 9 |

The distribution is six cases each at orders 6, 7, and 8, four at order 9,
and two at order 11. No order above 11 was required. The maximum balanced
condition number was 1.0686882116757186; balanced backward error
2.6187917511006806e-17; incident closure 3.675972359325286e-17; scattering
closure 2.4097473269235716e-16; force-decomposition residual
4.416726258397735e-17; and \(\max|F_z|=0\). The largest separately recorded
unbalanced physical residual was 0.13760124781707805 at the highest-order
linear case; it is not a gate quantity under the frozen protocol and did not
replace the well-scaled backward-error checks.

### Frozen predictions and conservative audit

| model/scope | points | RMSE log | MAE log | median factor | p90 factor | max factor | factor-2 fraction | Spearman |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| M1 all | 24 | 0.30428741758101263 | 0.23669628490332037 | 1.2046241646800737 | 1.6396919856038026 | 1.7739199789282882 | 1.0 | 0.94931497265256559 |
| M1 \(N=45\) | 12 | 0.2998410355482452 | 0.23502815152638434 | 1.1864058238291151 | 1.6349136797379165 | 1.731476807303588 | 1.0 | 0.92307692307692313 |
| M1 \(N=105\) | 12 | 0.30866975620815057 | 0.23836441828025645 | 1.223885454113065 | 1.6347465097734948 | 1.7739199789282882 | 1.0 | 0.9370629370629372 |
| P3 all | 24 | 0.3601285659000717 | 0.29504890872259426 | 1.243448773777605 | 1.9217733481305757 | 2.1182135331549041 | 0.95833333333333337 | 0.95739130434782604 |
| P3 \(N=45\) | 12 | 0.33540151142667074 | 0.27237639966812593 | 1.2040499398257227 | 1.8495175063754974 | 1.9365364006767223 | 1.0 | 0.965034965034965 |
| P3 \(N=105\) | 12 | 0.3832636091350044 | 0.31772141777706248 | 1.3109937094506292 | 1.9073497456343025 | 2.1182135331549041 | 0.91666666666666663 | 0.97202797202797209 |

For M1 at 1%, 5%, and 10%, the blind predicted-safe counts remained 6, 12,
and 18; observed-safe counts were 12, 18, and 24. False-safe counts were zero
at all thresholds, false-unsafe counts were six, and the worst observed errors
inside the predicted-safe sets were 0.0023131458140526213,
0.0095783633650698419, and 0.023309356560246079. The per-size predicted counts
were 3, 6, and 9 versus observed counts 6, 9, and 12, also with zero false
safe. P3 is diagnostic and did not intervene in the M1 decision.

### Local coupling and matched-size trend

Uniform dilation preserves each family's normalized local-coupling structure.
For \(N=45\),
\(\overline\Lambda/\Lambda_{\max}\) is 0.9705928566009828,
0.7574142282928984, and 0.5993644550204902 for linear, compact, and irregular;
the corresponding fractions above \(0.9\Lambda_{\max}\) are
0.9555555555555556, 0.4, and 0.022222222222222223. For \(N=105\), the ratios
are 0.9871539952996681, 0.8057047053727496, and 0.6524956283739974, with
fractions 0.9809523809523809, 0.4857142857142857, and
0.05714285714285714.

The 12 ratios \(R_{105/45}\), ordered by linear, compact, irregular and levels
1--4, are:

| family | L1 | L2 | L3 | L4 |
|---|---:|---:|---:|---:|
| linear | 0.9990943696202368 | 0.9993171534249062 | 0.9994060093099823 | 0.9993202907893857 |
| compact | 0.8591892364233161 | 0.8987970339956127 | 0.9249101055695469 | 0.9760737958144303 |
| irregular | 0.9284137282405170 | 0.9640561987628187 | 0.9692742569165456 | 1.0001866386186442 |

Their median is 0.9726740263654879, linear-method 90th percentile
0.9993974374579226, minimum 0.8591892364233161, and maximum
1.0001866386186442. Only one ratio exceeds 1; none exceeds 1.10 or 1.25. The
frozen descriptive classification is `NO_SYSTEMATIC_DETERIORATION`. The
combined 48-row table retains the \(N=15\to28\to45\to105\) sequence without
refitting M1.

### Cost, gate, audit, and determinism

The sequential campaign used one worker and one BLAS thread. It accumulated
23209.431251620874 case-seconds over 160 orders. The slowest case took
7431.875041276799 s and the slowest order 3384.133216622984 s. The maximum
final balanced dimension was 6930, peak process memory 7393840 KiB, and the
largest conservative memory estimate 9220780800 bytes. No resource precheck
failed.

Every sufficiency and scientific gate item passed. The exact result is:

```text
PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX
GO_T15_SYNTHESIS_AND_MANUSCRIPT
NO_SYSTEMATIC_DETERIORATION
```

The independent post-revelation audit recalculated and passed these six final
orders: `t14_1_n45_linear_level1`, `t14_1_n45_compact_level3`,
`t14_1_n45_irregular_level4`, `t14_1_n105_linear_level2`,
`t14_1_n105_compact_level4`, and `t14_1_n105_irregular_level1`. Two
analysis-only executions were byte-identical and did not call Model E. The
figure was rendered with the installed TeX backend because the analyzer source
was frozen before revelation; visual inspection confirmed six populated,
legible panels, distinct \(N=45\)/\(N=105\) encodings, visible reference
guides, no clipped axes, and no NaN or infinity. The final suite reported
464 passed in 89.03 s with warnings treated as errors.

The verification environment was Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0,
Matplotlib 3.11.1, scipy-openblas 0.3.33.112.0, an AMD Ryzen 9 9950X, and
186 GiB visible RAM. Commands included:

```bash
.venv/bin/python scripts/preregister_t14_1_large_n.py
.venv/bin/python scripts/run_t14_1_large_n.py --resume --workers 1 --blas-threads 1
.venv/bin/python scripts/run_t14_1_large_n.py --audit-existing --blas-threads 1
MPLCONFIGDIR=/tmp/t14_1_mpl .venv/bin/python scripts/run_t14_1_large_n.py --analyze-only
.venv/bin/python -m pytest -q -W error
git diff --check
sha256sum results/data/t14_1_*.csv results/figures/t14_1_large_n_validation.png
```

### Official T14.1 artifact hashes

The five blind phase-A hashes are:

```text
5b5bb22e92d26f2b74177f97e8b1b9f59857074f8b9309bf51edd59cf1c2daa5  results/data/t14_1_large_n_manifest.csv
8f3d93e2f84afaa96a487ac25ce6db25de57a16e47d925de518e5e17f4df0410  results/data/t14_1_local_coupling.csv
0708f1d31b34ae8d210c9f686e08e71cc1085d00027ebd7d92c42f90443d7041  results/data/t14_1_frozen_predictions.csv
b0586ec363a919586661a4f9f28e4427cf9d36e99c4f4beb54b1d07d6f7e1f72  results/data/t14_1_frozen_protocol.csv
b07a733f83ea5e3ca88bf4527a7481fd2996113a916a0e3ff37c974523ebd440  results/data/t14_1_prior_artifact_hashes.csv
```

The phase-B hashes are:

```text
da476135a66167f1907c99318f7132510d18de9fdaf208fc99e7080d86775f75  results/data/t14_1_model_e_convergence.csv
a03bc3ca119184a4a08c0fbd753a2c4074e8a159a90610e8df2fb1762f863df4  results/data/t14_1_forces.csv
655db6b204b4bb56733f721502f0ffcf32744c031f80f9ee059752fc1b7b3048  results/data/t14_1_case_summary.csv
57a266fafdfaee8b9b902b9b5a7ff163a1cf8ff568d532201d51e239e38422f0  results/data/t14_1_large_n_predictions.csv
f4bb337201b55d4f52d2613490dda372cfe5968299d254a5d050435735d53498  results/data/t14_1_metrics.csv
d4e168992f656ff2ad6cb1fab9af31febeac5f9d4aa75791f4c8f94fa45aa6be  results/data/t14_1_threshold_audit.csv
b52276cae020284c2a0dc5c04aea4be41dcb575af06b2ba8f33643a7fba73ce4  results/data/t14_1_matched_large_n_pairs.csv
424845ed783cb43f67a79def3e4a194056daee71207e2e06d5e32e18f8c38816  results/data/t14_1_combined_scale_sequence.csv
838c55e3a78431313dfcf5a4f5e2daa2d67f2c3b34ca902eaa9b85e0a8102752  results/data/t14_1_performance.csv
88f9b2730f92f938a2f58818fb3e9577a90dbb4ea1ce03a62d2a3672e278e55c  results/data/t14_1_gate.csv
5ae33fc8327e44f80edb7f3d65ade7652ded891bb9f2c4526e3756ee08923dbc  results/figures/t14_1_large_n_validation.png
```

All 95 artifacts recorded from T01--T14 retained their exact path, size, and
SHA-256. The five blind files also remained byte-identical after revelation.

### Scientific limitations

The confirmation covers only two deterministic large sizes, the three
prescribed planar template families, four coupling levels, \(ka=0.1\),
\(f_0=0\), \(f_1=0.8\), positive contrast, identical fixed fluid spheres,
an ideal unbounded fluid, and the approved complete Model-E interaction
force. It does not cover arbitrary geometries, negative contrast, other
frequencies, viscosity, streaming, walls, or dynamics; it is not a theorem or
a universal error guarantee. M1 and P3 were not recalibrated. T15 was not
started.

## T14/P0 — methodological freeze and paper pipeline

P0 started from commit
`e98080da520ec4f1f36b41ece2b76bc281df3a92`. The checkout already contained
completed T14/T14.1 evidence and a mixed user worktree with relocated task
records and two untracked PDFs. Those pre-existing changes remain outside P0.

P0 performed no acoustic solve, fit, recalibration or scientific campaign.
It created the four `docs/PAPER_*.md` products, `campaigns/` schemas and
templates, `src/acoustic_ms/paper_pipeline.py`,
`src/acoustic_ms/plot_style.py` and their two test modules.

The manifest validator checks structural rules plus \(ka=k\,radius\),
multipole-order ordering, unique case IDs, planned-only `TBD` hashes and
panel-source references
(`src/acoustic_ms/paper_pipeline.py::validate_manifest`). The graphic layer
uses STIX Two/STIX with a safe serif fallback, 89/183 mm widths, 9 pt axes,
8 pt ticks/legend, accessible color/marker redundancy and reversible
PDF/SVG/PNG export
(`src/acoustic_ms/plot_style.py::diagnostic_style`,
`save_diagnostic_figure`).

The baseline suite passed 464 tests in 87.93 s. The focused P0 suite passed 14
tests, including valid examples, deliberately invalid mutations, reversible
rc state, noninteractive exports and byte-identical repeated PDF/SVG/PNG
writes. The final full suite passed **478 tests in 88.04 s** with warnings
treated as errors. The explicit 95-artifact preregistration audit also passed,
and `results/` plus `papers/` had no diff.

There is no complete-dimer \(B_E\) aggregator, Model-E subset
inclusion–exclusion or \(\Phi_E^{(3..5)}\) API. Existing connected expansions
apply only to C or D. P1 is therefore the canonical dimer benchmark and P2
remains unopened.

Existing timing tables provide no dimer unit cost. Measured references are
18.30 s/order for T14 and 145.06 s/order for T14.1; they are not dimer bounds.
P1 must start with a timed preregistered pilot. The P0 recommendation is
`GO_P1_WITH_CONDITIONS`.

## P1.1 — canonical-dimer decisions frozen, cases disabled

P0 was finalized through GitHub PR #1. Its head remained
`893577fb9f2745ee6cd1d0d6deea9cb9276a6fa4`, the diff contained exactly 17
files and no path under `results/` or `papers/`, and GitHub created merge commit
`926e639fe2d327eacd09a2542208500891399687` on `main`.

The scientific review of draft PR #2 was applied on
`agent/p1-1-decision-record` while the 30 pre-existing local changes remained
unstaged and outside scope. `docs/P1_1_DIMER_DECISION_RECORD.md` now records
the approved physical grid, exact ordering, convergence/failure policy,
provisional resources, separate pilot and reproducibility audits. The
remaining gate is procedural: P1.2 and P1.3 must implement and audit future
`B_E` work, while P1.4 must generate final manifest hashes and authorize any
enablement.

Schema `1.1.0` is additive and keeps campaign constants at the manifest level
while requiring `ka`, `k_rad_m`, `material_id`, `material_model`, `f0`,
`f0_applicable`, `f1`, `distance_ratio` and `theta_rad` per case. The validator
routes versions explicitly; schema `1.0.0` and its example remain valid without
reinterpretation. Semantic checks cover per-case `ka=k a`, contiguous order,
planned-case disablement, rigid-sentinel handling and rotational-audit twins.

`campaigns/p1/campaign_manifest.yaml` freezes 102 unique ordered IDs: 96
primary cases from two `ka` values, six materials and eight separations at
zero angle, followed by six `pi/4` rotational audits linked to zero-angle
twins. Every case is disabled and the final hash remains `TBD` for P1.4. The
separate `campaigns/p1/pilot_manifest.yaml` contains one disabled
`development` rigid case at `ka=0.1`, `d/a=2.1`, `theta=0`; it is explicitly
excluded from future P1.6 scientific tables.

The frozen numerical policy evaluates orders 2--21, forbids stopping before
5 and requires two consecutive convergence passes in all applicable channels.
Resources remain provisional at one worker, one BLAS thread, 4 GiB/case,
30 min/case and 12 h total. Two no-solver regenerations and the six rotational
audits are assigned to later authorized execution. The common-order `B_E`
audit remains assigned to P1.3.

This revision implemented no `B_E`, called no solver, executed neither pilot
nor campaign, generated no force and changed no path under `results/` or
`papers/`. P1.2--P1.6 remain unopened. The focused manifest suite passed
**26 tests in 0.52 s**; the complete suite passed **493 tests in 87.08 s** with
warnings treated as errors.


## P1.2 — independently converged complete-dimer API

P1.1 was finalized from PR #2 after confirming its exact head
`24ec933f366cb4950ad4050d83ce804d89d4eb43`. The PR was marked ready and
merged into `main` as merge commit
`4a5b58408dc40302568758b2bdea54701beb4747`. P1.2 was then developed from
that updated `main` on `agent/p1-2-be-api` in an isolated clean worktree, so
the pre-existing changes in the original checkout remained untouched.

The new public entry point is:

```python
solve_model_be_nodal(
    positions_xyz, k, radius, energy_density, f0, f1, *,
    lmax_min=2, lmax_max=21, minimum_stop_lmax=5,
    convergence_tolerance=1.0e-5,
    solver=solve_model_e_nodal,
) -> ModelBEResult
```

For every lexicographically ordered pair \(i<j\), the API sends the original
two positions to Model E in that orientation and independently evaluates
orders until all applicable force channels have two-step convergence, never
before order 5. The final interaction-force rows are associated directly with
particles \(i\) and \(j\); no action--reaction assumption or frame rotation is
introduced.

The P1.2 follow-up fixes the distinction between historical and current
convergence. `confirmation_lmax` retains the first historical two-step pass,
whereas `confirmed` is true only when the two most recent changes at the
current order are both applicable and within tolerance. The solver stops only
when all applicable channels are simultaneously confirmed in that final
window. A nonmonotonic fake channel that first confirms at order 4, varies at
order 5 and stays fixed thereafter is unconfirmed at the cap 6; with one more
unchanged step it reconfirms and stops at order 7.

`ModelBEResult` returns the global \(B_E\) vector only when all pairs are
eligible. Its deterministic ledger retains each pair's individual dimer
forces, attempted/evaluated/final/failed orders, channel histories,
applicability, convergence, final numerical diagnostics and explicit failure
stage/reason. Pair failures do not suppress later ledger entries, but no
partial or imputed global force is exposed.

The established Model-E gate is centralized in
`evaluate_model_e_numerical_diagnostics` and reused without changing its
thresholds. Position, physical parameter, order, tolerance, rigid sentinel,
non-overlap and planarity checks run before pair evaluation. Historical Model
B in `comparison.py` and \(B_L\) in `transferability.py` were not modified.

The P1.2-specific suite uses only injected fake solvers and passed **33 tests
in 0.57 s**. It covers vector accumulation, original orientation,
lexicographic pair order, independent final orders, null-channel
applicability, deterministic repetition, local and malformed-solver failures,
all-or-nothing global eligibility, input validation, and every established
numerical gate. The complete repository suite passed **526 tests in 97.16 s**
with warnings treated as errors.

No pilot, campaign or production calculation ran. All confirmatory and
development cases remain disabled, hashes remain deferred to P1.4, and no
path under `results/` or `papers/` changed. Full physical dimer identities,
common-order sensitivity, rotation/reflection, action--reaction and
asymptotic-limit tests remain P1.3. The handoff decision is `GO_P1.3`; it does
not authorize P1.4 or campaign execution.
