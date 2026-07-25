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
