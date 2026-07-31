# Silva--Bruus multiscattering

Python tools for studying acoustic interaction forces between identical spheres in a pressure nodal plane. Models A--C provide the Silva--Bruus, matched Rayleigh pairwise, and global \(L_{\max}=1\) hierarchy; T07 adds balanced multipolar Model D. T11 adds Model E: exact lossless-fluid Mie coefficients, global multiple scattering, and the complete multipolar force including external--scattered and scattered--scattered channels. Earlier A--D APIs and results retain their approved meanings.

T12 audits the frozen T08 \(\rho_1\) law with 28 preregistered Model-E
sentinels in the calibration domain \(N\leq4\). It uses the complete
interaction force in three dimensions, does not recalibrate the predictor, and
does not evaluate the \(N=6,10\) holdout. Sentinel approval is not external
validation and is not a universal criterion.

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

## T06.1 scaling diagnosis

T06.1 post-processes the existing 1,920 quartet configurations without
re-running their force solves. The connected three-body amplitude scales
approximately to first order and the irreducible four-body amplitude to second
order in \(\eta=|f_1|(a/d_{\min})^3\) over the audited data. The collective
predictor \(\Lambda_{\max}\) improves the descriptive grouped collapse for the
three fixed geometry families, but it is exploratory and is not a universal
validity criterion. Model D, T07, higher multipoles, new force sweeps, and
\(N>4\) remain outside this result.

## T07: multipolar Model D

T07 adds planar multipolar Model D for identical nodal clusters with \(N\leq4\). It uses a balanced globally coupled system, preserves the external--scattered force observable, and studies total-force and connected-term convergence separately. Model D at \(L_{\max}=1\) is numerically equivalent to Model C. The implemented coefficients are the leading Rayleigh term at each multipole, not an exact finite-frequency T-matrix; viscosity, streaming, walls, dynamics, and scattered--scattered force terms remain out of scope.

## T08: transferability and frozen article data

T08 evaluates 312 deterministic configurations at \(ka=0.1\), spanning
\(N\in\{2,3,4,6,10\}\), three cluster families where applicable, four
positive dipole contrasts, and six separations. The matched pairwise baseline
\(B_L\) sums isolated Model-D dimers at the same multipolar order as the global
solution. Predictor selection and empirical thresholds use only \(N\leq4\);
\(N=6,10\) is held out until final evaluation. Cross-validation selected the
dipolar balanced-operator spectral radius \(\rho_1\), and the prespecified
diagnostic transferability criterion was supported within this sampled nodal
domain. This is not a universal validity theorem. T08 closes the computational
sweeps used for the article; the data are frozen.

## T09: analytical foundation of rho_1

T09 reduces the balanced nodal \(L=1\) operator to an exact \(N\times N\)
dipolar matrix,

\[
(K_b)_{ij}=\frac{f_1}{2}\left(\frac{a}{r_{ij}}\right)^3
e^{ikr_{ij}}(1-ikr_{ij}),\qquad i\ne j.
\]

This derives the inverse-cube scale, connects \(\rho_1\) to convergence of the
Neumann rescattering series, and explains why the leading pairwise correction
and connected three-body term are first order while the connected four-body
term begins at second order. The derivation does not make the empirical T08
fit or thresholds universal. Reproduce the symbolic, spectral, and
Neumann-series audits with:

```bash
.venv/bin/python scripts/analyze_t09_rho_foundation.py
```

The full derivation is in
[`TAREFA_T09_FUNDAMENTACAO_ANALITICA_RHO1.md`](TAREFA_T09_FUNDAMENTACAO_ANALITICA_RHO1.md).

## T10: exact isolated-sphere Mie coefficients

T10 adds the exact diagonal partial-wave response of a homogeneous lossless
fluid sphere and the analytic rigid-sphere limit. A deterministic campaign
compares orders \(\ell=0,\ldots,5\) with the existing leading Rayleigh
coefficients for \(10^{-3}\le ka\le0.1\). The new coefficients are deliberately
not connected to Model D: exact single-sphere Mie coefficients are not a
complete collective-force model. Global integration and the
`scattered--scattered` channel are implemented separately by T11 Model E.

```bash
.venv/bin/python scripts/analyze_t10_mie_rayleigh.py
```

## T11: complete-reference Model E

Model E uses the exact T10 diagonal T-matrix and evaluates the complete
partial-wave radiation force. T11.1 stabilizes its linear algebra by solving
\((I-D^{1/2}UD^{1/2})q=D^{1/2}a\), followed by \(d=D^{1/2}q\) and
\(b=a+Ud\), without forming \(D^{-1/2}\). Its interaction force is
split exactly into external--scattered and scattered--scattered terms. An
independent surface-stress quadrature validates the force normalization, and a
compact six-case campaign audits \(L_{\max}=2,\ldots,9\) without opening the
T12 sentinel campaign.

```bash
.venv/bin/python scripts/analyze_t11_model_e.py
.venv/bin/python scripts/analyze_t11_1_model_e_stability.py
```

The current scope remains lossless identical fluid spheres in an ideal fluid,
fixed in the nodal plane. Convergence internal to Model E does not validate the
empirical \(\rho_1\) thresholds; that transfer study belongs to T12--T14.

## T12--T12.1: sentinel audit and convergence diagnosis

T12 tested the frozen T08 \(\rho_1\) law against 28 preregistered Model-E
sentinels with \(N\leq4\). T12.1 extends exactly the ten incompletely resolved
cases through at most \(L_{\max}=21\), without changing Models A--E or any
earlier artifact. All 28 interaction-force references are now directly
confirmed; two scattered--scattered channels remain explicitly
`unconfirmed_at_21`.

Leave-\((N,\mathrm{family})\)-out diagnostics show that a recalibrated
\(\rho_1\) power law meets the preregistered descriptive targets on this small
sentinel set. This supports only a separate T12.2 recalibration study. It does
not validate a universal criterion, does not open the \(N=6,10\) holdout, and
does not authorize T13 or T14.

## T12.2: controlled rho1 recalibration

T12.2 performs the preregistered seven-fold leave-\((N,\mathrm{family})\)-out
recalibration of the same power law against the 28 confirmed Model-E
interaction-force errors. The candidate improves the frozen T08 calibration,
but retains one false-safe case at the 10% threshold. The resulting decision is
`NO_GO_T13_RHO1_NOT_QUANTITATIVE`; no alternative predictor, new force solve,
or \(N=6,10\) holdout case was evaluated.

## T12.3: grouped mechanistic validity criterion

T12.3 post-processes the same 28 confirmed \(N\leq4\) sentinels with nested
leave-\((N,\mathrm{family})\)-out validation. The primary candidate uses the
geometric coupling sum \(\Lambda_{\max}\); a pre-specified bivariate
\((\Lambda_{\max},\rho_1)\) model is diagnostic. Training-only residual
margins remove false-safe classifications while retaining the preregistered
minimum safe coverage, so the literal result is
`GO_T13_VALIDATE_LAMBDA_MAX`. This means only that \(\Lambda_{\max}\) is a
candidate for later external validation. T12.3 does not read the \(N=6,10\)
holdout, run a new force solve, start T13/T14, or establish a universal rule.

## T13: external validation of the frozen Lambda-max criterion

T13 preregistered and published a response-blind sample of 24 T08 holdout
clusters before any new Model-E solve: four target levels in each combination
of \(N\in\{6,10\}\) and the linear, compact, and irregular families. The
frozen M1 law in \(\Lambda_{\max}\) passed its literal external gate with all
24 references eligible, zero conservative false-safe cases at 1%, 5%, and
10%, global log-RMSE 0.465095692587546, factor-two fraction 0.875, and
Spearman 0.964347826086956.

This supports transfer of the frozen criterion only to the sampled planar
domain at \(ka=0.1\). It is not a universal theorem, does not refit M1, and
does not start T14. Reproduce the already frozen campaign analysis without
calling Model E using:

```bash
.venv/bin/python scripts/run_t13_external_validation.py --analyze-only
```
