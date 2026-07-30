# Silva--Bruus multiscattering

Python tools for studying acoustic interaction forces between identical spheres in a pressure nodal plane. Models A--C provide the Silva--Bruus, matched Rayleigh pairwise, and global \(L_{\max}=1\) hierarchy; T07 adds balanced multipolar Model D. T08 tests transferability through \(N=10\) using deterministic planar families, while preserving \(N\leq4\) as calibration and \(N=6,10\) as an untouched external holdout. The project remains restricted to the external--scattered force observable and leading Rayleigh coefficients at each retained multipole.

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
