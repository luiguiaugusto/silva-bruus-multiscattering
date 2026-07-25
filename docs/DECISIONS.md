# Decisions

- Python 3.11+ is the official implementation language.
- Scientific routines live in the importable `src/acoustic_ms` package.
- Notebooks, when introduced, are demonstrations only and will never contain
  the sole version of a scientific routine.
- T01 depends only on NumPy at runtime and pytest for development.
- Matplotlib is an optional `plot` dependency, used only by the reproducible Figure 2 script; the scientific package does not import it.
- The T02 corrected formula is a published two-particle benchmark, not a general multiple-scattering solution.
- T03 uses dense NumPy linear algebra at `Lmax=1` only; SciPy and SymPy are runtime dependencies for special functions and Gaunt coefficients.
- `Lmax=1` truncates multipolar order but not the number of rescattering events; no radiation-force API is added in T03.
- T04 implements Model C at Rayleigh level using Eq. (22)/(27) cross terms only, with no scattered--scattered products.
- The T03 production solver remains at `Lmax=1`; T04 uses local evaluation through `ell=2` only and reports no three-particle force results.
- T05 is restricted to canonical N=3 trimers. Model B and C deliberately share the T04 solver and observable, so C-B isolates multibody rescattering at Lmax=1.
- The scalar nodal-plane oracle is test-only. No zero-total-force constraint is imposed on the global scalene interaction observable, and T05 does not measure multipolar correction or introduce Model D.

- T05.1 defines numerical nullity relative to the global configuration scale, \(128\,\epsilon_{\mathrm{mach}}F_{\mathrm{scale}}\), without an absolute floor. Correction amplitudes are RMS vector magnitudes per particle, not component RMS; only derived metrics and corresponding artifacts changed, while A, B, C and their equations remained unchanged.
- T05.1a, T05.1b, and T05.1c are exclusively documentary. Binary determinism is assessed in the same numerical environment.
- T06 reports only planar \(N=4\) results at Rayleigh \(L_{\max}=1\); Model D and higher scattered multipoles remain out of scope.
- For \(N=4\), \(C-B\) is not an irreducible four-body contribution: it equals the embedded three-body sum plus \(\boldsymbol{\Phi}^{(4)}\).
- Every connected term is built exclusively from subsets solved by the same Model C. Model A is a comparison baseline and does not define \(\boldsymbol{\Phi}^{(3)}\) or \(\boldsymbol{\Phi}^{(4)}\).
- The decomposition is vectorial. No zero-sum constraint is imposed on the approved irregular-quartet observable.

- T06.1 is post-processing of the already versioned T05/T06 CSVs. It performs no new trimer or quartet force sweep; the only additional Model C evaluation is one centered dimer for the \(N=2,3,4\) comparison.
- \(\Lambda_{\max}\) is an exploratory geometric diagnostic. Its grouped improvement is reported descriptively and no universal validity threshold is defined.
- Within each fixed-shape dilation family, \(\Lambda_{\max}=C_g\eta\), so exponent and log-space fit quality are necessarily unchanged from the \(\eta\) fit.
- T06.1 changes no force model, solver, connected-body definition, or protected T03--T06 artifact. T07, Model D, higher multipoles, and new cluster families remain outside its scope.
