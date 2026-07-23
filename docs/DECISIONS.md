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
