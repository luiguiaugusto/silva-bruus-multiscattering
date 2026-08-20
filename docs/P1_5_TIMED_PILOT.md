# P1.5 timed resource pilot

Status: **single real attempt complete; `GO_P1.6A_BLIND_FREEZE`**.

This protocol authorizes exactly one execution of
`p1_pilot_rigid_ka010_d0210_t000`. It does not authorize P1.6 or any
confirmatory case. The two P1 manifests remain byte-identical to P1.4:

- confirmatory SHA-256:
  `9d360de6e61d901cff3f84c477f367773251103db12386dbb8156bd1ec2addca`;
- pilot SHA-256:
  `d8f56ce20f6f0821d84fd6f36e1f76c855f63f55d809ba9a7201ba52097a43bf`.

## Frozen execution

The runner validates the public locks before execution, requires the 102
confirmatory cases to remain disabled and accepts only the one enabled pilot.
The physical input is rigid, `ka=0.1`, `k=0.1 m^-1`, `a=1 m`, API sentinel
`f0=0` with `f0_applicable=false`, `f1=1`, `d/a=2.1` and `theta=0` at centered
positions `(-1.05,0,0)` and `(1.05,0,0)` metres.

`solve_model_be_nodal` receives the real `solve_model_e_nodal` and evaluates
`L=2,...,21`, with minimum stop 5, tolerance `1e-5`, the two-change final
window and all four applicable channels. There is one worker and one BLAS
thread. An internal `SIGALRM` enforces 1800 s; `RLIMIT_AS` conservatively
enforces 4 GiB while peak process RSS is recorded independently. No retry,
resume, retuning or overwrite is permitted.

The exact executed shell invocation was:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONHASHSEED=0 \
  timeout --signal=TERM --kill-after=30s 1860s \
  /home/luigui/Documents/silva-bruus-multiscattering/.venv/bin/python \
  scripts/run_p1_5_timed_pilot.py --execute
```

The 1860 s outer watchdog is infrastructure-only and allows 60 s to serialize
the internally enforced 1800 s controlled timeout.

## Timing definitions

For every attempted order, a scoped single-worker instrumentation around the
unchanged Model E records:

- `assembly_diagnostics_seconds`: Model E wall time outside the linear solve
  and complete-force calls, including modal/translation assembly and numerical
  diagnostics;
- `linear_solve_seconds`: wall time inside the existing
  `numpy.linalg.solve` call;
- `force_postprocess_seconds`: wall time inside the existing complete-force
  evaluations;
- `order_wall_seconds`, accumulated time, total case time and process peak RSS.

The instrumentation restores the original callables after every order and
does not change Models B, `B_L`, E or `B_E`.

## Atomic single-use artifacts

The four frozen paths are used unchanged:

- `campaigns/p1/pilot/data_raw.csv`;
- `campaigns/p1/pilot/data_derived.csv`;
- `campaigns/p1/pilot/data_plot.csv`;
- `campaigns/p1/pilot/failures.csv`.

`campaigns/p1/pilot/performance.csv` is the additional performance table.
All five files are written in a hidden sibling directory, flushed, and then
published by one atomic directory rename. If the final directory exists, the
runner refuses before calling any solver. A failed publication removes the
hidden incomplete directory.

Raw data retain every attempted order, particle and the `total`,
`interaction`, `external_scattered` and `scattered_scattered` channels,
together with diagnostics and final eligibility. Every artifact is
`classification=development`, and every force row has
`include_in_scientific_tables=false`. Derived and plot tables contain only
resource metrics, never force metrics.

After the single execution, `--verify-derived` generates derived/plot bytes
twice from raw/performance/failure bytes, compares both generations and then
compares them with the published files. It calls no solver.

## Frozen classification

- `GO_P1.6A_BLIND_FREEZE`: execution and serialization are complete and both
  resource limits pass. `unconfirmed_at_21` is explicitly accepted here with
  `eligible=false` and no scientific-force promotion.
- `NO_GO_P1.6_RESOURCE_LIMIT`: the 1800 s or 4 GiB limit is exceeded.
- `INCONCLUSIVE_P1.5`: infrastructure failure or incomplete artifacts.

The response-blind runner was committed and pushed before execution as
`a5a2a9c58f5e65b7986e24c7c64879246d946131`.

## Observed resource result

The only real invocation ran from `2026-08-19T20:17:27.345044Z` through
`2026-08-19T20:25:41.478313Z`. It evaluated every order 2--21 once, with no
retry. The final result was `converged=false`, `eligible=false`,
`failure_stage=convergence` and `stop_reason=unconfirmed_at_21`.

At `L=21`, total, interaction and external--scattered were confirmed in the
final two-change window. Scattered--scattered was applicable and had final
successive change `7.4980301857038341e-06`, but lacked two consecutive final
passes, so the all-channel stop correctly remained false. All Model E
numerical diagnostics passed.

Total case wall time was `494.13323493499774 s`. Summed per-order Model E wall
time was `494.1200086578028 s`, decomposed as
`494.0773162767291 s` assembly/diagnostics, `0.01332716818433255 s` linear
solve and `0.029365212889388204 s` force postprocessing. Peak RSS was
`311857152 bytes`. The final full/active mode counts were 484/231 per particle
and the final system dimension was 462. Both frozen resource limits passed.

Artifact SHA-256 values are:

- `data_raw.csv`:
  `a4416cae58654371ddcf680ce1a8470ab227c58760b8e1d507893e91883574da`;
- `data_derived.csv`:
  `ccd1f7a1aac92a25c51dfb822530fc55f2c27d77c0d85bbcc4215397a3bf2026`;
- `data_plot.csv`:
  `9a94fb1203ae122f89a3eb3f49074bea89f1fd7879224b58c6af7e9cafc38424`;
- `failures.csv`:
  `cd60af766e7340aa04e1b3a1fb2f4b7948f7901163ea75fb5ac42ef4e93e3e8f`;
- `performance.csv`:
  `9bb573c524a31a183856289610dc91478d14c8cacb254c9f3a85a2ad00048222`.

Two no-solver derivations were byte-identical to the published derived/plot
files. Focused tests passed `34 passed in 0.67s`; the complete suite with
warnings as errors passed `574 passed, 1 xfailed in 652.79s`. The sole xfail
is the unchanged P1.3 G7. No file under `papers/` or `results/` changed.

The frozen classification is `GO_P1.6A_BLIND_FREEZE`. This closes P1.5 and
does not authorize or execute P1.6.
