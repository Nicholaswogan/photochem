# Integration benchmarks

Local, opt-in experiments for stiff steady-state solvers. Each case owns its
physical inputs, initialization, and convergence target. Only the compiler,
Fortran comparison driver, and analysis helpers are shared. Generated results
are ignored by Git; long integrations are not part of CI.

```text
integration_benchmarks/
    common/       shared runner, Fortran driver, and comparison helpers
    archean/      original climate-initialized and saved-state case
    venus/        ported legacy Venus case, with fresh and saved-atmosphere starts
    run_ptc.py    entry point (--case archean or --case venus)
    compare_ptc.py
```

See [Archean](archean/README.md) and [Venus](venus/README.md) for physical controls,
inputs, and migration details. Results are discussed in each case's
`PTC_RESULTS.md`. The old root-level `benchmark1.py`, mechanisms, and checkpoint
now live under `archean/`.

Activate the Python/compiler environment used for the CMake build, then run from
the repository root:

```sh
cmake --build build --target photochem_test_api photochem_test_ptc -j 4
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --modes cvode cvode1 ptc --seconds 45
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --fresh --modes cvode cvode1 ptc --seconds 45
PYTHONPATH=. python integration_benchmarks/compare_ptc.py --case venus
```

Omitting `--case` retains Archean as the default. `--fresh` selects a newly
constructed atmosphere; otherwise the case's existing saved atmosphere is used.
`--reuse-inputs` reuses an already exported initial state and controls.
`--jacobian-method analytical` explicitly selects the analytical Jacobian for
every method, including when reusing older controls. Autodiff and finite
differences can also be selected for diagnostic comparisons. `--build`
selects a different CMake build directory. `--steps` bounds accepted steps
(default 7000); `--seconds` bounds integration CPU time (default 120). There is
also an outer process wall-time timeout. A step can exceed the CPU budget before
returning; the limit is checked between accepted steps.

Case preparation copies inputs and writes `initial.h5`, `controls.json`, and a
Fortran namelist under `<case>/results/{restart,fresh}/`. Model-state snapshots
omit several physical controls, so the driver reapplies that namelist after
restoration for every method. Input SHA-256 hashes are recorded. No preparation
or comparison writes to the original input fixtures. All runs use one OpenMP
thread; solver initialization is timed, while model construction, compilation,
and final diagnostics are excluded from integration CPU time.

Modes are `cvode`, `cvode1`, `ptc`, `ptc-residual`, and experimental sensitivity
variants. `ptc-floor12` uses a `1e-12` correction-scale abundance floor instead of
`1e-20`; it is not a recommended default. `cvode1-tight` tests `rtol=1e-6`,
`atol=1e-23`. The Archean-only `ptc-periodic` and `ptc-toa` diagnostics deliberately
perturb its persistent profile or TOA target. They are not same-physics accuracy
comparisons and should not be applied to Venus, which has fixed altitude profiles.

The PTC policy is unchanged across cases. Both variants retain the 25% scaled
correction safeguard. The residual variant uses a bounded ratio without an
additional growth multiplier. PTC is experimental and Fortran-only; production
CVODE defaults and the Python API are unchanged. Only the time-threshold
criterion is tested (`conv_longdy=0`). PTC time is pseudo-time, not an accurate
physical trajectory.

Each mode writes a CSV trace, summary with stopping reason, final HDF5 state,
and per-species column-budget CSV. PTC counters are cumulative. CVODE RHS,
Jacobian, and error-test counters refer only to the final solver segment because
robust restarts reset SUNDIALS statistics; unavailable linear-solve counts are
`-1`. Normalized budget reconciliation uses the sum of absolute integrated
process contributions, so signed trace contributions do not cancel its scale.
This is an RHS/diagnostic consistency check, not independent conservation proof.

`compare_ptc.py` writes snapshot differences and a progress PDF below the case's
`results/` directory. It labels whether each run reached its target and skips
direct comparison of different grids. Differences involving unconverged runs
are descriptive only, not accuracy estimates.
