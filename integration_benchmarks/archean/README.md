# Archean integration benchmark

These are opt-in numerical experiments, separate from CI smoke tests. The first
case compares default CVODE, order-one CVODE, and experimental linearly implicit
Euler PTC. All methods use the same atmosphere and physical configuration.

## Inputs

- `mechanism1.yaml` and `settings1.yaml` define the case.
- `state1.pkl` is the existing saved atmospheric state used by the restart case.
- `benchmark1.py` constructs the climate-based fresh initial atmosphere. Its
  current direct invocation selects order-one CVODE; the comparison driver
  explicitly chooses each method independently of that setting.

The runner never overwrites `state1.pkl`. It copies inputs and exports an HDF5
initial state under `archean/results/restart/` or `archean/results/fresh/`. Input SHA-256 hashes
are recorded when preparing a case. Add the original inputs to Git along with
these scripts when preserving the benchmark. Generated outputs are ignored.

Activate the Python/compiler environment used for the CMake build, then run from
the repository root:

```sh
cmake --build build --target photochem_test_api photochem_test_ptc -j 4
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case archean --modes cvode cvode1 ptc --steps 7000 --seconds 120
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case archean --fresh --modes cvode cvode1 ptc --steps 7000 --seconds 60
PYTHONPATH=. python integration_benchmarks/compare_ptc.py --case archean
```

The driver is compiled using CMake's existing compiler and link configuration.
Use `--build PATH` for a different build tree. Use `--reuse-inputs` to reuse an
already exported initial state instead of regenerating it. Each mode runs in a
fresh process with one OpenMP thread. `--seconds` limits integration CPU time;
there is also a more generous process wall-time timeout. Compilation and input
preparation are excluded from integration timings.

Additional local experiments:

```sh
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case archean --reuse-inputs --modes ptc-residual --steps 1000 --seconds 45
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case archean --reuse-inputs --modes ptc-periodic ptc-toa --steps 1000 --seconds 35
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case archean --fresh --reuse-inputs --modes ptc-floor12 --steps 3000 --seconds 60
```

`ptc-periodic` switches to periodic profile maintenance and deliberately changes
the prescribed temperature by 10% after setup. `ptc-toa` changes the target TOA
pressure and tightens its maintenance thresholds. These exercise real shared
resynchronization transactions; their final atmospheres are not reference-state
comparisons. `ptc-floor12` changes the step-control abundance floor from
`1e-20` to `1e-12` for a sensitivity experiment, not as a recommended default.
`cvode1-tight` is an optional order-one tolerance sensitivity run
(`rtol=1e-6`, `atol=1e-23`).

## What is measured

All runs use `conv_longdy=0` and `equilibrium_time=1e17` s. Only accepted steps
advance PTC pseudo-time. Its clock is not an accurate physical trajectory.

Outputs include per-step CSV traces, an explicit stopping reason, CPU time,
accepted steps, rejection counts, callback counts, maintenance updates, final
HDF5 snapshots, and final tendency/budget diagnostics. PTC reports cumulative
linear solves and callbacks. CVODE callback and error-test counts currently
refer only to its final solver segment because robust restarts reset SUNDIALS
counters; they are not directly comparable to PTC totals. Unavailable CVODE
linear-solve counts are marked `-1`.

PTC uses the existing approximate Jacobian, rebuilt at every accepted state.
Both timestep policies use the same 25% scaled-change safeguard. The default
uses ATMOS-style growth factors; the alternative uses the bounded residual ratio
without an extra growth multiplier. Residual ratios use weights frozen at the
accepted state. The adapter rejects time-dependent photon and rate callbacks.

The atmospheric adapter uses a `1e-20` abundance floor for correction scaling.
It rejects significant negative trials. Negative values below `1e-4` of the
correction scale may be projected to zero before reevaluating the trial RHS;
projection counts and maximum scaled corrections are recorded. With the default
floor and nonnegative starting state, such a crossing is limited to trace
abundances below approximately `1e-24`. Initial trace negatives are checked
against that same floor and projected before starting. This is an explicit
numerical allowance, not a particle or species-specific repair.

The prototype supports the time-threshold stopping rule only. General
history-based convergence, public backend selection, Python bindings, and
accurate time-dependent evolution remain deferred. Do not interleave a PTC
session with another integrator or external model edits; reinitialize after an
external model change. PTC reinitializes its own storage after shared maintenance
changes the atmosphere and retains cumulative pseudo-time and cost counters.

`budget_reconciliation_error` compares each column-integrated full RHS with the
process-resolved production/loss diagnostic, normalized by the sum of absolute integrated production and
loss contributions. This checks diagnostic consistency, not independent conservation proof.
Elemental column tendencies are also reported; these describe an open atmosphere
with boundary exchanges and must not be interpreted as closed-system errors.

`compare_ptc.py` compares same-grid final profiles and significant species
columns against order-one CVODE, writes `comparison.json`, and produces
`archean/results/progress.pdf`. It skips direct comparisons when grids differ. Bounded,
unconverged runs and deliberately modified physical cases are not accuracy
references even when a numerical difference can be calculated.
