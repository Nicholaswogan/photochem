# Venus comparisons with the analytical Jacobian

**Every result in this report uses `jacobian_method=3`**, for CVODE and PTC.
The current Venus case defaults to the analytical Jacobian. Earlier autodiff
runs, including the interrupted longer run, are archived separately under
ignored `results/legacy-autodiff/` and are not mixed into these tables.

## Setup

The nominal legacy Venus chemistry setup has 114 evolved species on 100 layers
(11,400 equations). T and Kzz are prescribed on altitude. There is no climate or
haze-microphysics iteration. Both starts retain the original physical controls,
cloud radii/optics, UV absorber, and time target `3.15e13 s`.

- **Saved atmosphere:** the legacy nominal text output, restored with current
  thermodynamic and initialization rules. Each method gets up to 45 integration
  CPU seconds.
- **Fresh start:** VIRA T/P, the supplied Kzz profile, and the prescribed initial
  gas composition. Each method gets **300 integration CPU seconds**.
- Both batches have a 100,000-step safety limit so the CPU budget, rather than
  the earlier 7000-step limit, can determine the result.
- All methods in each batch use byte-identical initial HDF5 inputs, the same
  physical settings, upwind molecular diffusion, `rtol=1e-3`, and `atol=1e-17`.
  History-based early convergence is disabled (`conv_longdy=0`).
- Default PTC retains the unchanged Archean policy: a 25% scaled-change limit
  and a `1e-20` correction-scale abundance floor. The `floor12` variant changes
  only that floor to `1e-12`. The residual variant remains the pure bounded
  residual-ratio update, without an extra growth multiplier.

Each method runs in an independent process with one OpenMP thread. Runs in a
batch were launched concurrently. Times below are per-process integration CPU
times; compilation, model construction, and final diagnostics are excluded,
while solver initialization is included. Wall times are not comparable to CPU
times under concurrent execution. These are individual numerical experiments,
not precision speed benchmarks.

The original mechanism, reproducible continuity repair, and UV-opacity migration
are documented in [README](README.md). These necessary port differences remain
in effect for every method.

## Saved-atmosphere results: 45-second CPU budget

| Method | Reached target? | Steps | CPU seconds | Final time / pseudo-time (s) |
|---|---|---:|---:|---:|
| Default CVODE | Yes | 119 | 10.99 | 3.3903e13 |
| Order-one CVODE | Yes | 93 | 9.50 | 3.2195e13 |
| Default change-based PTC | No, CPU limit | 252 | 45.10 | 6.4116e4 |
| Pure residual-ratio PTC | No, CPU limit | 292 | 45.04 | 1.5920 |
| PTC, abundance floor `1e-12` | Yes | 41 | 6.86 | 3.5248e13 |

Default PTC is repeatedly limited by upper-atmosphere trace chemistry, including
N2D in the top layer. The larger floor allows this warm start to converge much
more readily. That sensitivity is a warning about the current scaling policy,
not evidence for selecting a different floor separately for every atmosphere.

Among the runs reaching the target:

| Comparison against order-one CVODE | Max significant column difference | Max profile difference above `1e-6` gas mixing ratio |
|---|---:|---:|
| Default CVODE | 0.063% | 0.028% |
| PTC with `1e-12` floor | 0.041% | 0.020% |

Significant columns exceed a `1e-10` column fraction. Smaller trace species have
larger differences: the larger-floor PTC maximum profile difference is about
0.69% above a `1e-10` mixing ratio and 8.8% above `1e-20`. Its maximum absolute
mixing-ratio tendency is `2.81e-13 /s`, compared with `3.00e-14 /s` for order-one
CVODE. Reaching the time threshold does not establish identical residual
accuracy or an independently exact root.

## Fresh-start results: five CPU minutes per method

**None reached the `3.15e13 s` target.** Every method used at least 300 CPU
seconds and stopped at the CPU limit. Neither the step cap nor a solver failure
ended these runs.

| Method | Steps | CPU seconds | Final time / pseudo-time (s) | Rejected attempts / final-segment CVODE error-test failures |
|---|---:|---:|---:|---:|
| Default CVODE | 7706 | 300.28 | 3.0274e8 | 3 |
| Order-one CVODE | 7276 | 300.49 | 8.2239e8 | 0 |
| Default change-based PTC | 2113 | 300.15 | 3.3629e3 | 388 |
| Pure residual-ratio PTC | 2216 | 300.05 | 1.6599e-10 | 46 |
| PTC, abundance floor `1e-12` | 2182 | 300.08 | 6.6488e6 | 136 |

The rejection columns have different meanings: PTC counts all rejected trial
steps, while CVODE's counter covers only its final robust-stepper segment. They
must not be treated as equal measures of solver effort.

The final maximum absolute mixing-ratio tendencies were:

| Method | Maximum absolute mixing-ratio tendency (`s^-1`) |
|---|---:|
| Default CVODE | 4.73e-11 |
| Order-one CVODE | 1.09e-12 |
| Default change-based PTC | 7.90e-7 |
| Pure residual-ratio PTC | 1.88e-6 |
| PTC, abundance floor `1e-12` | 2.08e-9 |

Order-one CVODE advanced farthest in time and had the smallest reported tendency
in this batch, but it still did not meet the requested stopping rule. The
snapshots are not converged reference solutions. Differences between them are
not estimates of steady-state accuracy.

Default PTC initially spends many steps controlling the growth of newly
produced O, then encounters aerosol and other trace-species restrictions. Its
last reported limiter was S8aer. The larger-floor variant progresses much
farther, but remains behind both CVODE runs and is still limited by sulfur
chemistry. The pure ratio policy barely grows its timestep from the initial
stiff transient, leaving its pseudo-time near zero despite thousands of steps.

## Checks, artifacts, and reproduction

The runner verified `jacobian_method=3` in all ten completed native summaries.
All five fresh-start summaries report at least 300 CPU seconds and
`reason=cpu_limit`. Initial-state hashes match across methods within each batch.
Setup checks verify both starts, initial UV optical depth, radii, physical
controls, and the selected analytical Jacobian. The standalone PTC numerical and
atmospheric smoke test also passed.

Current analytical results are under `results/restart/` and `results/fresh/`.
Each contains controls and run metadata, input hashes, per-method traces,
summaries, budgets, final HDF5 states, and a comparison JSON. The progress plot is
`results/progress.pdf`. All generated outputs are ignored by Git. Archived
autodiff data are kept separately; they are not used in the current plot or
comparison JSON.

To repeat from the repository root in the build environment:

```sh
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --reuse-inputs --jacobian-method analytical --modes cvode cvode1 ptc ptc-residual ptc-floor12 --steps 100000 --seconds 45
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --fresh --reuse-inputs --jacobian-method analytical --modes cvode cvode1 ptc ptc-residual ptc-floor12 --steps 100000 --seconds 300
PYTHONPATH=. python integration_benchmarks/compare_ptc.py --case venus
```

The standard runner runs methods sequentially; the completed batches used
independent concurrent processes with the same per-process limits. Omit
`--reuse-inputs` when preparing the case for the first time on a new checkout.

## Interpretation

Venus strengthens the case for improving PTC's global scaling and startup/
acceptance policy before public integration. Change-based PTC with a larger
floor works well for a warm start, but does not solve the fresh-start problem
within five CPU minutes. Pure residual-ratio PTC is ineffective here. A residual
policy with deliberate timestep growth remains untested in these comparisons.
Keep the production default unchanged and avoid treating one tuned warm-start
success as evidence of general robustness.
