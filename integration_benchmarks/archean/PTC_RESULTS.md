# First Euler PTC experiment

## Outcome

The Fortran-only linearly implicit Euler prototype works and gives a useful
saved-state result. It does **not** yet establish a generally better solver than
CVODE. Keep CVODE as default and defer public backend selection.

The ATMOS-style change-based policy was substantially more effective than a
pure residual-ratio policy on the saved-state case. However, the same default
PTC policy performed poorly on a climate-generated fresh atmosphere. The
step-control policy and scaling still need work before broader adoption.

## Implementation exercised

- `src/photochem/photochem_ptc.f90`: banded correction, change-based or
  residual-ratio timestep adaptation, bounded rejection, and callback rollback.
- `src/photochem/photochem_evoatmosphere_ptc.f90`: existing RHS and analytical
  banded Jacobian, gas-density-based scaling, accepted pseudo-time, and shared
  pressure/profile maintenance. Python bindings are unchanged.
- `tests/photochem/test_ptc.f90`: small numerical checks and a bounded atmosphere
  test with a real TOA grid transaction.

Both timestep policies retain a maximum 25% scaled correction. PTC rebuilds and
factors the banded Jacobian each step; there is no reuse optimization. Trace
negative projections are explicitly bounded and counted. The default abundance
floor is `1e-20`; the allowed trace projection is at most approximately `1e-24`
in abundance for a nonnegative state. The largest scaled projection in the
saved-state default run was `8.31e-5`, below the `1e-4` allowance.

## Saved-state comparison

All methods started from the same saved state and physical configuration, with
`conv_longdy=0`, `equilibrium_time=1e17`, and up to 7000 accepted steps. Runs used
one OpenMP thread. Times below are integration CPU seconds, excluding setup of
the model, compilation, and final budget/output calculations; solver
initialization is included. These are individual local runs, not precision
performance measurements.

| Method | Reached time threshold? | Accepted steps | CPU seconds | Final time / pseudo-time (s) |
|---|---|---:|---:|---:|
| Default CVODE | Yes | 4074 | 51.24 | 1.0005e17 |
| Order-one CVODE | Yes | 124 | 1.96 | 1.6841e17 |
| Change-based PTC | Yes | 47 | 2.49 | 3.8767e17 |
| Residual-ratio PTC, bounded to 500 steps | No | 500 | 24.19 | 53.94 |

The full default-CVODE run matters: shorter runs appeared stuck, but this restart
**did eventually converge under the requested time criterion**. Do not describe
this result as default CVODE failing within 7000 steps.

Default PTC used 48 linear solves, 47 Jacobian evaluations, 95 RHS evaluations,
and one rejected attempt. It needed fewer accepted steps than order-one CVODE,
but more Jacobian factorizations and slightly more CPU time. Order-one CVODE is
therefore still a very competitive and simpler solution for this restart.

The two successful PTC/CVODE1 states were close:

- Maximum relative column difference for species exceeding a `1e-10` column
  fraction: approximately `5.45e-4` (0.055%).
- Maximum relative profile difference above a `1e-20` gas mixing ratio:
  approximately `2.88e-3` (0.29%). Above `1e-6`, approximately 0.15%.
- Maximum temperature difference: approximately `2.1e-6` K.
- Maximum absolute mixing-ratio tendency: `2.04e-17 /s` for PTC versus
  `2.58e-13 /s` for order-one CVODE.
- Normalized column RHS/production-loss reconciliation error: `2.01e-10` for
  PTC versus `2.75e-9` for order-one CVODE.

The time criterion alone does not imply equal residual accuracy, and neither
solution is an independently exact reference. Column-budget reconciliation
checks that the diagnostics reconstruct the RHS; it is not an independent proof
of elemental conservation. Elemental column tendencies are reported separately
for this open atmosphere.

A higher step-control floor (`1e-12`) also succeeded from the saved state: 45
steps and about 2.32 CPU seconds. An additional order-one CVODE sensitivity run
with `rtol=1e-6`, `atol=1e-23` did not reach the time threshold in 1200 steps.
Neither experiment is a reason to tune each individual atmosphere.

## Fresh-start comparison

All methods started from the same climate-generated atmosphere. Each run was
bounded to approximately 60 integration CPU seconds; **none reached the time
threshold** during that budget.

| Method | Accepted steps | Final time / pseudo-time (s) |
|---|---:|---:|
| Default CVODE | 3898 | 4.55e9 |
| Order-one CVODE | 4776 | 5.91e6 |
| Default change-based PTC | 1131 | 4.03e3 |
| Change-based PTC, abundance floor `1e-12` | 1088 | 1.66e5 |

PTC initially limits steps on newly produced trace species, then on aerosol
concentrations during condensation. A larger abundance floor helps early
progress but does not fix the overall result. The default policy is much less
effective here than CVODE. These bounded runs do not establish eventual
convergence or failure; they do establish that the saved-state advantage does
not automatically carry over to a fresh start.

## Maintenance and correctness checks

Two deliberately perturbed saved-state experiments exercised actual atmospheric
transactions:

- Periodic profile mode with a 10% change in prescribed temperature: one
  resynchronization, then continued stepping to about `4.20e5` s within a
  35-second CPU budget.
- A changed TOA target with tighter pressure thresholds: one grid
  resynchronization, then continued stepping to about `1.72e7` s within the same
  budget.

Neither reached the time threshold within its budget. These runs exercise
maintenance compatibility, not full coupled convergence. An initial version
that merely tolerated tiny negative trace densities could become trapped at its
positivity threshold. Bounded trace projection removed that failure in the
periodic run; significant negative trials still reject.

The focused numerical smoke checks cover a known banded linear solution, stiff
oscillatory damping, nonlinear relaxation, exhausted retries, failed trial
callbacks and rollback, and the trace-projection allowance. The atmospheric
smoke checks TOA maintenance and preservation of accepted pseudo-time. Native
`test_api`, `test_jacobian`, and `test_production_loss` also passed during this
pass. Long solver comparisons remain outside CI. The original `state1.pkl` was
not modified (SHA-256
`ce9e3a2dafc5b7627c31827e0a72b1bae611fed4b045f0366242cbd07e3eb87d`).

## Recommendation

Keep the prototype for further experiments, with the change-based policy as the
initial candidate. Do not add a public `method="ptc"` interface or change the
production default yet. The next useful numerical work is a better global
acceptance/growth policy for fresh atmospheres, tested across both starting
conditions with unchanged controls. A pure residual-ratio update without an
extra growth mechanism is not a good choice for this case.

Before public integration, add history-based convergence support, a fuller
backend lifecycle contract, and additional difficult atmospheres. The current
Fortran adapter is a time-threshold experimental driver, not a finished
replacement for the robust CVODE stepper.
