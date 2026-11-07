# Photochem steady-state integration with Euler PTC

## Objective

Develop an optional linearly implicit Euler pseudo-transient-continuation (PTC) backend for finding photochemical steady states. The goal is to solve cases where CVODE stalls or oscillates without requiring users to tune absolute tolerance, relative tolerance, or integration order for each atmosphere. A small runtime improvement alone is not the motivation: practical convergence versus failure is the main comparison.

Keep CVODE as the default. PTC is initially a steady-state method, not an alternative accurate transient integrator. Leave `evolve` on CVODE.

This plan expands roadmap item 7.5. No PTC implementation or public API change has been made as part of writing this plan.

## Evidence and references already available

A local benchmark motivated this investigation. Starting from the saved atmosphere, default higher-order CVODE showed persistent numerical chatter and limited time advancement. Order-one CVODE advanced beyond $10^{17}$ s in about 124 accepted steps in the local restart experiment. A separate order-one run with a timestep cap also satisfied the history-based convergence criterion, supporting the interpretation that the chatter was numerical rather than a substantial physical cycle.

Tightening relative tolerance helped the history-based criterion but did not resolve the time-advancement problem in the corresponding time-only experiment. Freezing the prescribed pressure-temperature profile or switching the chemistry Jacobian did not eliminate the problem. These observations justify a prototype; they do not establish that PTC will outperform order-one CVODE.

For this benchmark, use `conv_longdy = 0.0` and `equilibrium_time = 1.0e17`, as requested. Reapply configuration after restoring a model state. Preserve the original checkpoint and run comparisons from copies. Establish tracked, reproducible benchmark inputs before relying on these experiments for long-term comparisons.

Useful references:

- [ATMOS code](https://github.com/VirtualPlanetaryLaboratory/atmos): one linearly implicit Euler correction per attempted step, a banded solve, and change-based timestep adaptation. Its species-specific exclusions, concentration repairs, and absolute-value operation should not be copied.
- `src/clima/clima_ptc.f90`: existing callback interfaces, band-storage conversion, and linear algebra. It is already fairly general, but has more machinery than needed here and uses an unscaled residual norm by default.
- `src/photochem/photochem_evoatmosphere_integrate.f90`: current CVODE stepping, robust-stepper maintenance, convergence, and restart behavior.

Leave the climate PTC implementation unchanged. Start a small independent implementation and borrow useful numerical patterns rather than refactoring both solvers at once.

## Numerical method

For the autonomous system $dy/dt = F(y)$, compute one correction at the current accepted state:

$$
\left(\frac{I}{\Delta\tau}-J(y_n)\right)\delta y = F(y_n),
\qquad y_{n+1}=y_n+\delta y.
$$

Here $\Delta\tau$ is a pseudo-timestep and $J$ is the existing approximate banded Jacobian. This is **linearly implicit Euler**, equivalent to one Newton correction for backward Euler starting at $y_n$. It does not iterate the nonlinear backward Euler equation to convergence. Keep that distinction explicit when comparing it with order-one CVODE.

Use the existing banded linear algebra; do not introduce a dense atmospheric Jacobian. Check the exact matrix storage conversion: the Jacobian callback layout and LAPACK factorization layout need not have the same diagonal row.

### First timestep policy: solution-change control

Start with an ATMOS-style policy because it is simple and avoids choosing a residual norm for timestep adaptation. Measure a dimensionless maximum correction using a consistently defined scale, schematically

$$
E = \max_i\frac{|\delta y_i|}{\max(|y_i|,s_i)}.
$$

The adapter supplies the scales $s_i$ in the units of the evolved variables. In Photochem these should correspond to an abundance floor converted to number density, rather than one number-density floor for every atmospheric layer. The exact floor and growth policy remain experimental choices. This floor is a step-control scale, not a claim of transient integration accuracy or a replacement convergence criterion.

Use the ATMOS policy as the first candidate: reject corrections above 25%, shrink the timestep by 0.7, and increase it by bounded factors when changes are small. Its existing growth factors range from 1.1 to 5. Keep adaptation separate from the linear solve so that the table or a smoother rule can be compared without changing the method.

Do not introduce species- or altitude-specific exemptions. Watch for trace species restricting timestep growth, and for small steps being mistaken for progress. A strict relative-change bound can also impede the transition to Newton-like steps at large pseudo-timesteps; measure this rather than assuming the ATMOS policy is optimal.

### Alternative to compare locally: residual-based control

Compare change-based adaptation with a bounded residual-ratio update, for example

$$
\Delta\tau_{n+1}
=\Delta\tau_n\frac{\|F(y_n)\|_W}{\|F(y_{n+1})\|_W}.
$$

Define and document the residual scaling before using this policy. Avoid an unscaled norm dominated by dense layers or abundant species. Handle zero or nonfinite norms and limit timestep growth and shrinkage. Do not require monotonic residual decrease at every step: legitimate transients can initially increase it.

Neither policy is selected as a universal default yet. Prefer one fixed, understandable policy across the eventual benchmark set over per-case tuning.

### Acceptance and failure behavior

- Check callback status, linear-solve status, finiteness, and trial-state validity.
- Handle unacceptable negative concentrations by reducing the step or timestep. Do not use absolute values or silently repair substantial negative abundances. Define any roundoff allowance explicitly.
- On rejection, restore both the accepted solution and atmospheric work state affected by trial evaluations.
- Bound retries and report a useful failure reason instead of shrinking indefinitely.
- Advance pseudo-time and accepted-step history only after acceptance. Track rejected attempts separately.
- Begin with fresh RHS and Jacobian evaluations at accepted states. Defer reuse optimizations until correctness and profiling justify them.

## Standalone solver and Photochem adapter

### Standalone numerical module

Add the standalone numerical module in `src/photochem/photochem_ptc.f90`. Keep its numerical method independent of atmospheric physics and the EvoAtmosphere adapter. It should own the accepted iterate, pseudo-time, timestep, banded linear-solve workspace, and attempt counters, with initialization, one-step, and cleanup operations.

Accept RHS and banded-Jacobian callbacks plus a small mechanism for user context and step validation/scaling. Keep atmospheric preparation, boundary conditions, chemistry, pressure maintenance, and scientific convergence outside the numerical module. Avoid building a general ODE framework or adding dense-solver support unless needed by a concrete use case. Keep Fortran `err` arguments last.

### Photochem adapter

Add the EvoAtmosphere adapter in `src/photochem/photochem_evoatmosphere_ptc.f90` following the current owner/submodule structure. It should:

1. Call the existing RHS and approximate banded Jacobian.
2. Supply correction scales consistent with Photochem's evolved number densities and gas-density convention.
3. Convert the band layout for the linear solve.
4. Prepare the accepted atmospheric state after a successful step.
5. Restore derived work arrays and mutable profile-dependent quantities after rejected trials. An RHS call can change more than its output vector.

Keep a single authoritative accepted chemical state. Do not duplicate scientific state or serialize numerical workspace into `EvoAtmosphereState`; restarting PTC should initialize fresh solver history.

### Shared robust-stepper behavior

Dispatch the numerical step through the selected backend while sharing TOA maintenance, pressure-temperature-eddy profile resynchronization, accepted-state history, convergence checks, and give-up policy where appropriate. Avoid copying the entire robust-stepper implementation.

Audit CVODE-specific assumptions, especially convergence checks that currently require a CVODE memory pointer and restart paths that initialize CVODE directly. Explicit grid/profile resynchronization must reinitialize the selected backend. Scheduled restarts used for CVODE history should be assessed separately; PTC has no multistep history to discard.

Destroying or switching a backend must release its workspace and invalidate its initialization status cleanly.

## Pseudo-time and convergence

For autonomous steady-state problems, accumulated accepted pseudo-time can support the requested $10^{17}$ s benchmark stopping rule. It is not an accurate physical trajectory or an independent guarantee that $F(y)=0$. Also inspect residuals, final composition, and budgets when validating the method.

Retain the existing scientific convergence options where meaningful and report the stopping reason. A small correction caused by a tiny timestep must not independently imply convergence.

Time-dependent stellar forcing or callbacks require an explicit policy before public support. Do not silently evaluate them at accumulated pseudo-time as though it were accurate physical evolution. The first prototype should use autonomous inputs, with unsupported time-dependent configurations rejected or kept outside its documented scope.

## Proposed public interface, after validation

Add an optional method selector with CVODE as the default:

```python
pc.find_steady_state(method="ptc")

pc.initialize_robust_stepper(pc.wrk.usol, method="ptc")
give_up, converged = pc.robust_step()
```

The initialized backend determines subsequent stepping and destruction; callers should not repeat the method on every step. Validate unknown methods immediately. Decide whether the basic `initialize_stepper` should expose PTC after the robust/steady-state path works; it is not required for the first public pass.

Implement matching Fortran entry points and C/Cython bindings in the repository's existing style. Document which controls apply to each backend. CVODE tolerances and integration-order settings must not appear to control PTC when they do not. Avoid exposing a large new set of tuning knobs before benchmark evidence supports them.

## Implementation and validation sequence

1. **Create benchmarks for testing.** Add reproducible benchmark scripts and tracked inputs, including mechanisms, settings, initial states, and configuration. Start with a case where CVODE stalls or oscillates and establish default CVODE and order-one CVODE baselines. Keep generated outputs separate from input fixtures.
2. **Build the standalone Fortran prototype.** Test the banded correction, rejection/rollback, and bounded timestep adaptation on small stiff systems with known steady states. Include an oscillatory decay problem to check damping and a nonlinear case to exercise rejected trials.
3. **Connect a Fortran-only Photochem adapter.** Restore the benchmark through HDF5 and compare the two timestep policies against both CVODE baselines. Do not wire the public Python API yet.
4. **Broaden local evaluation.** Include a chatter case, an abundance-tolerance-sensitive case, a case with condensation or pressure feedback, and an ordinary case where CVODE already works well. Use one PTC policy and unchanged defaults across cases. Additional difficult atmospheres are needed before claiming broad robustness, not before starting the prototype.
5. **Integrate backend selection.** If the prototype demonstrates a repeatable practical advantage, share robust-stepper maintenance, add method selection and bindings, and test switching, destruction, restarts, and grid/profile resynchronization.
6. **Document the result.** Record the chosen policy, remaining failures, applicability limits, and whether PTC should remain opt-in. Leave CVODE as default unless broader evidence supports changing it.

For each local comparison, record convergence and stopping reason, pseudo-time or physical time as appropriate, accepted/rejected steps, RHS/Jacobian evaluations, linear solves, wall time, and final atmospheric differences. Check elemental budgets and boundary fluxes where applicable. Reaching the time threshold by enormous steps is useful only if the resulting state is also scientifically consistent with a trustworthy solution.

Keep CI coverage small: standalone numerical smoke cases and a bounded Photochem backend smoke test, including default-CVODE behavior and initialization/cleanup. Keep long convergence comparisons and policy experiments in local opt-in benchmark scripts. Do not promote PTC based on a local benchmark alone.

## Deferred work

Fully nonlinear backward Euler, higher-order Rosenbrock methods, accurate PTC transient output, advanced Jacobian reuse, a shared Clima/Photochem solver refactor, and changing the default backend are outside the first implementation. They can be revisited if the simple prototype demonstrates a clear need.
