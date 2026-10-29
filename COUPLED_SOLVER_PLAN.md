# Coupled photochemistry–climate solvers and rollback

## Purpose and current baseline

Develop the coupled model and atmospheric state support together on `plan/evoatmosphere-state-checkpoints`. The existing atmospheric-state design is in `EVOATMOSPHERE_STATE_PLAN.md` on that branch; this document adds coupling architecture and rollback requirements.

The initial implementation alternates photochemical steady-state integration and radiative-convective climate solves with under-relaxation. The local `tests/photochem/test_coupled.py --full` case converged in 10 coupled iterations in about four minutes. Every chemical solve converged, with no integration or resync errors. Default test execution remains constructor-only smoke coverage; the full coupled calculation is local and outside CI.

## Physical model and public API

Keep `CoupledPhotochemistryClimate` as the owner of `pc` (`EvoAtmosphere`) and `c` (`AdiabatClimate`). Both remain accessible for configuration and inspection. Exclude aerosol optics for the initial implementation.

The public methods are:

- `solution()`: export the current photochemical atmosphere without changing it or claiming convergence.
- `initial_guess(...)`: solve initial climate, initialize photochemistry, and return an independent coupled state. This modifies the model and invalidates any previous solver.
- `model_state()`: capture an independent coupled state without modifying the models or claiming convergence.
- `initialize_solver(initial_guess=state, ...)`: restore the atmosphere and initialize the selected coupling algorithm from a coupled state.
- `iterate()`: perform one useful unit of coupling work, returning `(give_up, converged)`.
- `find_steady_state(...)`: initialize and repeatedly iterate until convergence or give-up.

Other methods, including the fixed-point map, grid construction, and profile exchange, are implementation helpers and should use leading underscores. Keep shared prescribed Kzz inputs on the coupling model; avoid duplicating planetary and atmospheric state already owned by `pc` or `c`. Climate condensate pressures use live photochemical fixed-pressure boundaries when present and lowest-layer mixing ratios times surface pressure otherwise. Climate determines condensate vertical profiles.

Use separate solver objects behind `model.solver`:

| Solver | Owned controls and runtime state |
| --- | --- |
| `FixedPointSolver` | Temperature/log-pressure tolerances, relaxation, iteration limit, trial vectors, residuals, and history. |
| `QuasiStaticSolver` | Outer chemical time, coupling interval, timestep-control strategy and tolerances, accepted/rejected-step counts, trial snapshots, and history. |

Both expose common convergence/give-up and chemical-convergence diagnostics. Keep method-specific diagnostics with their solver rather than adding unrelated fields to the coupling model.

Intended usage:

```python
guess = model.initial_guess(edd=1e5)
model.initialize_solver(
    initial_guess=guess,
    method="fixed_point",
    solver_options={"relaxation": 0.5},
)
while True:
    give_up, converged = model.iterate()
    if give_up or converged:
        break
atmosphere = model.solution()
```

Only `method="fixed_point"` is currently implemented. The future `method="quasi_static"` will select adaptive time-evolving coupling. One fixed-point iteration evaluates one chemistry-to-climate map. One quasi-static iteration will complete one accepted outer time interval, including internal retries.

## Coupled warm starts

`initial_guess()` and `model_state()` return the same versioned nested dictionary:

- `schema_version`: currently 1.
- `photochem`: `pc.model_state().to_dict()`.
- `edd`, `edd_pressure`: independent arrays containing original prescribed Kzz knots.
- `climate`: only `pressure` and `temperature` arrays (surface followed by layer centers), and layer `convecting_with_below`. Profile shapes and values are checked against the destination climate grid; climate configuration remains the destination's responsibility.

The climate section is captured and restored through private helpers, leaving room for a future cloud–climate component to own its restart data. The model retains it in `_climate_state`, updated after each successful climate solve or restoration; capture always copies that stored representation. After restoration it holds valid restart profiles even before live `c` outputs are rebuilt. Fixed-point vectors are derived from the stored profiles rather than serving as shared physical state.

```python
state = model.model_state()
other.initialize_solver(initial_guess=state)
# Or restore and run a new solve:
other.find_steady_state(initial_guess=state)
```

Initialization restores photochemistry directly, including its installed profile mismatch, and retains climate profiles as restart guesses without solving climate. Live `c` atmospheric outputs are rebuilt by the first map evaluation; before then, another `model_state()` captures the retained restart data. A new solver resets history and convergence flags and verifies convergence through iteration.

Destination boundary conditions, callbacks, numerical settings, and climate configuration remain in effect. Compatible static mechanisms and optical data are the caller's responsibility; metadata checks do not prove complete physical equivalence. Snapshots exclude time, integrator history, and solver bookkeeping. They support warm restart, not exact trial rollback or identical iteration sequences. File persistence and quasi-static time/rollback guarantees remain later work.

CI runs the smoke test only. Local `tests/photochem/test_coupled.py --restart` adds state round-trip and validation checks with an initial climate solve; `--full` also checks coupled convergence and restarting a converged state in a fresh model.

## State representations and final-result policy

The photochemical model is authoritative for the atmospheric solution. The fixed-point vector is a separate representation: log10 surface pressure and temperatures on a fixed-size, equally spaced log-pressure grid between the surface and climate top. Surface pressure defines temperature coordinates; it is not imposed on the chemical model. Climate top defaults to 1 dyn/cm² and can differ from the photochemical top.

The fixed-point map prepares the proposed T(P) with `set_press_temp_edd_profile`, advances chemistry, solves climate using the resulting composition, and returns the climate vector without applying it back to photochemistry. Use the unrelaxed residual for convergence. Never perform an extra final temperature update that invalidates chemical convergence.

Current chemical defaults are 7,000 accepted steps and `conv_longdy=1e-3`, adjustable through `pc.var`. If chemical convergence is optional, permit continuation only after exhausting the accepted-step budget with usable state. Report chemical and coupled convergence separately. Other chemistry failures remain errors; climate must converge.

## Quasi-static adaptive coupling

The intended approximation is chemistry evolving in physical time with climate equilibrating rapidly to the changing composition:

- Hold climate fixed while advancing chemistry over an outer interval.
- Solve climate at the interval endpoint and update the thermal state.
- Adapt the outer interval using either a change limit or a coupling-error estimate, independently of CVODE's internal chemistry tolerances.

No artificial temperature-relaxation ODE or chosen thermal response timescale is required. This is a quasi-static approximation, not a model of finite thermal inertia.

### Two timestep-control strategies

Both strategies use the same chemistry segments, climate updates, accepted-state snapshots, and rollback machinery. Select the controller within `QuasiStaticSolver`; separate coupling classes are unnecessary.

| Strategy | Acceptance and interval control | Tradeoff |
| --- | --- | --- |
| Change limited | Limit the temperature change over one trial interval; reduce and retry when excessive, and grow the next interval gradually when changes are small. | One chemistry segment and climate solve per trial, but no estimate of coupling error. |
| Step doubling | Compare one full interval with two half intervals; accept and adapt using the endpoint difference. | More expensive, but provides an estimate of coupling error. |

### Change-limited stepping and rejection

1. Capture the accepted physical state and time.
2. Advance chemistry over one trial interval, then solve climate.
3. Compare old and new temperatures on common pressure coordinates using a combined absolute and relative change limit:

   ```python
   change = max(abs(T_new - T_old) / (T_atol + T_rtol * abs(T_old)))
   ```

4. Accept the updated state when `change <= 1`; otherwise restore the accepted state, shorten the interval, and retry.
5. Adjust the next interval using the measured change, with bounded growth and shrinkage. Consider additional surface-pressure or composition limits if temperature alone misses important evolution.

These tolerances limit physical change per interval; they do not estimate timestep accuracy or guarantee stability. Check interval refinement when validating this strategy. Implement this simpler controller first, then reuse its trial and rollback framework for step doubling.

### Step doubling and rejection

1. Capture the accepted physical state and time.
2. Trial one full interval, followed by a climate update.
3. Restore the same starting state and time.
4. Trial two half intervals, with climate updates at the midpoint and endpoint.
5. Compare the two endpoints on common pressure coordinates, using temperature, surface pressure, and suitable composition measures with relative and absolute tolerances.
6. Accept the two-half-step result if the estimated error satisfies tolerance; otherwise restore the accepted state, reduce the interval, and retry.

Establish the method order and corresponding error/step-size controller explicitly.

### Shared safeguards and termination

For both controllers, bound interval growth, rejection retries, and minimum interval size. Recoverable trial chemistry/climate failures can reject a step; invalid configuration and failed restoration are hard errors.

Preserve absolute chemical time across segments and rejections, including the time seen by flux/rate callbacks. Rejected trial time is not accepted progress. Accurate outer time stepping does not itself prove steady state; final termination checks must address both chemical stationarity and climate consistency.

## Atmospheric and climate snapshots

Implement independent in-memory atmospheric capture/restore before the adaptive solver. Preserve composition, current T/Kzz, grid, particle radii, planetary geometry, tropopause, and installed profile/TOA configuration. Regenerate derived quantities, but do not resynchronize the stored P–T–Kzz profile or remove its existing periodic mismatch on restore. CVODE history need not be serialized; restarting it must preserve the outer absolute time.

Reconcile this with the existing checkpoint plan, which deliberately retains destination boundary conditions and callbacks and excludes time/history. Coupling rollback must either preserve mutable boundary/condensation controls or ensure trial steps never change them. Store outer time and any mutable callback state separately; arbitrary external callback side effects cannot be assumed to roll back automatically.

Capture climate state as well, including atmospheric fields and any solver guesses/configuration modified during trial climate solves. The second branch must not inherit unrecorded state from the first. Define failure guarantees honestly and validate that both physical models remain usable after restoration.

HDF5 persistence is a later phase and is not needed for in-memory trial branching. The state design should also accommodate flat NumPy-array dictionary export/import for parameter-grid workflows, using version metadata, independent copies, and explicit validation. Grid aggregation requires consistent keys, shapes, and dtypes across points.

## Future cloud and haze coupling

Add this after the basic solvers and state support are complete. Allow a replaceable cloud–climate component that internally converges climate and cloud/haze microphysics together. Its inputs may include composition, particle production rates from photochemistry, and atmospheric temperature, pressure, altitude, and density. Keep climate-specific exchange behind a few helpers now; defer the detailed component interface.

Proposed constructor API:

```python
model = CoupledPhotochemistryClimate(
    ...,
    climate_model=HazyClimateModel,
)
```

Here `climate_model` is a class constructed with shared model settings. The default implementation and custom components would follow the same small interface:

- `solve(inputs)`: converge climate and any internal cloud/microphysics coupling; report convergence.
- `profiles()`: return atmospheric profiles and optional haze opacity.
- `capture_state()`: return an independent snapshot for rollback.
- `restore_state(state)`: restore a captured state.

Use an adapter around `AdiabatClimate`, allowing a custom component to own both climate and microphysics models. Finalize constructor arguments, method signatures, input/output fields, units, and failure guarantees when implementing this phase.

Initially, microphysics could supply climate opacity only while photochemistry retains its own particle treatment. This approximation can give the two models different haze shielding and material removal. A later option should exchange microphysical opacity with photochemistry so that haze also affects photolysis consistently. Define production units, material accounting, and responsibility for particle transport/removal before implementing that exchange.

The feedback would be:

```text
photochemistry ↔ (climate ↔ clouds/haze microphysics)
```

The full opacity array need not become an explicit fixed-point variable. Quasi-static coupling can exchange temperature and opacity as model state between chemical time segments. For a fixed-point solve, however, temperature and surface pressure alone do not define the map if it also depends on opacity left over from the previous evaluation. Resolve that feedback through an inner haze–photochemistry iteration or additional compact variables, such as a haze-production profile; decide after testing. Include haze consistency in convergence checks.

The replaceable component must expose its resulting profiles and convergence status, and support capture/restore of microphysics state, opacity, and solver guesses for rejected outer timesteps. Its internal solver need not be prescribed by the outer coupling API.

## Implementation order

1. Implement and test in-memory `EvoAtmosphereState` capture/restore and the absolute-time restart path needed by coupling. Keep expensive integrations out of routine CI.
2. Separate the existing fixed-point algorithm into its solver object while preserving its tested behavior and public coupling API; rename implementation helpers with underscores.
3. Implement climate snapshot/restore and a combined accepted-state transaction for trials and rollback.
4. Add bounded chemistry time segments and fixed outer-interval quasi-static coupling first.
5. Add change-limited adaptive intervals and acceptance/rejection first, then add step doubling with scaled error measurement using the same framework.
6. Verify rollback, time-dependent callbacks, behavior under interval refinement, and final chemical/climate status for both controllers. Compare representative solutions and cost with the existing fixed-point method.
7. Add HDF5 coupled-state persistence and flat dictionary interchange separately. The initial public coupled warm-start API is described above; stronger rollback guarantees remain separate.
8. Add the replaceable cloud–climate component and cloud/haze feedback options described above, including rollback and convergence checks for the additional state.
