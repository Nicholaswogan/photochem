# Coupled hydrocarbon haze and prescribed surface temperature

## Purpose and scope

Extend `photochem/extensions/coupled.py` with two composable physical features:

1. External steady-state hydrocarbon haze microphysics, exchanged with chemistry and climate through a single-sweep outer fixed point.
2. A prescribed surface temperature, achieved by adjusting the photochemical CO2 pressure boundary.

These give four supported modes, sharing initialization, state handling, profile exchange, and numerical solvers. Implement these concrete cases rather than a general interface for arbitrary climate-model plugins.

This document is a future implementation plan, not a description of completed features. It supplements [COUPLED_SOLVER_PLAN.md](COUPLED_SOLVER_PLAN.md) and supersedes that document's earlier generic cloud–climate component proposal. The existing no-external-haze fixed-point solver remains the baseline. `QuasiStaticSolver` remains deferred.

Initial microphysics follows the simple McKay coagulation–sedimentation approach: one characteristic particle size per pressure level, steady state, and no evaporation. Sulfuric acid clouds, S8 aerosols, fractal aggregates, and time-dependent particle evolution are later extensions.

The first implementation uses one chemical steady-state calculation, one climate update, and one microphysical calculation per outer evaluation. It does not require an internally converged climate–microphysics fixed point. A separate operation can later repeat the same climate–microphysics sweep to equilibrium for quasi-static time stepping. Climate is substantially cheaper than chemistry, so repeated climate evaluations are acceptable in the prescribed-temperature scalar search.

## Four modes

| Mode | Climate update | Feedback supplied to photochemistry |
| --- | --- | --- |
| Haze-free | Ordinary radiative-convective equilibrium | Temperature profile and surface-pressure coordinate |
| Haze-free, prescribed surface T | Solve CO2 pressure to obtain the target surface temperature | Temperature profile, surface-pressure coordinate, and CO2 pressure boundary |
| Haze-coupled | One climate solve with trial haze optics, followed by one microphysical solve | Temperature, surface-pressure coordinate, and external haze optical properties |
| Haze-coupled, prescribed surface T | Solve CO2 pressure with trial haze optics fixed, followed by one microphysical solve | Temperature, surface-pressure coordinate, external haze optical properties, and CO2 pressure boundary |

Surface pressure is exchanged in all four modes to define the prescribed profile coordinates. It is not imposed as a total-pressure boundary condition on photochemistry; chemistry determines its own surface pressure, which must become consistent with the exchanged value at convergence.

“Haze-free” means no new externally coupled hydrocarbon haze. Existing photochemical particles and their chemical/transport treatment remain available. Native optical contributions must be controlled explicitly when external haze optics represent the same material.

Use one climate–microphysics sweep with these explicit branches. Keep the physical formulation separate from numerical iteration: named blocks define the exchanged unknowns, and the selected solver proposes their next values. Avoid four independent coupling implementations. The shared operation always runs climate before microphysics, so microphysics receives the climate-adjusted condensable composition as well as consistent temperature, density, molecular weight, and geometry.

## Shared coordinates and units

Let `nc` be the climate layer count and `M = nc + 1`. Retain `M` samples for each coupling profile initially.

The coupling pressure coordinates are

```python
s = log10(P_surface)
u = linspace(s, log10(target_pressure), M)
P = 10**u
```

Pressure is in dyn/cm². `target_pressure` is the prescribed photochemical TOA pressure, not the current approximately maintained top pressure. Climate retains its own layer grid and top pressure. Pressure nodes are derived coordinates, not independent vector entries.

All coupling profiles use these coordinates. Outputs from native climate, photochemical, or haze grids must be mapped to the coupling grid. A profile must always be interpreted with the surface pressure that defines its coordinates. Remapping between pressure intervals must preserve the intended integrated source or particle column; changing the modeled surface boundary is not merely geometric stretching.

Temperature is in K. Particle radius `r` is in cm. Define the exchanged particle abundance as

```text
A(u) = particle column per unit log10-pressure
units: particles/cm² per pressure decade
```

For a downward pressure interval, its particle column is `integral A(u) du`. The volume number density needed by microphysics is derived from the current atmospheric geometry:

```text
n_particle = A(u) / abs(dz/du)
```

This preserves particle column under stretching/compression between fixed pressure boundaries. It does not prevent physical changes from production, coagulation, or settling when microphysics is solved again.

Use hydrostatic pressure as the default coordinate for shared prescriptions, consistent with climate and the coupled model's default persistent thermal-profile setting. A standalone photochemical setter may support gas-pressure mapping as an explicit alternative. Do not silently interpret a gas-pressure source profile as hydrostatic pressure; transform through its supplied atmospheric geometry when needed.

## Exact outer fixed-point variables

Use named blocks and explicit pack/unpack helpers. For positive radii, define `R_i = log10(r_i)`. For the CO2 boundary, define `c = log10(P_CO2_BC)`, with pressure in dyn/cm². Store particle column `A_i` directly with numerical scaling so that a genuinely haze-free profile can be exactly zero.

`T_0` is the surface temperature; indices `1` through `M-1` are the remaining coupling samples. In prescribed-temperature modes, `T_0 = T_target` is imposed on the trial profile and is not an independent unknown. The inverse climate solve must actually satisfy the target; overwriting an unconverged climate surface temperature with the target is not a valid solve.

| Mode | Packed outer vector, in order | Length |
| --- | --- | ---: |
| Haze-free | `[s, T_0, ..., T_(M-1)]` | `M + 1` |
| Haze-free, prescribed surface T | `[s, T_1, ..., T_(M-1), c]` | `M + 1` |
| Haze-coupled | `[s, T_0, ..., T_(M-1), R_0, ..., R_(M-1), A_0, ..., A_(M-1)]` | `3*M + 1` |
| Haze-coupled, prescribed surface T | `[s, T_1, ..., T_(M-1), R_0, ..., R_(M-1), A_0, ..., A_(M-1), c]` | `3*M + 1` |

For the initial single-population model, radius and column are one-dimensional profiles. Generalization to multiple materials/populations is deferred.

The following are not additional outer unknowns:

- Pressure-grid nodes: derived from `s` and `target_pressure`.
- Spectral optical properties: derived from particle profiles and fixed optical tables.
- Gas composition, density, molecular weight, and chemical production: outputs of chemistry and climate/profile construction.
- Altitude and layer thickness: reconstructed from atmospheric structure.
- Haze production profile and its column integral: outputs of each chemical evaluation, fixed during its climate–microphysics sweep.
- Climate convective classification: a warm-start guess, not a physical coupling unknown.

The absence of these quantities from the vector requires that each component closes its own dependencies. Required inputs must not be taken silently from an unrelated previous evaluation.

## Single-sweep outer fixed point

The map is a sequential update, with one outer solver handling all thermal, particle, and optional CO2 feedback:

```python
trial_opacity = haze_opacity(x.particles)  # Omitted in haze-free modes.

composition, production = photochemistry(
    thermal_profile=x.thermal_profile,
    co2_boundary=x.co2_boundary,          # Only a solved control in fixed-T modes.
    opacity=trial_opacity,
)

atmosphere_new = climate_update(
    photochemical_inputs=composition,
    opacity=trial_opacity,
    target_temperature=T_target,         # None for free surface temperature.
)

particles_new = microphysics(             # Omitted in haze-free modes.
    atmosphere=atmosphere_new,
    production=production,
)

return pack(atmosphere_new, particles_new)
```

In detail:

1. Unpack the trial thermal profile, surface-pressure coordinate, optional particle profiles, and optional CO2 boundary control. Reconstruct the prescribed surface temperature when applicable.
2. Set the trial CO2 pressure boundary before chemical integration in fixed-T modes. The trial total surface pressure defines thermal/particle coordinates; it is not imposed on chemistry.
3. Derive the trial particle optical prescription and install it in photochemistry. Suppress the matching native particle optical contribution while retaining the intended native chemical and transport processes.
4. Apply the thermal/Kzz prescription and converge chemistry once. Build dry-gas profiles, surface constraints, and the pressure-based production prescription from that state.
5. Run climate to radiative-convective equilibrium with the trial haze optical prescription fixed. In fixed-T modes, replace this single climate calculation with the climate-only CO2 root search described below.
6. Construct the climate-adjusted atmosphere over the entire haze domain, including the extension above climate top. Run steady McKay microphysics once, using the chemical production prescription from step 4.
7. Return the new thermal profile, surface pressure, optional particle profiles, and optional CO2 control on the appropriate coupling coordinates. Check all active unrelaxed residual blocks and form the next damped/Anderson proposal only if another iteration is needed.

Climate may reshape condensable gas profiles. Microphysics uses that adjusted gas composition, not photochemical molecular weight or vapor profiles left over from a different thermal state.

The new particles do not trigger another climate solve in this evaluation. They supply optics to chemistry and climate on the next outer evaluation. Likewise, the newly solved CO2 control is applied to chemistry on the next evaluation. At a fixed point, these input/output differences vanish to the selected tolerances, so every process sees consistent exchanged fields.

Each free-temperature haze evaluation costs one chemical steady-state solve, one climate solve, and one microphysical solve. Each prescribed-temperature haze evaluation costs one chemical solve, several climate solves, and one microphysical solve. Here “one climate solve” means convergence of its existing radiative-convective calculation; “one microphysical solve” means the complete steady McKay column calculation.

Photochemistry remains authoritative for the exported atmospheric solution. Do not apply an extra final thermal, opacity, or CO2 update that invalidates chemical convergence. Export the last chemical atmosphere with its residual diagnostics rather than silently replacing it with an unevaluated outer proposal.

## One fast sweep versus a converged fast subsystem

Use two operations with distinct meanings:

```python
climate_microphysics_sweep(...)
# One converged climate update, followed by one steady microphysical solve.

solve_climate_microphysics(...)
# Repeat the same sweep until the fast subsystem is self-consistent.
```

`FixedPointSolver` uses the first operation once per chemical evaluation. Implement this path first. It has no additional climate–haze iteration limit or convergence requirement hidden inside the outer map; the outer solver resolves the combined feedback.

The second operation is intended for future quasi-static endpoint updates. Hold the endpoint chemical exchange inputs and production fixed while repeating the sweep. Its optional inner vector is

```text
h = [s, R_0, ..., R_(M-1), A_0, ..., A_(M-1)]
```

Include the CO2 control `c` in that fast-state consistency check when prescribed temperature is enabled. Temperatures and climate-adjusted composition are derived outputs of each climate update. In a haze-free free-temperature case, one climate solve already completes the fast update.

Specify separate tolerances and iteration limits when implementing this converged wrapper. It may use damping/acceleration, but must reuse the same physical sweep rather than introduce another version of the model. Comparing single-sweep and inner-converged strategies later is an optimization/validation exercise, not a prerequisite for the initial fixed-point implementation.

## Prescribed surface temperature: climate-only CO2 search

After the outer chemical solve, freeze the trial haze optical prescription and the chemical exchange inputs other than the released CO2 constraint. Define

```text
f(c) = RCE_surface_temperature(c; fixed_trial_optics, fixed_non_CO2_inputs) - T_target
c = log10(P_CO2_BC)
```

Use a bounded, preferably bracketed scalar root solve. Each scalar evaluation converges climate only. It does not rerun chemistry or microphysics. Keep the same optical prescription at the same absolute pressure coordinates across all scalar trials; changing the trial climate surface pressure must not redefine or slide that prescription onto a new normalized grid.

Once a root is accepted:

1. Retain or reconstruct the actual climate atmosphere corresponding to that root, including its adjusted condensable gas profiles.
2. Run microphysics once with that atmosphere and the production prescription from the completed chemical solve.
3. Return the new particles and CO2 control as outer outputs. The next outer iteration handles their effect on chemistry and climate.

The new particles may change the temperature produced at the same CO2 control on the next evaluation. This is expected outer feedback, not a reason to run microphysics inside the scalar search. When the outer state converges, the optics used by the root search agree with those derived from the returned particles.

The climate calculation releases the old CO2 constraint while retaining the other supplied constraints. It must not also fix the old CO2 vertical abundance profile. The CO2 control must have the same physical meaning in climate and the photochemical lower boundary: distinguish gas partial pressure from a reservoir-equivalent pressure when condensation is active. Explicitly resolve that mapping rather than assuming the existing climate `P_i` argument always has the desired boundary interpretation.

Check the actual computed surface-temperature mismatch independently of the outer profile residuals. `T_0 = T_target` is imposed only when constructing the trial chemical profile; it must not hide a failed inverse climate solve. Reuse valid climate warm starts, but return the state at the accepted root, not whichever scalar trial executed last.

No numerical formulation guarantees a root for every setup. Report inability to bracket a target or a failed climate calculation explicitly. Multiple roots and runaway atmospheres are physical possibilities; the first implementation is a practical strategy, not a claim of universal convergence.

## Custom optical-property prescriptions: first implementation stage

Provide similar Python/Fortran capabilities in climate and photochemistry. Exact parameter names can be decided during implementation; specify the interpretation explicitly and retain existing per-cm behavior by default.

| Prescription | Input quantity | Conversion to layer optical depth |
| --- | --- | --- |
| Per altitude | Extinction coefficient, cm⁻¹ | Integrate over altitude; existing midpoint approximation uses `alpha * dz` |
| Per log10-pressure | `d tau / d log10(P)`, dimensionless per pressure decade | Integrate over the layer's log10-pressure interval |

Inputs include increasing wavelength in nm, decreasing pressure in dyn/cm², and arrays shaped `(n_pressure, n_wavelength)` for extinction/optical-depth density, single-scattering albedo, and asymmetry. Define validation, copying, interpolation, extrapolation, and unset behavior consistently. In particular, explicitly choose the optical prescription outside its supplied pressure range; do not inherit an undocumented interpolator default.

Use pressure-based optical-depth injection for aerosol coupling. For a single population,

```text
d tau / du = Qext(r, wavelength) * pi*r² * A(u)
u = log10(P)
```

`w0` and `g0` come from the same particle optical model. Geometry is not needed by this conversion when `A` is already a particle column per log10-pressure.

Integrate pressure prescriptions onto receiving layer boundaries so grid splitting/merging preserves total optical depth. Use extinction-weighted scattering albedo and scattering-weighted asymmetry for layer averages. Supply actual pressure widths/boundaries through climate's radiative pipeline: its duplicated midpoint pressures cannot be differenced to recover radiative layer widths.

Photochemistry must retain the prescribed profiles separately from the sampled layer arrays and evaluate them during atmospheric preparation, including after TOA regridding. Retain the existing approximate banded Jacobian; do not finite-difference all of preparation or introduce microphysical updates during RHS calls. Define the custom-channel precedence relative to legacy layer arrays and explicitly control native aerosol optical contributions.

Include enough prescription data in state support to restore or reinstall the same pressure-based fields before atmospheric preparation needs them. Warm restart and future rollback must not silently depend on destination opacity left over from a previous run.

Acceptance checks: unchanged per-cm results, optical-column invariance under stretching at fixed pressure boundaries, conservative pressure regridding, correct wavelength mapping, and usable chemical integration with smooth and sharply varying prescribed haze.

## Standalone McKay microphysics: second stage

Implement a public function with private numerical helpers; no class is required for this stateless calculation. Reference the existing code in:

- `/Users/nicholas/Documents/Research_local/NPP/RaviHaze/calcs/haze.py`.
- `/Users/nicholas/Documents/Research_local/NPP/RaviHaze/calcs/codex_reference/McKay`.

Specify inputs before porting equations:

- Host pressure and temperature profiles, atmospheric number density and molecular weight, gravity, and altitude/layer geometry, with CGS units explicit.
- A nonnegative production shape per unit log10-pressure plus its column mass-production rate in g/cm²/s.
- Particle material density, seed/top-boundary assumptions, charging parameter, and viscosity prescription.
- Internal-grid resolution and optional upper extension.

For normalized source shape `w(u)`, derive volumetric production from the current geometry:

```text
integral w(u) du = 1
q_mass(z) = column_production * w(u) * abs(du/dz)
```

Implement this conservatively in pressure intervals. At fixed production inputs, changes in temperature and geometry must not change integrated production. Do not renormalize a truncated source in a way that silently relocates production below a moving surface.

Retain an internal haze grid independent of the host/coupling grids. It may extend above the host domain to place a specified seed boundary above the production region. Define temperature/composition extension and derive density/altitude consistently; do not extrapolate atmospheric fields independently into an inconsistent state. Additional space with zero production does not recover an omitted source above the chemical domain.

Return both internal-grid diagnostics and host-grid particle profiles, including pressure, radius, number density, mass density, and downward mass flux. Convert host results to particle column per log10-pressure for coupling.

Derive the mass and particle-number assumptions explicitly, including how production supplies seeds or particle mass. Reconcile collision-kernel versus number-loss coefficient conventions and the original/Python transport formulas rather than assuming they are identical. Handle zero production and the source-region upper boundary without using arbitrary floors as physical production.

Acceptance criterion: roughly reproduce the reference Titan radius and number-density profiles with documented input choices. Also verify mass-flux balance, nonnegative outputs, grid convergence, and sensitivity to the seed boundary and upper extension. Do not claim an exact reproduction from visually similar profiles alone.

## Table-based haze optics: third stage

Implement a dedicated class whose constructor loads the selected precomputed particle optical tables from `data/` completely into memory. Runtime evaluation performs interpolation only, not new Mie calculations.

Its evaluation accepts pressure-based particle radii and columns and returns the optical prescriptions required by climate and photochemistry. Support the two models' wavelength grids while using the same underlying particle material and table conventions. Keep this class separate from the microphysical equations.

Before choosing tables, confirm material identity, shape/size-distribution assumptions, radius definition and units, spectral coverage, and the meaning of stored extinction/scattering efficiencies. Define explicit out-of-range behavior. Do not assume spherical Mie tables also represent fractal aggregates.

Verify known table nodes, interpolation, units, zero-haze behavior, and the identity between derived extinction times geometry and pressure-integrated optical depth.

## Material bookkeeping and model limitations

For the initial coupling, photochemistry owns gas/aerosol chemistry, condensation, evaporation where already present, and material transport/budgeting. External McKay microphysics diagnoses the hydrocarbon particle population used for radiation. It is solved outside CVODE's RHS and has no history-dependent radius update during chemical trial evaluations.

Disabling native haze optics does not disable native chemical reactions or transport. Do not overwrite native particle radii merely to install external optics: that would also change the photochemical particle conversion and settling calculations.

The external radiative population can therefore differ from the native chemically transported population. This is an explicit approximation. Record differences in aerosol mass columns/profiles and determine whether they are acceptable. Do not claim exact agreement of the two aerosol representations merely because native chemical elemental budgets close.

## Residuals, numerical controls, and diagnostics

Keep physical blocks separate from their packed numerical representation. Retain unrelaxed residuals, plain damping by default, and optional Anderson history confined to the numerical solver.

Define separate tolerances/scales for:

- Relative/absolute temperature consistency.
- Total log10-surface-pressure consistency.
- Log10-CO2 boundary consistency when enabled.
- Radius consistency where particles materially contribute.
- Particle-column consistency, with an absolute tolerance allowing zero/trace haze.
- Actual target-temperature satisfaction in prescribed-temperature modes.

Specify how profile residuals are compared across differing trial/output pressure coordinates. Radius in a particle-free region is not a meaningful convergence condition. Bound proposals to valid pressures, positive radii/temperatures, and nonnegative columns. Do not reuse the current assumption that every packed entry after surface pressure is a temperature.

Report outer residuals alongside chemical convergence, climate convergence, scalar-root status, and evaluation counts. Include an inner climate–haze convergence status only when using the future converged fast-subsystem wrapper. Distinguish failed physical solves from exhausting an outer iteration budget. Each outer evaluation must use the supplied trial particles and CO2 control, not hidden mutable values from a previous call.

Test the update ordering explicitly: chemistry and climate use the trial particle optics; microphysics uses the returned climate atmosphere; no extra climate solve follows the microphysical update. In fixed-T modes, scalar root trials must not call chemistry or microphysics, and their haze prescription must remain unchanged. Numerical tests should count these component calls as well as checking residuals and physical bounds.

## State, initialization, and later quasi-static coupling

Preserve the public workflow: `initial_guess()` constructs a state; `model_state()` captures it; `initialize_solver(initial_guess=...)` and `find_steady_state(initial_guess=...)` start fresh numerical sessions.

Extend the versioned state to capture the selected mode, relevant thermal/particle profiles and coordinates, the prescribed temperature constraint, and the solved CO2 boundary control. Preserve optical prescriptions or deterministically regenerate them from captured particles and compatible optical configuration. Aim for NumPy-array leaf values, including scalar schema metadata, while retaining nested physical organization.

General destination boundary conditions remain destination settings. The solved CO2 pressure is a special coupled control and must be explicitly restored/applied. Validate material/table compatibility and all active profile shapes. Reset Anderson history, counters, and convergence flags; do not confuse a warm start with exact integrator continuation.

Future `QuasiStaticSolver` will advance chemistry over an outer physical time interval, then use `solve_climate_microphysics` to converge the fast subsystem at the endpoint. The wrapper repeats the same sweep used by `FixedPointSolver`; it does not require another physical implementation. In fixed-T mode, each repeated sweep uses a climate-only CO2 root search with that sweep's haze prescription fixed, followed by one microphysical update. Endpoint acceptance requires consistency of the final particles, climate, and optional CO2 control.

Quasi-static stepping additionally needs exact accepted-state rollback, including particle/component state, installed opacity, and CO2 boundary control, plus absolute chemical time. The steady McKay assumption treats particle adjustment as rapid; finite-time aerosol evolution is a separate extension. Prescribing surface temperature by adjusting CO2 is a boundary-control assumption, not a predicted carbon-cycle evolution. A single fast sweep at a time endpoint would be an additional approximation requiring interval-refinement validation, not the initial quasi-static design.

## Work order and validation policy

1. Implement and validate custom optical prescriptions in climate and photochemistry, including pressure-coordinate injection and native-opacity controls.
2. Implement standalone McKay microphysics and its Titan comparison.
3. Implement the in-memory table-based haze optical-properties class.
4. Introduce named coupling blocks, expanded state handling, and single-sweep haze coupling: chemistry, then climate with trial optics, then microphysics. Preserve baseline behavior and benchmark damping/Anderson locally.
5. Add prescribed surface temperature for haze-free and haze-coupled modes. The CO2 search repeats climate only, holding trial haze optics fixed; microphysics runs once after the accepted root. Validate CO2 boundary consistency, target residuals, and scalar-root failures.
6. Validate all four single-sweep fixed-point modes locally. Later, add the converged fast-subsystem wrapper for quasi-static endpoint updates, followed by adaptive stepping and exact rollback. Benchmark optional repeated fast sweeps for steady-state solving only if useful.

CI remains smoke-only. Keep deterministic local tests for interpolation/conservation, particle/opacity reconstruction, state independence/restoration, and numerical failure paths. Keep Titan comparisons, coupled integrations, four-mode physical examples, and timing benchmarks explicitly opt-in. Record approximation choices and test configurations so comparisons remain reproducible.
