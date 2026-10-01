# EvoAtmosphere atmospheric state and HDF5 persistence

## Status and purpose

This is a future implementation plan. This branch records the design only; it does not implement the API.

Provide a lightweight atmospheric state that can initialize a fresh or existing `EvoAtmosphere` without loading a legacy atmosphere text file or rebuilding the atmospheric composition from a climate profile. Support both an independent in-memory state object and an HDF5 file readable from Fortran and Python.

The state is an atmospheric snapshot, not a full model or CVODE checkpoint. Static reaction and optical data are loaded by the destination constructor. General numerical settings, callbacks, integration history, time, and convergence counters are outside the saved state.

## Intended Python API

```python
from photochem import EvoAtmosphere, EvoAtmosphereState

pc1 = EvoAtmosphere(
    "reactions.yaml",
    "settings.yaml",
    "star.txt",
    "atmosphere.txt",
)
pc1.var.max_resync_failures = 5
pc1.find_steady_state()

state1 = pc1.model_state()
state1.save_to_h5("my_state.h5")

state2 = EvoAtmosphereState("my_state.h5")
pc2 = EvoAtmosphere(
    "reactions.yaml",
    "settings.yaml",
    "star.txt",
)
pc2.initialize_from_state(state2)
```

`pc2` receives the saved atmospheric fields and persistent-profile configuration. `pc2.var.max_resync_failures` retains the destination's value. A new stepper must be initialized before taking integration steps. The same restore operation also accepts `state1` directly, without a file.

Prefer the single restore name `initialize_from_state`; the exploratory `initialize_to_state` spelling in `example.py` is not a second operation.

## State contents

Introduce a public Fortran type named `EvoAtmosphereState`. Keep the existing private `AtmosphereState` as the internal container for candidate construction, derived arrays, and rollback. The new type owns copies of its data and contains no pointers into a live model.

| Saved field | Purpose |
| --- | --- |
| `bottom_atmos`, `top_atmos`, `z`, `dz` | Preserve domain boundaries and layer geometry in cm. |
| `temperature` | Current layer temperatures in K, including any periodic-profile mismatch. |
| `edd` | Current layer eddy diffusion in cm²/s. |
| `particle_radius` | Current particle-radius profiles in cm. |
| `usol` | Evolved species densities in the existing model convention, including numerical trace values. |
| `trop_alt` | Current tropopause altitude in cm. |
| `planet_radius` | Radius used by the atmosphere; gas-giant initialization can change this from the constructor value. |
| `press_temp_edd_profile` | Enablement, mode, hydrostatic-pressure choice, tropopause pressure, prescribed P-T-Kzz arrays, and mismatch tolerances. |
| `toa_pressure_maintenance` | Enablement, target pressure, and normal/extreme pressure factors. |
| Compatibility metadata | Species names/order, species and particle counts, layer count, and planet mass. |

Do not save `grav`, `trop_ind`, `xs_x_qy`, `particle_xs`, or `gibbs_energy`. Recompute these from the saved fields and destination static data. Prepared densities, pressures, transport coefficients, chemical rates, and other workspace arrays are also recomputed.

General robust-stepper settings such as `max_resync_failures`, restart intervals, and convergence limits remain destination settings. The stored TOA and prescribed-profile configuration belong to the atmospheric state because they define how that atmosphere is maintained.

## Restore semantics and compatibility

- Require the same mechanism and compatible static configuration. Check species ordering, particle ordering/counts, array dimensions, and planet mass explicitly. Species names alone do not prove that two reaction mechanisms or optical databases are identical; document that limitation rather than claiming complete compatibility detection.
- Initially require the destination layer count to match. Regridding to a different layer count is a separate feature.
- Preserve current temperature, Kzz, and composition directly. Do not normalize or interpolate the saved densities merely to construct an initial guess.
- Prepare the restored state with `KeepCurrentProfile`. Restoring a periodic profile must not perform a fresh P-T-Kzz synchronization. Continuous synchronization resumes through the normal preparation/RHS policy after restore.
- Keep destination boundary conditions and callbacks. Preparation applies its usual boundary and clipping rules, so identical atmospheric results require compatible boundary conditions. Do not promise that the entire destination model is identical to the source.
- Preserve saved numerical trace abundances using the existing preparation conventions; do not introduce a new positivity floor or reject every negative trace solely for serialization.
- Snapshotting and saving must not modify the source model, its integrator, or its profile configuration.
- A successful restore leaves no active ordinary or robust stepper. Restart integration using the existing initializer APIs.
- Missing required fields, unsupported schemas, or incompatible states fail explicitly. Do not add silent legacy defaults.

## Fortran implementation

### Pass 1: State type and in-memory restore

1. Define `EvoAtmosphereState` and public interfaces alongside `EvoAtmosphere`, with implementation details in the initialization submodule or a dedicated state submodule if it improves readability.
2. Add `model_state` to copy the initialized atmospheric state and compatibility metadata. Return an error for an uninitialized source.
3. Add `initialize_from_state` for both fresh and initialized destinations. Validate the state and configuration before modifying the destination or destroying an existing stepper.
4. Allocate a local internal `AtmosphereState`, copy the primary fields, compute gravity, and use `finalize_atmosphere_state` to regenerate temperature-dependent and particle-optical properties.
5. Commit the candidate and prepare the atmosphere with `KeepCurrentProfile`. Reuse existing state-copy and rollback machinery rather than duplicating array assignments throughout the API.
6. Define failure handling explicitly: validation/candidate failures preserve the old atmosphere and stepper; failures after commit or stepper teardown must report whether the old state remains usable. Do not label the entire operation failure atomic unless every failure path supports that promise.

Relevant existing code:

- `src/photochem/photochem_evoatmosphere.f90`: internal `AtmosphereState` and public interfaces.
- `src/photochem/photochem_evoatmosphere_init.f90`: allocation, `copy_model_to_state`, `copy_state_to_model`, initialization, and `finalize_atmosphere_state`.
- `src/photochem/photochem_evoatmosphere_grid.f90`: gravity/grid construction and recovery-status handling.
- `src/photochem/photochem_evoatmosphere_rhs.f90`: atmosphere preparation and profile synchronization policies.
- `src/photochem/photochem_vars.f90`: stored profile and TOA-maintenance types.

Keep interface documentation in `photochem_evoatmosphere.f90`, following the project's existing convention. Attach methods to the model only when public or needed across implementation files.

### Pass 2: HDF5 read/write

1. Implement Fortran `save_to_h5` and HDF5 construction/loading on the state type using the project's existing HDF5 infrastructure. The Python layer delegates to these routines, so both languages use one file format.
2. Use an explicit integer schema version and readable groups such as `/compatibility`, `/atmosphere`, `/press_temp_edd_profile`, and `/toa_pressure_maintenance`.
3. Store primary arrays as double precision without decimal rounding. Specify dataset dimensions, species string representation, units, and Fortran/Python array orientation in the format documentation.
4. Handle zero-particle states and disabled/unallocated profile arrays explicitly. Persist enablement flags even when their associated arrays are absent.
5. Load into a local candidate state and validate schema, required datasets, types, dimensions, and values before replacing an existing state object.
6. Prefer uncompressed datasets initially to avoid filter/plugin requirements. Compression can be added later if there is a demonstrated need.
7. Close HDF5 resources on every error path. Write to a temporary sibling file and replace the destination only after a successful close, so failed writes do not destroy an existing saved state. Decide and document overwrite behavior before implementing it.

HDF5 does not contain the mechanism, stellar flux, cross sections, callbacks, or executable solver state. It must not require a live model merely to load an `EvoAtmosphereState` object.

### Pass 3: Python bindings

1. Add C API ownership functions for state allocation/loading, deletion, snapshot creation, restoration, and saving.
2. Wrap the owned Fortran object as `EvoAtmosphereState` and export it from `photochem`.
3. Expose the intended constructor and methods: `EvoAtmosphereState(filename)`, `state.save_to_h5(filename)`, `pc.model_state()`, and `pc.initialize_from_state(state)`.
4. Ensure the state remains valid after the source model is destroyed or changed. Free the owned object exactly once and translate Fortran errors into the established Python exception type.
5. Begin with an opaque Python state object unless there is a concrete need to edit its fields. Public editable array views would require additional validation and ownership rules.

### Pass 4: Verification and documentation

Extend the existing Fortran API and Python suites with meaningful round-trip cases:

- Capture and restore an atmosphere into a fresh destination.
- Restore into an already initialized destination and verify stepper invalidation on success.
- Confirm source/state/destination independence after the source changes or is destroyed.
- Cover particle-free and particle-bearing states, including nondefault particle-radius profiles.
- Cover disabled, continuous, and periodic profiles and TOA-maintenance configurations.
- For periodic mode, deliberately save a valid profile mismatch and verify restore does not resynchronize it.
- Compare in-memory and HDF5 restore results, including regenerated optical and thermodynamic fields.
- Verify destination general solver settings, including `max_resync_failures`, are retained.
- Reject incompatible species/layer counts, malformed datasets, unsupported schemas, and invalid geometry before disturbing an existing model.
- Check applicable failure paths for rollback, state usability, HDF5 handle cleanup, and preservation of an existing file after a failed save.

Use the `photochem` conda environment. Build the Fortran library and Python extensions, run the relevant existing suites, and document the API and file-format scope. Run a short integration after restore to check that regenerated state is usable; bit-for-bit reproduction of the original CVODE trajectory is not an acceptance criterion.

## Deferred work

- CVODE restart/history serialization and preservation of convergence counters or logical time.
- Full model configuration, custom boundary conditions, callback serialization, and reaction/optical data packaging.
- Restore with a different mechanism or layer count.
- Automatic conversion of earlier HDF5 schemas.
- Replacing the gas-giant dictionary API. Once the base state API is stable, the extension could use it for its atmospheric component while saving its climate-grid metadata separately.

## Completion criteria

The example API works for a fresh compatible destination, in-memory and HDF5 restores agree, omitted derived fields are regenerated correctly, periodic-profile mismatch is preserved, general destination solver settings are retained, and failure semantics are documented and tested. None of this implementation is required for the current commit.
