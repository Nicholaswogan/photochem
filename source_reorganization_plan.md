# Source Reorganization Plan

> **Temporary development document.** This file records the design and
> implementation checklist for roadmap item 8.1. It may be committed and
> updated while the work is in progress, then deleted once the reorganization
> is complete and the lasting architectural decisions are reflected in the
> source tree and permanent documentation.

## Status

- [ ] Phase 8.1a: reorganize the existing Photochem source, tests, bindings,
      and Python package boundary
- [ ] Phase 8.1b: bring Clima source into this repository
- [ ] Phase 8.1c: bring Equilibrate source into this repository
- Optional after completion: extract demonstrated shared code into neutral
  modules if it is clearly worthwhile
- [ ] Remove this temporary plan and close roadmap item 8.1

## Objective

Organize the Photochem implementation, tests, bindings, and package boundaries
around explicit ownership before adding the Clima and Equilibrate sources to
this repository. The final repository should make it easy to answer:

1. Which module constructs and maintains each major model object?
2. Which routines coordinate a complete model lifecycle?
3. Which routines are standalone calculations that can be reused?
4. Which functionality might optionally be shared by multiple models after
   their existing implementations are safely in one repository?

The organizing rules are:

```text
Owner modules
  Construct and maintain one owned type

EvoAtmosphere modules
  Coordinate dat + var + wrk
  Manage lifecycle and transactions
  Implement public model methods

Standalone modules
  Pure equations, numerical kernels, constants, and enums
```

The reorganization should improve ownership without changing numerical or
public API behavior unless a change is explicitly identified, tested, and
committed as a separate refactor.

## Architectural invariants

### PhotochemData

`PhotochemData` is the immutable, resolved model definition. Once constructed,
its contents should not change during atmospheric initialization or model
integration.

Its constructor owns:

- the chemical mechanism;
- species, atoms, reactions, and reaction indices;
- thermodynamic and Henry-law data;
- photolysis, particle-optical, and Rayleigh input data;
- immutable planet and model choices resolved from settings;
- immutable array and Jacobian dimensions.

Construction helpers for `PhotochemData` must not use `PhotochemVars` as input
or scratch storage. Values such as `data_dir` should be passed explicitly.

### PhotochemVars

`PhotochemVars` is persistent configuration and prepared model state that can
change between integrations without rereading the mechanism and static data.

A constructed `PhotochemVars` is valid as **configured model storage**, but it
does not yet describe a physically initialized atmosphere. Atmosphere-dependent
arrays may be allocated and left unset during construction. They must not be
interpreted until the owning `EvoAtmosphere` has
`atmosphere_initialized == .true.`.

Its constructor owns:

- boundary conditions;
- condensation and rainout configuration;
- stellar photon flux and radiative controls;
- integration and convergence settings;
- allocation of arrays whose shapes depend on `PhotochemData` and `nz`.

It also owns internal maintenance routines for its prepared, persistent fields,
including temperature-dependent cross sections and Gibbs energies. These
routines should continue accepting explicit arrays where needed so they can
operate on transactional candidate state.

### PhotochemWrk

`PhotochemWrk` owns scratch arrays, current integration state, convergence
history, SUNDIALS resources, and other runtime bookkeeping. The current
`PhotochemWrk` base type and `PhotochemWrkEvo` extension should be consolidated
into this single concrete type; the extension no longer represents a useful
architectural distinction.

Construction should use an ordinary constructor with explicit integer
dimensions rather than depending on complete `PhotochemData` or
`PhotochemVars` objects:

```fortran
wrk = PhotochemWrk(nsp, np, nq, nz, nrT, kj, nw)
```

`wrk%usol` and other prepared runtime quantities become meaningful only after
successful atmospheric preparation.

### EvoAtmosphere

`EvoAtmosphere` is the aggregate root and lifecycle authority. It owns
`PhotochemData`, `PhotochemVars`, and `PhotochemWrk`, and it alone determines
whether those components collectively describe an initialized atmosphere.

It is responsible for:

- sequencing component construction;
- atmosphere initialization;
- candidate-state validation and transactional commits;
- rollback after failed preparation;
- public state mutations;
- RHS and Jacobian evaluation;
- integration and stepper lifecycle;
- model-level diagnostics and output.

`AtmosphereState` remains an internal transaction type used to construct and
validate atmospheric changes without copying all of `PhotochemVars`.

## Construction dependency graph

The desired construction flow is:

```text
PhotoSettings
    |
    +--------------------------+
    |                          |
    v                          v
PhotochemData ------------> PhotochemVars

explicit dimensions ------> PhotochemWrk

PhotochemData + PhotochemVars + PhotochemWrk
                         |
                         v
                  EvoAtmosphere
```

In approximate Fortran:

```fortran
settings = PhotoSettings(settings_file, err)
if (allocated(err)) return

dat = PhotochemData(mechanism_file, settings, data_dir, err)
if (allocated(err)) return

var = PhotochemVars(dat, settings, flux_file, err)
if (allocated(err)) return

wrk = PhotochemWrk(dat%nsp, dat%np, dat%nq, var%nz, &
                   dat%nrT, dat%kj, dat%nw)
```

`PhotochemVars` may depend on both `PhotoSettings` and the completed
`PhotochemData`. Mutable settings should not be copied through
`PhotochemData` merely to make the diagram strictly linear.

## Proposed source tree

Keep one public facade for each eventual model at the top of `src/`, with a
flat implementation directory for each model:

```text
src/
├── CMakeLists.txt
├── photochem.f90
├── clima.f90                         # Added in Phase 8.1b
├── equilibrate.f90                   # Added in Phase 8.1c
├── photochem_version.f90.in
├── dependencies/
│
├── photochem/
│   ├── CMakeLists.txt
│   ├── photochem_const.f90
│   ├── photochem_enum.f90
│   ├── photochem_eqns.f90
│   ├── photochem_radtran.f90
│   ├── photochem_evoatmosphere_chemistry.f90
│   ├── photochem_settings.f90
│   ├── photochem_data.f90
│   ├── photochem_vars.f90
│   ├── photochem_wrk.f90
│   ├── photochem_evoatmosphere.f90
│   ├── photochem_evoatmosphere_init.f90
│   ├── photochem_evoatmosphere_profiles.f90
│   ├── photochem_evoatmosphere_grid.f90
│   ├── photochem_evoatmosphere_boundary.f90
│   ├── photochem_evoatmosphere_rhs.f90
│   ├── photochem_evoatmosphere_integrate.f90
│   └── photochem_evoatmosphere_output.f90
│
├── clima/                            # Flat Clima implementation
└── equilibrate/                      # Flat Equilibrate implementation
```

An optional `src/shared/` directory may be introduced later if obvious reuse
is worth extracting. It is not part of the required initial layout or the
definition of completion for roadmap item 8.1.

The names under `photochem/` retain the `photochem_` prefix so module names and
filenames remain easy to associate in compiler output, build files, and code
searches.

## Bindings and Python package layout

Before importing Clima and Equilibrate, establish the final component-first
bindings layout using the existing Photochem wrappers:

```text
bindings/
├── CMakeLists.txt
├── photochem/
│   ├── CMakeLists.txt
│   ├── c_api/
│   └── cython/
├── clima/                            # Added with the Clima import
│   ├── c_api/
│   └── cython/
└── equilibrate/                      # Added with the Equilibrate import
    ├── c_api/
    └── cython/
```

Each component owns its complete C API and Cython stack. Shared binding code
may be extracted later only if concrete reuse appears. Compiled extensions
must retain their existing private Python names: `photochem._photochem`,
`photochem._clima`, and `photochem._equilibrate`.

After the Photochem wrappers move, the repository-root `photochem/` directory
is only the installable Python package. Its current public organization is a
reasonable initial target:

```text
photochem/
├── __init__.py
├── clima.py
├── equilibrate.py
├── io.py
├── extensions/
└── utils/
```

Audit this package boundary before importing the other components, but avoid a
broad Python API redesign. The thin Clima and Equilibrate facades, extension
modules, utility imports, and existing public import paths should remain stable
during the source import.

## Test layout

Establish component test ownership before importing external tests:

```text
tests/
├── CMakeLists.txt
├── README.md
├── photochem/
│   ├── CMakeLists.txt
│   ├── test_*.f90
│   ├── test_python.py
│   └── test fixtures
├── clima/                            # Added with the Clima import
├── equilibrate/                      # Added with the Equilibrate import
└── integration/                      # Cross-component behavior only
```

Keep each component directory flat initially. Use component-prefixed CMake
target names so targets remain globally unique. Test fixtures must be resolved
independently of the build directory; do not preserve assumptions that the
build tree is a sibling of the source tree. Move existing tests mechanically
before splitting truly cross-component coverage into `tests/integration/`.
Continue invoking the test executables and scripts directly; introducing CTest
is a separate optional improvement rather than part of this reorganization.

## Target Photochem modules

### `photochem_settings.f90`

Owns:

- `PhotoSettings`;
- `SettingsBC`;
- `SettingsParticle`;
- `CondensationParameters`;
- the `PhotoSettings` constructor interface;
- settings YAML parsing and validation;
- private string, list, and duplicate-checking helpers used only by settings.

This absorbs the settings portions of:

- `photochem_types.f90`;
- `photochem_types_create.f90`.

### `photochem_data.f90`

Owns the following types:

- `PhotochemData`;
- `Reaction`;
- `BaseRate` and all derived rate types;
- `Efficiencies`;
- `MultiArrheniusRate`;
- `ProdLoss`;
- `ThermodynamicData`;
- `XsectionData`;
- `ParticleXsections`.

Provides a constructor resembling:

```fortran
dat = PhotochemData(mechanism_file, settings, data_dir, err)
```

Owns the current static-data routines for:

- parsing mechanisms, species, particles, atoms, and reactions;
- parsing reaction equations and rate parameters;
- resolving species and production/loss indices;
- validating duplicate reactions and short-lived species constraints;
- reading thermodynamic data;
- reading Henry-law data;
- reading wavelength bins;
- reading photolysis cross sections;
- reading particle optical data;
- reading Rayleigh data;
- computing immutable Jacobian storage dimensions.

The current `var` arguments used only to obtain `var%data_dir` should be
replaced by an explicit `data_dir` argument. Other unnecessary `var` arguments
should be removed.

`parse_reaction` belongs here because it is part of mechanism interpretation,
even if it remains public for focused Fortran tests.

### `photochem_vars.f90`

Owns the following types and interfaces:

- `PhotochemVars`;
- `PressureTempEddProfile`;
- `TOAPressureMaintenance`;
- `binary_diffusion_fcn`;
- `time_dependent_flux_fcn`;
- `time_dependent_rate_fcn`;
- `time_dependent_rate_fcns`.

Provides a constructor resembling:

```fortran
var = PhotochemVars(dat, settings, flux_file, err)
```

Owns routines for:

- applying mutable settings;
- allocating boundary-condition and model-grid arrays;
- resolving and validating boundary conditions against `PhotochemData`;
- initializing particle condensation parameters;
- reading the stellar photon flux;
- refreshing temperature-dependent prepared variables;
- computing Gibbs energies on a supplied temperature grid;
- interpolating gas photolysis data onto a supplied state;
- interpolating particle optical data onto a supplied state.

The constructor establishes configured storage, but does not set
`EvoAtmosphere%atmosphere_initialized`.

### `photochem_wrk.f90`

Owns:

- `PhotochemWrk`;
- `SundialsData`;
- `SundialsDataFinalizer`;
- the `PhotochemWrk` constructor and work-array allocation;
- SUNDIALS cleanup and finalization;
- runtime bookkeeping reset operations.

It should depend only on low-level constants, intrinsic modules, and SUNDIALS,
not on complete `PhotochemData` or `PhotochemVars` objects.

The current base type and Evo-specific extension should be merged. All fields
currently in `PhotochemWrkEvo`, including surface pressure, hydrostatic arrays,
mixing ratios, and robust-stepper bookkeeping, become fields of
`PhotochemWrk`. Construction should be:

```fortran
wrk = PhotochemWrk(nsp, np, nq, nz, nrT, kj, nw)
```

### `photochem_evoatmosphere.f90`

Defines the aggregate `EvoAtmosphere` type and its public API. It contains:

- `EvoAtmosphere`;
- generic constructor declarations;
- type-bound procedure declarations;
- private transaction types such as `AtmosphereState`;
- public result types such as `ProductionLoss`;
- lifecycle state and its documented invariants.

`AtmosphereState` is defined privately in this parent module so every
EvoAtmosphere submodule can use it through host association.

Implementation remains in focused submodules.

### `photochem_evoatmosphere_init.f90`

Owns aggregate and atmosphere construction:

- `create_EvoAtmosphere`;
- `create_EvoAtmosphere_static`;
- `initialize_from_atmosphere_file`;
- `initialize_atmosphere_z`;
- `initialize_atmosphere_p`;
- `copy_model_to_state`;
- `copy_state_to_model`;
- rollback helpers;
- `require_atmosphere_initialized`;
- atmosphere-file reading;
- altitude- and pressure-profile mapping;
- profile interpolation helpers;
- `finalize_atmosphere_state`.

This absorbs the atmosphere-specific functionality currently in
`photochem_input_after_read.f90`. The unused nontransactional
`finalize_atmosphere_initialization(dat, var, err)` should be removed if a final
usage audit confirms that it has no callers. Mapping, file-reading, and
interpolation routines are private helpers; `finalize_atmosphere_state`
remains a private parent-module procedure because grid operations also use it.

### `photochem_evoatmosphere_profiles.f90`

Owns public and internal profile mutations:

- `set_temperature`;
- `set_press_temp_edd`;
- `map_press_temp_edd`;
- `set_press_temp_edd_profile`;
- `clear_press_temp_edd_profile`;
- `reset_press_temp_edd_profile`;
- `apply_press_temp_edd_profile`;
- local transaction types such as `TemperatureState` and
  `PressTempEddState`.

### `photochem_evoatmosphere_grid.f90`

Owns vertical-grid transactions:

- `update_vertical_grid`;
- `TOA_at_pressure`;
- `VerticalGridWork`;
- altitude and pressure residual functions;
- candidate-grid construction;
- candidate temperature/Kzz seeding and validation;
- candidate composition and particle-radius mapping;
- grid commit and rollback logic.

### `photochem_evoatmosphere_boundary.f90`

Owns:

- `set_lower_bc`;
- `set_upper_bc`;
- `set_rate_fcn`;
- `gas_fluxes`;
- model-level boundary-condition validation.

### `photochem_evoatmosphere_rhs.f90`

Continues to own:

- preparation of the current atmosphere for RHS evaluation;
- reaction, diffusion, condensation, and boundary RHS terms;
- chemistry-only and full Jacobians;
- production/loss diagnostics;
- RHS-specific lower-boundary application.

### `photochem_evoatmosphere_integrate.f90`

Continues to own:

- CVODE callbacks;
- stepper creation, configuration, restart, and destruction;
- `evolve`;
- `step`;
- robust integration;
- convergence checking;
- model-top pressure maintenance orchestration;
- steady-state search.

Restart-file helpers may remain here initially when they are tightly coupled to
the integration workflow.

### `photochem_evoatmosphere_output.f90`

Owns:

- `out2atmosphere_txt`;
- `out2atmosphere_txt_base`;
- output-only formatting and file helpers.

### `photochem_evoatmosphere_chemistry.f90`

Owns the chemistry and atmospheric-process kernels used privately by the
EvoAtmosphere RHS workflow. They retain explicit `PhotochemData` and
`PhotochemVars` interfaces rather than taking the complete model object:

- `reaction_rates`;
- `photorates`;
- `chempl`;
- `chempl_sl`;
- `chempl_t`;
- `rainout`;
- `gas_saturation_density`;
- `molec_per_particle`.

`out2atmosphere_txt_base` should move out of this module because it is not a
chemistry kernel.

### Existing standalone modules

Retain their current focused purposes:

- `photochem_const.f90`: numerical and physical constants;
- `photochem_enum.f90`: enum-like integer definitions;
- `photochem_eqns.f90`: pure or nearly pure primitive equations;
- `photochem_radtran.f90`: the self-contained UV two-stream solver.

`photochem_eqns` should remain independent of model-owned derived types.
Type-specific operations such as thermodynamic-data evaluation belong with
their owning type; primitive numerical and physical equations belong in
`photochem_eqns`.

## Expected removal of current files

Once their contents have moved and callers have been updated, remove:

```text
src/photochem_types.f90
src/photochem_types_create.f90
src/photochem_common.f90
src/input/photochem_input.f90
src/input/photochem_input_read.f90
src/input/photochem_input_after_read.f90
src/evoatmosphere/photochem_evoatmosphere_utils.f90
```

The `src/input/` and `src/evoatmosphere/` directories should then disappear.
EvoAtmosphere remains a meaningful module grouping expressed through filenames
and Fortran submodules rather than a unique nested directory.

## Import strategy for Clima and Equilibrate

After the Photochem source, bindings, and Python package boundaries have been
established, copy each external component's source and bindings directly into
its final locations with as little modification as possible:

```text
src/clima/
bindings/clima/
src/equilibrate/
bindings/equilibrate/
```

Duplicate code is acceptable at this stage. The initial objective is to make
this repository own and build the three existing implementations, not to
redesign all three simultaneously. Apart from flattening their source layouts,
updating build integration, and making unavoidable path or module-order
adjustments, preserve the imported code as-is.

The three complete models should remain independently owned aggregates. Do not
attempt to create universal `Data`, `Vars`, `Wrk`, or model-state types shared
by all three codes.

The import must also account for these existing repository assumptions:

- Record the imported Clima and Equilibrate versions and source revisions for
  provenance, and retain their GPL license provenance. Preserve their version
  modules during the mechanical import, then remove the independent component
  release versions once both components are owned by this repository.
- Centralize external dependency targets in this repository. Clima currently
  creates `h5fortran`, `finterp`, `dop853`, and `minpack` targets, while
  Photochem already consumes `finterp` and `minpack` indirectly. Preserve
  LAPACK discovery and the shared OpenMP policy without defining duplicate
  CMake targets.
- Preserve the required build direction: Clima provides
  `clima_saturationdata` and `clima_useful` modules used by Photochem.
  Equilibrate is consumed through its Python extension and higher-level
  Photochem utilities rather than by the Photochem Fortran library.
- Preserve wrapper dependency declarations, generated FYPP sources, extension
  install destinations, and the
  `photochem_clima_data` default-data behavior.
- Import test fixtures required by the component suites, including the Clima
  settings/species/atmosphere inputs and Equilibrate thermodynamic data. Keep
  those assets under their component test directories when they are not
  runtime package data.
- Defer standalone examples, notebooks, and conversion/development scripts.
  A future top-level `docs/` tree will provide examples and documentation for
  all three codes after the combined package is stable.

## Optional shared implementation

A `src/shared/` directory is optional and is not required to complete roadmap
item 8.1. Consider it only after the in-tree Photochem, Clima, and Equilibrate
implementations build and test successfully in their initially imported form.

If the combined source reveals an especially clear opportunity, small neutral
concepts or computational kernels may be extracted:

```text
                     shared
                   /    |    \
                  v     v     v
          photochem   clima   equilibrate
```

Possible candidates include:

- numerical kinds and truly universal physical constants;
- saturation-vapor-pressure data and evaluation;
- thermodynamic polynomial data and evaluation;
- chemical composition and elemental bookkeeping;
- hydrostatic or equation-of-state primitives;
- narrowly reusable radiative-transfer numerical kernels.

For example, Photochem currently imports `SaturationData` from Clima. It may
eventually be cleaner for both models to depend on a neutral saturation module,
but leaving the initial imported ownership unchanged is acceptable.

Do not create `shared/` merely to make the final tree appear deduplicated. Any
shared module that is introduced should have:

- a domain-neutral name;
- at least two concrete consumers;
- no dependency on a complete model aggregate;
- a focused and testable API;
- enough benefit to justify modifying otherwise working imported code.

Existing external libraries such as `futils` and `fortran-yaml-c` should
continue providing generic interpolation, file, and YAML facilities rather than
being duplicated internally.

## Implementation checklist

Each pass should compile and run the relevant tests before being committed.
Mechanical moves and behavioral changes should be separate whenever practical.

### Pass 1: Settings ownership

- [x] Create `photochem_settings.f90`.
- [x] Move settings types out of `photochem_types.f90`.
- [x] Move the `PhotoSettings` constructor and parsing implementation out of
      `photochem_types_create.f90`.
- [x] Update module imports without changing behavior.
- [x] Add or retain focused settings-construction tests.
- [x] Build and run the full Fortran test suite.

### Pass 2: PhotochemData ownership

- [x] Create `photochem_data.f90`.
- [x] Move `PhotochemData` and its subordinate data types.
- [x] Move mechanism, reaction, thermodynamic, Henry, wavelength, cross-section,
      particle-optical, and Rayleigh readers into the owner module.
- [x] Replace construction-time `var%data_dir` access with explicit
      `data_dir` arguments.
- [x] Remove all unnecessary `PhotochemVars` arguments from data readers.
- [x] Split immutable settings application from mutable vars settings.
- [x] Add the isolated `PhotochemData` constructor.
- [x] Add focused constructor tests.
- [x] Build and run the full Fortran test suite.

### Pass 3: PhotochemVars ownership

- [x] Create `photochem_vars.f90`.
- [x] Move `PhotochemVars` and its subordinate types and callback interfaces.
- [x] Move mutable settings and boundary-condition setup into this module.
- [x] Move stellar-flux reading into this module.
- [x] Move model-grid allocation into this module.
- [x] Move temperature-dependent prepared-variable maintenance into this
      module while preserving explicit-array interfaces.
- [x] Add the isolated `PhotochemVars` constructor.
- [x] Document configured validity versus initialized-atmosphere validity.
- [x] Add focused constructor tests.
- [x] Build and run the full Fortran test suite.

### Pass 4: Work ownership

- [x] Create `photochem_wrk.f90`.
- [x] Merge the current `PhotochemWrk` and `PhotochemWrkEvo` into one
      `PhotochemWrk` type.
- [x] Move `PhotochemWrk` and the SUNDIALS owner types into the new module.
- [x] Replace `%init` allocation with the explicit-dimension constructor
      `wrk = PhotochemWrk(...)`.
- [x] Update `EvoAtmosphere` and all call sites to use the consolidated type.
- [x] Preserve reliable SUNDIALS finalization and cleanup semantics.
- [x] Update imports and tests.
- [x] Build and run the full Fortran test suite.

### Pass 5: Aggregate static construction

- [x] Replace `setup_static` and `read_static_files` orchestration with the new
      owner constructors.
- [x] Make `create_EvoAtmosphere_static` visibly compose settings, data, vars,
      and work objects.
- [x] Confirm that successful static construction leaves
      `atmosphere_initialized == .false.`.
- [x] Confirm that all allocated-but-uninitialized atmosphere fields remain
      inaccessible through guarded public operations.
- [x] Remove obsolete joint-mutator routines.
- [x] Build and run the full Fortran and Python test suites.

### Pass 6: Atmosphere initialization ownership

- [x] Move atmosphere-file reading and profile mapping into
      `photochem_evoatmosphere_init.f90`.
- [x] Define `AtmosphereState` privately in the parent EvoAtmosphere module.
- [x] Keep mapping, file-reading, and interpolation routines private to the
      initialization implementation.
- [x] Keep `finalize_atmosphere_state` available privately to EvoAtmosphere
      submodules.
- [x] Preserve transactional initialization and rollback behavior.
- [x] Remove `finalize_atmosphere_initialization` if the final caller audit
      confirms that it is unused.
- [x] Import `parse_reaction` tests directly from `photochem_data`.
- [x] Remove the obsolete `photochem_input` module and source directory.
- [x] Build and run the full Fortran and Python test suites.

### Pass 7: Split EvoAtmosphere operations by responsibility

- [x] Create the profiles submodule and move profile operations.
- [x] Create the grid submodule and move vertical-grid operations.
- [x] Create the boundary submodule and move boundary operations.
- [x] Create the output submodule and move output operations.
- [x] Keep RHS and integration in their existing focused submodules.
- [x] Remove `photochem_evoatmosphere_utils.f90`.
- [x] Build and run the full Fortran and Python test suites.

### Pass 8A: Standalone and compatibility cleanup

- [x] Move `ProductionLoss` into the EvoAtmosphere parent module and export it
      through the `photochem` facade.
- [x] Update wrappers and remove `photochem_types.f90`.
- [x] Rename `photochem_common.f90` to `photochem_evoatmosphere_chemistry.f90`.
- [x] Move `out2atmosphere_txt_base` to the output submodule.
- [x] Build and run the full Fortran and Python test suites.

### Pass 8B: Flat Photochem layout

- [x] Move all Photochem implementation files into the flat
      `src/photochem/` directory.
- [x] Update CMake source lists, FYPP inputs, module output paths, and install
      rules.
- [x] Remove the empty `src/evoatmosphere/` directory.
- [x] Perform a clean build in the `photochem` conda environment.
- [x] Run all Fortran and Python tests locally.
- [ ] Confirm CI on supported compilers and Python versions after pushing the
      branch.

### Pass 8C: Component-oriented test layout

- [x] Move existing tests and fixtures into `tests/photochem/`.
- [x] Make the root test CMake file coordinate component test directories.
- [x] Use component-prefixed CMake target names.
- [x] Remove build-location assumptions from test fixture paths.
- [x] Preserve current test behavior and run the full Fortran and Python suites.

### Pass 8D: Photochem bindings layout

- [x] Create `bindings/photochem/c_api/` and `bindings/photochem/cython/`.
- [x] Move the existing Photochem Fortran/C and Cython wrappers out of the
      Python package without changing wrapper behavior.
- [x] Move component binding CMake ownership under `bindings/photochem/` and
      connect it through `bindings/CMakeLists.txt`.
- [x] Preserve the installed extension name `photochem._photochem`.
- [x] Build the wheel and run the full Fortran and Python test suites.

### Pass 8E: Python package boundary audit

- [x] Confirm that the root `photochem/` directory contains only installable
      Python-package code after the bindings move.
- [x] Preserve `photochem.clima`, `photochem.equilibrate`, `extensions`,
      `utils`, and all existing public import paths.
- [x] Update packaging manifests and source-distribution inputs as needed.
- [x] Avoid broad Python module moves until all three components are in-tree.
- [x] Build and test the installed package.

### Pass 9A: Import Clima

- [x] Copy nearly all existing Clima source into a flat `src/clima/`
      implementation directory.
- [x] Record the imported version, source revision, and license provenance.
- [x] Integrate Clima's external dependencies without duplicate CMake targets,
      preserving LAPACK and OpenMP behavior.
- [x] Copy the Clima C API and Cython bindings into `bindings/clima/` while
      preserving `photochem._clima` and installed climate-data behavior.
- [x] Import the Clima tests and required fixtures into `tests/clima/` without
      importing examples or notebooks.
- [x] Build and test the in-tree Clima before removing its CPM dependency.

### Pass 9B: Import Equilibrate

- [x] Copy nearly all existing Equilibrate source into a flat
      `src/equilibrate/` implementation directory.
- [x] Record the imported version, source revision, and license provenance.
- [x] Copy the Equilibrate C API and Cython bindings into
      `bindings/equilibrate/` while preserving `photochem._equilibrate`.
- [x] Import the Equilibrate tests and required thermodynamic fixtures into
      `tests/equilibrate/` without importing examples or development scripts.
- [x] Build and test the in-tree Equilibrate before removing its CPM dependency.

### Pass 9C: Remove external component copies and stabilize

- [ ] Keep the imported implementations unchanged except for flattening,
      necessary build integration, and unavoidable path or module-order fixes.
- [ ] Accept duplicate code between the three model directories during the
      initial import.
- [ ] Preserve separate public facade modules.
- [ ] Update CMake so this repository builds all three components directly.
- [ ] Remove CPM/external-source copies once the in-tree builds are stable.
- [x] Remove the independent Clima and Equilibrate release versions; expose
      the unified package version only as `photochem.__version__`.
- [x] Verify the Clima and Equilibrate public facade exports after version
      consolidation.
- [ ] Verify runtime data defaults, wheel contents, and source distributions.
- [ ] Ensure macOS deployment tags produce installable Apple Silicon wheels.
- [ ] Preserve the existing Python package API during the transition.
- [ ] Stabilize clean local builds and CI.

### Optional follow-up: Extract demonstrated reuse

This work is explicitly optional and does not block completion of roadmap item
8.1:

- [ ] If worthwhile, inventory obvious cross-model duplication.
- [ ] If worthwhile, identify types or kernels with at least two real
      consumers.
- [ ] If worthwhile, move `SaturationData` to a neutral owner.
- [ ] If worthwhile, evaluate shared thermodynamic data and polynomial
      evaluation for Photochem and Equilibrate.
- [ ] Add focused tests for any shared module that is actually introduced.
- [ ] Ensure any shared dependency direction remains
      `shared -> model consumers`.

### Pass 10: Completion

- [ ] Audit module visibility and remove unintentionally public helpers.
- [ ] Audit for circular or unnecessary module dependencies.
- [ ] Verify installed Python modules for Photochem, Clima, and Equilibrate.
- [ ] Verify source distributions and conda builds do not depend on deleted
      external source trees.
- [ ] Confirm all three in-tree model implementations build and test even if
      no `shared/` directory was introduced.
- [ ] Update permanent developer documentation where necessary.
- [ ] Mark roadmap item 8.1 complete.
- [ ] Delete this temporary plan.

## Verification requirements

At the end of every material pass:

- configure and build using the `photochem` conda environment;
- run `make -j` from the configured `build/` directory;
- run relevant focused Fortran tests;
- run the complete Fortran test suite before committing an ownership change;
- install the Python package into the `photochem` environment when wrappers or
  source lists change;
- run the Python tests after public-module or packaging changes;
- inspect `git diff` to distinguish intended source moves from behavioral
  changes.

At the end of Phase 8.1a, the existing Photochem test results should match the
pre-reorganization baseline.

## Decisions to preserve

- Keep `src/photochem.f90` as the single public Photochem entry point.
- Prefer flat per-model source directories.
- Group subordinate types with their owning major type rather than creating a
  separate file for every small type.
- Treat `EvoAtmosphere` as the aggregate lifecycle authority.
- Treat constructed `PhotochemVars` as configured storage, not as proof that an
  atmosphere has been initialized.
- Consolidate `PhotochemWrk` and `PhotochemWrkEvo` into one `PhotochemWrk` and
  construct it with `wrk = PhotochemWrk(...)`.
- Keep one authoritative prepared temperature/Kzz state even when persistent
  profiles allow it to change during integration.
- Keep `AtmosphereState` as an internal transactional handoff rather than
  copying all of `PhotochemVars`.
- Initially copy Clima and Equilibrate into flat directories with minimal
  modification; duplicate code is acceptable.
- Treat `shared/` extraction as optional follow-up, not as a completion
  requirement.

## Open questions

Record decisions here as implementation reveals more information:

- [ ] Should Photochem's UV two-stream implementation remain model-specific or
      share numerical kernels with Clima radiative transfer?

These questions do not block Phase 8.1a unless a concrete dependency requires
an answer.
