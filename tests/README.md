# Test suite

The tests are divided by purpose so that correctness checks, broad memory
exercises, and longer integrations can evolve independently. Some overlap is
intentional: focused API tests also run under Valgrind because they exercise
allocation-heavy code paths that should retain memory-safety coverage.

CI uses one rolling, recent conda Python and GNU compiler toolchain and records
the resolved versions in its log. The full Fortran correctness suite and
installed Python smoke tests run normally. Valgrind covers both Clima tests,
the Equilibrate test, and Photochem's focused API and broad memory tests.

CI deliberately uses bounded solver exercises rather than requiring complete
photochemical steady states or a converged radiative-convective-equilibrium
calculation. The focused suites take real rocky-planet and gas-giant steps, and
the Clima adiabat test performs a partial RCE solve. Longer end-to-end models
remain appropriate for documentation examples and release validation when the
documentation is built.

Test inputs are resolved from the source tree, so the Fortran executables and
Python smoke test may be launched from any working directory.
Generated Fortran test files are written under the component build directory
and removed by the test that creates them; Python tests use temporary
directories managed by the standard library.

## Test programs

### Clima

- `test_adiabat` exercises atmospheric-profile construction, surface-temperature
  solves, column inventories, radiative-convective equilibrium, and custom
  particle and mixing-ratio inputs.
- `test_radtran` exercises standalone radiative transfer, atmosphere-file
  interpolation, opacity reporting, and custom optical properties.
- `test_python.py` checks the installed Clima wrapper's public state, nested
  radiative-transfer views, solver controls, and rebin utilities.
- The experimental Clima `Climate` model is intentionally excluded from this
  repository, so its corresponding upstream test is not imported.

### Equilibrate

- `test_equilibrate` checks YAML and legacy thermodynamic input, equilibrium
  compositions, metallicity solves, the Sonora species set, repeated
  construction, and NASA7 thermodynamic functions.
- `test_python.py` checks the installed `photochem.equilibrate` wrapper against
  a known equilibrium composition, including public mixture properties,
  control round trips, and copy semantics.

### Photochem

- `test_input` contains focused input parsing and validation tests.
- `test_api` contains focused correctness and error-behavior tests for public
  `EvoAtmosphere` operations.
- `test_jacobian` compares the analytical and automatic-differentiation
  chemistry Jacobians over focused numerical cases.
- `test_production_loss` characterizes the species production-and-loss
  diagnostic and checks reactions, rainout, condensation and evaporation,
  vertical transport, boundary fluxes, distributed sources, custom rates, and
  Zahnle escape. For evolved species, it verifies that all reported
  contributions reconstruct the full right-hand side.
- `test_memory` broadly exercises `EvoAtmosphere` workflows. It checks returned
  errors, but does not comprehensively validate numerical outputs; its primary
  purpose is execution under Valgrind.
- `test_evolution` is a longer end-to-end evolution exercise. It is built in CI
  to prevent source drift, but is intended for manual runs and is not executed
  in CI.
- `test_python.py` is a smoke test for the installed Python wrappers.

The analytical chemistry Jacobian is the production default. Automatic
differentiation remains selectable as an independent correctness oracle. Any
new composition-dependent chemistry tendency must receive a matching
analytical derivative and analytical-versus-autodiff coverage in
`test_jacobian`; unsupported terms must not be silently omitted from the
analytical path.

## Production-and-loss diagnostic contract

`production_and_loss` is a user-facing explanation of the tendency of one
selected species. Every reported profile has units of molecules/cm^3/s, is
nonnegative, and is accompanied by a reaction equation or process label.
Column-integrated values have units of molecules/cm^2/s. A signed contribution
is represented by its positive part in `production` and the magnitude of its
negative part in `loss`.

The diagnostic accounts for chemical reactions, rainout,
condensation and evaporation, internal vertical transport, lower and upper
boundary fluxes, distributed sources, custom rates, and Zahnle hydrogen
escape. Internal transport excludes exchange across the model boundaries, so
its vertically integrated net contribution must vanish apart from roundoff.
Lower and upper boundary exchange is reported simply as `lower boundary flux`
and `upper boundary flux`, independent of the configured boundary-condition
parameterization. Fixed-density and fixed-pressure lower boundaries receive an
implicit reservoir contribution that makes the reported bottom-cell tendency
match the replaced zero right-hand-side equation.

For an evolved species, the sum of all production columns minus the sum of all
loss columns must reconstruct its full right-hand-side tendency. For a
short-lived species, the result instead describes its algebraic chemical
balance because it has no transport ODE equation. Reaction entries with a
stoichiometric multiplicity greater than one are consolidated into one labeled
contribution with the appropriate multiplier.

## Running locally

After configuring and building the project:

```sh
./build/tests/clima/test_adiabat
./build/tests/clima/test_radtran
./build/tests/equilibrate/test_equilibrate
./build/tests/photochem/test_input
./build/tests/photochem/test_api
./build/tests/photochem/test_jacobian
./build/tests/photochem/test_production_loss
./build/tests/photochem/test_memory
```

The longer evolution test can be run manually when needed:

```sh
./build/tests/photochem/test_evolution
```

On a system with Valgrind, check all three components. The Photochem commands
cover both the broad memory harness and the focused API tests:

```sh
valgrind --error-exitcode=1 --leak-check=full ./build/tests/clima/test_adiabat
valgrind --error-exitcode=1 --leak-check=full ./build/tests/clima/test_radtran
valgrind --error-exitcode=1 --leak-check=full ./build/tests/equilibrate/test_equilibrate
valgrind --error-exitcode=1 --leak-check=full ./build/tests/photochem/test_memory
valgrind --error-exitcode=1 --leak-check=full ./build/tests/photochem/test_api
```

After installing the Python package, run its smoke test:

```sh
python tests/equilibrate/test_python.py
python tests/clima/test_python.py
python tests/photochem/test_python.py
```
