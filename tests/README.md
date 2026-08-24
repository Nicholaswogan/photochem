# Test suite

The tests are divided by purpose so that correctness checks, broad memory
exercises, and longer integrations can evolve independently. Some overlap is
intentional: focused API tests also run under Valgrind because they exercise
allocation-heavy code paths that should retain memory-safety coverage.

The Fortran executables use paths relative to `build/` and must be run with
that directory as the working directory.

## Test programs

- `test_input` contains focused input parsing and validation tests.
- `test_api` contains focused correctness and error-behavior tests for public
  `EvoAtmosphere` operations.
- `test_jacobian` compares the analytical and automatic-differentiation
  chemistry Jacobians over focused numerical cases.
- `test_production_loss` characterizes the species production-and-loss
  diagnostic and checks that its currently reported reaction and rainout
  contributions reconstruct the corresponding chemistry tendency. As the
  diagnostic is completed, this test will expand to reconcile all reported
  contributions with the full right-hand side.
- `test_memory` broadly exercises `EvoAtmosphere` workflows. It checks returned
  errors, but does not comprehensively validate numerical outputs; its primary
  purpose is execution under Valgrind.
- `test_evolution` is a longer end-to-end evolution exercise. It is built in CI
  to prevent source drift, but is intended for manual runs and is not executed
  in CI.
- `test_python.py` is a smoke test for the installed Python wrappers and is run
  from `tests/`.

The analytical chemistry Jacobian is the production default. Automatic
differentiation remains selectable as an independent correctness oracle. Any
new composition-dependent chemistry tendency must receive a matching
analytical derivative and analytical-versus-autodiff coverage in
`test_jacobian`; unsupported terms must not be silently omitted from the
analytical path.

## Production-and-loss diagnostic contract

`production_and_loss` is a user-facing explanation of the tendency of one
selected species. Every reported
profile has units of molecules/cm^3/s, is nonnegative, and is accompanied by a
reaction equation or process label. Column-integrated values have units of
molecules/cm^2/s. A signed contribution is represented by its positive part in
`production` and the magnitude of its negative part in `loss`.

The completed diagnostic will account for chemical reactions, rainout,
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
stoichiometric multiplicity greater than one will ultimately be consolidated
into one labeled contribution with the appropriate multiplier.

## Running locally

After configuring and building the project:

```sh
cd build
./tests/test_input
./tests/test_api
./tests/test_jacobian
./tests/test_production_loss
./tests/test_memory
```

The longer evolution test can be run manually when needed:

```sh
./tests/test_evolution
```

On a system with Valgrind, check both the broad memory harness and the focused
API tests:

```sh
valgrind --error-exitcode=1 --leak-check=full ./tests/test_memory
valgrind --error-exitcode=1 --leak-check=full ./tests/test_api
```

After installing the Python package, run its smoke test from the repository's
`tests/` directory:

```sh
cd tests
python test_python.py
```
