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
- `test_memory` broadly exercises `EvoAtmosphere` workflows. It checks returned
  errors, but does not comprehensively validate numerical outputs; its primary
  purpose is execution under Valgrind.
- `test_evolution` is a longer end-to-end evolution exercise. It is built in CI
  to prevent source drift, but is intended for manual runs and is not executed
  in CI.
- `test_python.py` is a smoke test for the installed Python wrappers and is run
  from `tests/`.

## Running locally

After configuring and building the project:

```sh
cd build
./tests/test_input
./tests/test_api
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
