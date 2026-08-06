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
- `benchmark_block_thomas` compares the standard band solver with the
  experimental block-Thomas solver on a perturbed Modern Earth atmosphere. It
  is a manual performance benchmark, not a CI correctness test.
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

Run both linear solvers to a fixed physical time with:

```sh
./tests/benchmark_block_thomas
```

The default final time is `1.0e10` seconds. Supply another positive value to
change the benchmark duration, for example:

```sh
./tests/benchmark_block_thomas 1.0e11
```

The benchmark doubles the Modern Earth lower-boundary CH4 flux, runs both
solvers from the same initial atmosphere, and reports wall time, CVODE work
counters, and final-state differences.

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
