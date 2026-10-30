import argparse
import json
from copy import deepcopy
from types import SimpleNamespace
import time

import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory

from photochem.extensions import coupled

from photochem.utils import stars
from photochem.utils import zahnle_rx_and_thermo_files

TEST_DIR = Path(__file__).resolve().parent

def input_files(output_dir):
    output_dir = Path(output_dir)
    zahnle_rx_and_thermo_files(
        atoms_names=['H','N','O','C'],
        rxns_filename=str(output_dir / 'mechanism_coupled.yaml'),
        thermo_filename=None
    )
    stars.solar_spectrum(
        outputfile=str(output_dir / 'flux_coupled.txt'),
        stellar_flux=1361.0
    )

def initialize():
    with TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        input_files(output_dir)
        model = coupled.CoupledPhotochemistryClimate(
            mechanism_file=str(output_dir / 'mechanism_coupled.yaml'),
            settings_file=str(TEST_DIR / 'settings_coupled.yaml'),
            flux_file=str(output_dir / 'flux_coupled.txt')
        )
    return model


def test_smoke():
    model = initialize()
    assert model.solver is None
    assert not model.pc.atmosphere_initialized
    assert model.c.P_top == 1.0
    solver = coupled.FixedPointSolver(model)
    assert solver.niterations == 0
    assert not solver.converged and not solver.give_up
    assert solver.chemistry_converged is None
    assert solver.solver_history == []

    # These must fail before climate solving or atmospheric initialization.
    for options in ({'method': 'quasi_static'},
                    {'solver_options': {'relaxation': 0.0}},
                    {'solver_options': {'anderson_memory': -1}}, {}):
        try:
            model.initialize_solver(initial_guess={}, **options)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid solver controls were accepted.')
        assert model.solver is None
        assert not model.pc.atmosphere_initialized


def run_solver_checks():
    """Check the extracted solver locally without chemistry/climate solves."""
    model = coupled.CoupledPhotochemistryClimate.__new__(coupled.CoupledPhotochemistryClimate)
    model.solver = None
    initial = np.array([6.0, 300.0, 200.0])
    restored_states = []
    trials = []
    model.pc = SimpleNamespace(initialize_from_state=restored_states.append)
    state = {'photochem': {'toa_target_pressure': np.array(0.1)}}
    model._validate_state = lambda value: (value, np.ones(2), np.array([10., 1.]), {})
    model._fixed_point_state_from_climate = lambda climate: initial.copy()

    def fixed_point_map(x, require_chemistry_to_converge, *, chemistry_diagnostics):
        assert require_chemistry_to_converge
        trials.append(x.copy())
        chemistry_diagnostics.update(converged=True, steps=1)
        return x + np.array([0.2, 20.0, 10.0])

    model._fixed_point_map = fixed_point_map
    model.initialize_solver(
        initial_guess=state,
        solver_options={'relaxation': 0.25, 'max_iterations': 2,
                        'require_chemistry_to_converge': True}
    )
    solver = model.solver
    assert isinstance(solver, coupled.FixedPointSolver)
    assert restored_states == [state]
    assert model.iterate() == (False, False)
    np.testing.assert_allclose(solver.fixed_point_x, initial)
    np.testing.assert_allclose(solver.fixed_point_gx, initial + [0.2, 20.0, 10.0])
    np.testing.assert_allclose(solver.temperature_residual, 20.0/300.0)
    np.testing.assert_allclose(solver.log10pressure_residual, 0.2)
    assert solver.chemistry_converged
    assert solver.chemistry_diagnostics == {'converged': True, 'steps': 1}
    assert len(solver.solver_history) == 1
    assert model.iterate() == (True, False)
    np.testing.assert_allclose(trials[1], initial + [0.05, 5.0, 2.5])
    assert model.iterate() == (True, False)
    assert len(trials) == solver.niterations == 2

    def identity_map(x, require_chemistry_to_converge, *, chemistry_diagnostics):
        chemistry_diagnostics.update(converged=False, steps=7000)
        return x.copy()

    model._fixed_point_map = identity_map
    assert model.find_steady_state(initial_guess=state)
    assert model.solver is not solver
    assert model.solver.niterations == 1
    assert not model.solver.chemistry_converged
    assert model.iterate() == (False, True)
    assert model.solver.niterations == 1

    def failing_map(x, require_chemistry_to_converge, *, chemistry_diagnostics):
        chemistry_diagnostics.update(converged=True, steps=3)
        raise coupled.ClimaException('Local test climate failure.')

    model._fixed_point_map = failing_map
    model.initialize_solver(initial_guess=state)
    try:
        model.iterate()
    except coupled.ClimaException:
        pass
    else:
        raise AssertionError('Map failure did not propagate.')
    assert model.solver.chemistry_converged
    assert model.solver.chemistry_diagnostics['steps'] == 3
    assert model.solver.niterations == 0
    try:
        model.iterate()
    except RuntimeError:
        pass
    else:
        raise AssertionError('A failed solver session was reused.')
    print('Local fixed-point solver checks passed.', flush=True)


def assert_state_equal(actual, expected):
    """Compare nested state dictionaries without coercing metadata arrays."""
    assert actual.keys() == expected.keys()
    for key in actual:
        if isinstance(actual[key], dict):
            assert_state_equal(actual[key], expected[key])
        else:
            np.testing.assert_equal(actual[key], expected[key])


def run_restart_checks():
    """Local state checks; requires one initial climate solve, no chemistry."""
    model = initialize()
    model.c.verbose = False
    guess = model.initial_guess(edd=1e5)
    assert set(guess['climate']) == {'pressure', 'temperature', 'convecting_with_below'}
    assert model.solver is None
    assert guess['edd'].shape == guess['edd_pressure'].shape == (2,)
    assert_state_equal(model.model_state(), guess)

    # Reinitialization preserves atmospheric fields and resets bookkeeping.
    model.initialize_solver(initial_guess=guess)
    assert model.solver.niterations == 0
    assert not model.solver.converged
    assert model.solver.chemistry_converged is None
    restored = model.model_state()
    assert_state_equal(restored, guess)
    restored['climate']['temperature'][0] += 10.0
    restored['edd'][0] *= 2.0
    assert_state_equal(model.model_state(), guess)

    fresh = initialize()
    fresh.initialize_solver(initial_guess=guess)
    assert_state_equal(fresh.model_state(), guess)
    np.testing.assert_equal(fresh.pc.wrk.usol, model.pc.wrk.usol)
    # A second capture before any climate rebuild must use restored guesses.
    round_trip = fresh.model_state()
    fresh.initialize_solver(initial_guess=round_trip)
    assert_state_equal(fresh.model_state(), guess)

    # Validation failures must leave the existing atmosphere and solver alone.
    saved_solver = fresh.solver
    invalid_states = []
    for key, value in [('schema_version', 999), ('edd', np.array([-1., 1.]))]:
        invalid = deepcopy(guess)
        invalid[key] = value
        invalid_states.append(invalid)
    for key, value in [('pressure', np.array([1e6, 1.])),
                       ('temperature', np.array([np.nan])),
                       ('convecting_with_below', np.zeros(50, dtype=int))]:
        invalid = deepcopy(guess)
        invalid['climate'][key] = value
        invalid_states.append(invalid)
    for factor in (1.0, 2.0):
        invalid = deepcopy(guess)
        invalid['photochem']['toa_target_pressure'] = np.array(
            factor * guess['climate']['pressure'][0])
        invalid_states.append(invalid)
    invalid = deepcopy(guess)
    del invalid['photochem']['temperature']
    invalid_states.append(invalid)
    for invalid in invalid_states:
        try:
            fresh.initialize_solver(initial_guess=invalid)
        except (ValueError, coupled.PhotoException):
            pass
        else:
            raise AssertionError('Invalid coupled state was accepted.')
        assert fresh.solver is saved_solver
        assert_state_equal(fresh.model_state(), guess)
    print('Local coupled restart checks passed.', flush=True)


def run_grid_checks():
    """Check shared-grid coordinates and climate guesses without model solves."""
    model = coupled.CoupledPhotochemistryClimate.__new__(coupled.CoupledPhotochemistryClimate)
    model.c = SimpleNamespace(P_top=1.0, T=np.zeros(4))
    maintenance = SimpleNamespace(target_pressure=0.1)
    model.pc = SimpleNamespace(var=SimpleNamespace(toa_pressure_maintenance=maintenance))

    grid = model._make_pressure_grid(6.0)
    assert grid.size == model.c.T.size + 1
    np.testing.assert_allclose(grid, [6., 4.25, 2.5, 0.75, -1.])
    # For a linear T(log10 P), interpolation has an analytical answer.
    x = np.concatenate(([6.0], 200.0 + 10.0*grid))
    T_init = model._climate_temperature_guess(x)
    np.testing.assert_allclose(T_init, [260., 252.5, 237.5, 222.5, 207.5])
    assert T_init.size == model.c.T.size + 1
    assert not np.allclose(T_init[1:], 0.5*(x[1:-1] + x[2:]))

    # Equal domain tops recover the original midpoint averaging exactly.
    maintenance.target_pressure = model.c.P_top
    grid = model._make_pressure_grid(6.0)
    x = np.concatenate(([6.0], 200.0 + 10.0*grid))
    np.testing.assert_allclose(model._climate_temperature_guess(x),
                               np.concatenate(([x[1]], 0.5*(x[1:-1] + x[2:]))))

    # A shallower chemical domain uses constant endpoint extrapolation for
    # the remaining climate layers.
    maintenance.target_pressure = 10.0
    grid = model._make_pressure_grid(6.0)
    x = np.concatenate(([6.0], 200.0 + 10.0*grid))
    np.testing.assert_allclose(model._climate_temperature_guess(x),
                               [260., 252.5, 237.5, 222.5, 210.])

    # Derive the vector after the saved chemical top has been restored.
    maintenance.target_pressure = 0.01
    climate = {'pressure': np.array([1e6, 1e5, 1e3, 1e1, 2.]),
               'temperature': np.array([300., 280., 240., 200., 180.])}
    vector = model._fixed_point_state_from_climate(climate)
    assert vector.size == model.c.T.size + 2
    np.testing.assert_allclose(vector[[0, 1, -1]], [6., 300., 180.])
    np.testing.assert_allclose(model._make_pressure_grid(vector[0])[-1], -2.0)
    assert maintenance.target_pressure == 0.01
    print('Local coupling-grid checks passed.', flush=True)


def run_acceleration_checks():
    """Local inexpensive checks of Anderson convergence and fallbacks."""
    def make_solver(mapping, memory, **options):
        calls = []

        def evaluate(x, require_chemistry_to_converge, *, chemistry_diagnostics):
            calls.append(x.copy())
            chemistry_diagnostics['converged'] = True
            return mapping(x)

        physical_model = SimpleNamespace(
            _fixed_point_map=evaluate, c=SimpleNamespace(P_top=1.0),
            pc=SimpleNamespace(var=SimpleNamespace(
                toa_pressure_maintenance=SimpleNamespace(target_pressure=0.1))))
        solver = coupled.FixedPointSolver(
            physical_model, anderson_memory=memory, max_iterations=500, **options)
        return solver, calls

    root = np.array([6.0, 300.0, 200.0])
    mapping = lambda x: root + np.array([0.3, 0.9, 0.7]) * (x - root)
    counts = []
    for memory in (0, 3):
        solver, calls = make_solver(mapping, memory, temperature_tol=1e-8,
                                    log10pressure_tol=1e-8)
        initial = np.array([6.1, 320., 210.])
        solver._initialize(initial)
        assert solver.iterate() == (False, False)
        np.testing.assert_allclose(solver._next_x, initial + 0.5*(mapping(initial)-initial))
        while True:
            give_up, converged = solver.iterate()
            assert len(solver._anderson_history) <= memory + 1
            if give_up or converged:
                break
        assert converged and not give_up
        np.testing.assert_allclose(solver.fixed_point_x, root, rtol=2e-7, atol=2e-7)
        assert len(calls) == solver.niterations
        solver.iterate()
        assert len(calls) == solver.niterations
        counts.append(solver.niterations)
    assert counts[1] < counts[0]

    # Singular residual differences must fall back to damping.
    solver, _ = make_solver(lambda x: x + [0., 10., 10.], 3)
    solver._initialize(root)
    solver.iterate()
    solver.iterate()
    assert solver.solver_history[-1]['anderson_reset'] == 'dependent_history'
    assert solver.solver_history[-1]['step_method'] == 'relaxation'

    # The secant fit predicts a negative temperature; never evaluate that trial.
    solver, calls = make_solver(lambda x: np.array([6., 2*x[1]+10., 200.]), 3)
    solver._initialize(np.array([6., 100., 200.]))
    solver.iterate()
    solver.iterate()
    assert solver.solver_history[-1]['anderson_reset'] == 'invalid_proposal'
    assert np.all(solver._next_x[1:] > 0.)
    assert len(calls) == 2

    solver, _ = make_solver(lambda x: np.array([6., 5*x[1], 200.]), 3)
    solver._initialize(root)
    solver.iterate()
    solver.iterate()
    assert solver.solver_history[-1]['anderson_reset'] == 'residual_growth'
    print(f'Local Anderson checks passed (affine map: {counts[0]} damped vs {counts[1]} accelerated evaluations).',
          flush=True)


def run_benchmark(output=None):
    """Compare damping and Anderson from one saved atmosphere, outside CI.

    Both runs use fresh compatible models and identical solver/chemistry
    settings. Shared initial-guess construction and model loading are excluded
    from timings; state restoration and all coupling evaluations are included.
    """
    source = initialize()
    source.pc.var.verbose = 0
    source.pc.var.atol = 1e-18
    source.c.verbose = False
    guess = source.initial_guess(edd=1e5)
    results = {}
    solutions = {}
    for label, memory in [('damped', 0), ('anderson', 3)]:
        model = initialize()
        model.pc.var.verbose = 0
        model.pc.var.atol = 1e-18
        model.c.verbose = False
        print(f'Starting {label} benchmark from the shared initial state.', flush=True)
        started = time.monotonic()
        model.initialize_solver(
            initial_guess=guess,
            solver_options={'relaxation': 0.5, 'anderson_memory': memory,
                            'require_chemistry_to_converge': True, 'verbose': True})
        total_steps = total_errors = total_resyncs = 0
        while True:
            give_up, converged = model.iterate()
            chemistry = model.solver.chemistry_diagnostics
            total_steps += chemistry['steps']
            total_errors += chemistry['integration_errors']
            total_resyncs += chemistry['resync_failures']
            if give_up or converged:
                break
        seconds = time.monotonic() - started
        solver = model.solver
        assert converged and not give_up, f'{label} did not converge.'
        assert solver.chemistry_converged
        solution = model.solution()
        assert all(np.all(np.isfinite(values)) for values in solution.values())
        assert np.all(solution['temp'] > 0.)
        assert np.all(solution['density'] > 0.)
        solutions[label] = solution
        results[label] = {
            'anderson_memory': memory,
            'iterations': solver.niterations,
            'seconds': seconds,
            'temperature_residual': solver.temperature_residual,
            'log10pressure_residual': solver.log10pressure_residual,
            'chemistry_converged': solver.chemistry_converged,
            'chemical_steps': total_steps,
            'integration_errors': total_errors,
            'resync_failures': total_resyncs,
            'climate_surface_temperature': float(model.c.T_surf),
            'photochemical_surface_pressure_bar': float(model.pc.wrk.surface_pressure),
            'history': solver.solver_history,
        }
        print(f'{label}: converged in {solver.niterations} evaluations, {seconds:.1f} seconds.', flush=True)

    # Compare photochemical temperatures at common physical pressures.
    damped, accelerated = solutions['damped'], solutions['anderson']
    low = max(damped['pressure'][-1], accelerated['pressure'][-1])
    high = min(damped['pressure'][0], accelerated['pressure'][0])
    logP = np.linspace(np.log10(low), np.log10(high), 100)
    temperatures = [np.interp(logP, np.log10(sol['pressure'][::-1]), sol['temp'][::-1])
                    for sol in (damped, accelerated)]
    results['comparison'] = {
        'runtime_speedup': results['damped']['seconds']/results['anderson']['seconds'],
        'map_evaluations_saved': results['damped']['iterations']-results['anderson']['iterations'],
        'max_relative_photochemical_temperature_difference': float(np.max(
            np.abs(temperatures[1]-temperatures[0])/temperatures[0])),
        'log10_photochemical_surface_pressure_difference': float(abs(np.log10(
            results['anderson']['photochemical_surface_pressure_bar']/
            results['damped']['photochemical_surface_pressure_bar']))),
    }
    summary = json.dumps(results, indent=2)
    if output is not None:
        Path(output).write_text(summary + '\n')
    print(summary, flush=True)
    return results


def run_full():
    """Run the expensive coupled convergence check locally, outside CI."""
    model = initialize()
    model.pc.var.verbose = 0
    model.pc.var.atol = 1.0e-18
    model.c.verbose = False
    started = time.monotonic()
    guess = model.initial_guess(edd=1e5)
    model.initialize_solver(initial_guess=guess, method="fixed_point",
                            solver_options={"relaxation": 0.5})
    solver = model.solver
    print('Initialized coupled solver with default controls.', flush=True)
    while True:
        print(f'Starting coupled iteration {solver.niterations + 1}.', flush=True)
        give_up, converged = model.iterate()
        chemistry = solver.chemistry_diagnostics
        print(
            f'Iteration {solver.niterations}: '
            f'T residual={solver.temperature_residual:.6g}, '
            f'log10P residual={solver.log10pressure_residual:.6g} dex, '
            f'chemistry converged={solver.chemistry_converged}, '
            f'chemical steps={chemistry["steps"]}, '
            f'integration errors={chemistry["integration_errors"]}, '
            f'resync failures={chemistry["resync_failures"]}, '
            f'elapsed={time.monotonic() - started:.1f} s',
            flush=True,
        )
        if give_up or converged:
            break

    print(
        f'Final: coupled converged={converged}, '
        f'chemistry converged={solver.chemistry_converged}, '
        f'climate T_surf={model.c.T_surf:.3f} K, '
        f'photochemical surface pressure={model.pc.wrk.surface_pressure:.6g} bar',
        flush=True,
    )
    assert converged and not give_up, 'Coupled solver exhausted its iteration limit.'
    assert solver.temperature_residual <= solver.temperature_tol
    assert solver.log10pressure_residual <= solver.log10pressure_tol
    solution = model.solution()
    assert all(np.all(np.isfinite(values)) for values in solution.values())
    assert np.all(solution['temp'] > 0.0)
    assert np.all(solution['pressure'] > 0.0)
    assert np.all(solution['density'] > 0.0)
    # Restart a fresh compatible model from the converged atmosphere.
    state = model.model_state()
    fresh = initialize()
    fresh.pc.var.verbose = 0
    fresh.pc.var.atol = 1.0e-18
    fresh.c.verbose = False
    fresh.initialize_solver(initial_guess=state)
    assert_state_equal(fresh.model_state(), state)
    assert fresh.solver.niterations == 0
    while True:
        give_up, restarted_converged = fresh.iterate()
        if give_up or restarted_converged:
            break
    assert restarted_converged and not give_up, 'Converged warm start failed to reconverge.'
    assert fresh.solver.chemistry_converged
    print('Full coupled-model and converged restart tests passed.', flush=True)


def main():
    parser = argparse.ArgumentParser(description='Run the coupled-model smoke test by default.')
    parser.add_argument('--full', action='store_true',
                        help='Run the expensive full coupled convergence test locally.')
    parser.add_argument('--restart', action='store_true',
                        help='Run local restart checks with one climate solve and no chemistry integration.')
    parser.add_argument('--benchmark', action='store_true',
                        help='Compare plain damping and Anderson on the coupled test, outside CI.')
    parser.add_argument('--benchmark-output', help='Optional JSON file for benchmark results.')
    parser.add_argument('--grid', action='store_true',
                        help='Run inexpensive coupling-grid interpolation checks without model solves.')
    args = parser.parse_args()
    if args.grid:
        run_grid_checks()
    elif args.benchmark:
        run_acceleration_checks()
        run_benchmark(args.benchmark_output)
    elif args.full:
        run_grid_checks()
        run_acceleration_checks()
        run_solver_checks()
        run_restart_checks()
        run_full()
    elif args.restart:
        run_grid_checks()
        run_solver_checks()
        run_restart_checks()
    else:
        test_smoke()
        print('Coupled-model smoke test passed.')


if __name__ == '__main__':
    main()
