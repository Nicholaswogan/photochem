import argparse
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
                    {'solver_options': {'relaxation': 0.0}}):
        try:
            model.initialize_solver(edd=1e5, **options)
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
    initialization_inputs = []
    trials = []

    def initial_guess(edd, **kwargs):
        initialization_inputs.append((edd, kwargs))
        return initial.copy()

    def fixed_point_map(x, require_chemistry_to_converge, *, chemistry_diagnostics):
        assert require_chemistry_to_converge
        trials.append(x.copy())
        chemistry_diagnostics.update(converged=True, steps=1)
        return x + np.array([0.2, 20.0, 10.0])

    model._initial_guess = initial_guess
    model._fixed_point_map = fixed_point_map
    model.initialize_solver(
        edd=1e5, mode='periodic',
        solver_options={'relaxation': 0.25, 'max_iterations': 2,
                        'require_chemistry_to_converge': True}
    )
    solver = model.solver
    assert isinstance(solver, coupled.FixedPointSolver)
    assert initialization_inputs[0][0] == 1e5
    assert initialization_inputs[0][1]['mode'] == 'periodic'
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
    assert model.find_steady_state(edd=1e5)
    assert model.solver is not solver
    assert model.solver.niterations == 1
    assert not model.solver.chemistry_converged
    assert model.iterate() == (False, True)
    assert model.solver.niterations == 1

    def failing_map(x, require_chemistry_to_converge, *, chemistry_diagnostics):
        chemistry_diagnostics.update(converged=True, steps=3)
        raise coupled.ClimaException('Local test climate failure.')

    model._fixed_point_map = failing_map
    model.initialize_solver(edd=1e5)
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


def run_full():
    """Run the expensive coupled convergence check locally, outside CI."""
    model = initialize()
    model.pc.var.verbose = 0
    model.pc.var.atol = 1.0e-18
    model.c.verbose = False
    started = time.monotonic()
    model.initialize_solver(edd=1e5, method="fixed_point",
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
    print('Full coupled-model test passed.', flush=True)


def main():
    parser = argparse.ArgumentParser(description='Run the coupled-model smoke test by default.')
    parser.add_argument('--full', action='store_true',
                        help='Run the expensive full coupled convergence test locally.')
    args = parser.parse_args()
    if args.full:
        run_solver_checks()
        run_full()
    else:
        test_smoke()
        print('Coupled-model smoke test passed.')


if __name__ == '__main__':
    main()
