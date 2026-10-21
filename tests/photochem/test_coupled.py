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
    assert not model.pc.atmosphere_initialized
    assert model.c.P_top == 1.0


def run_full():
    """Run the expensive coupled convergence check locally, outside CI."""
    model = initialize()
    model.pc.var.verbose = 0
    model.pc.var.atol = 1.0e-18
    model.c.verbose = False
    started = time.monotonic()
    model.initialize_solver(edd=1e5)
    print('Initialized coupled solver with default controls.', flush=True)
    while True:
        print(f'Starting coupled iteration {model.niterations + 1}.', flush=True)
        give_up, converged = model.iterate()
        chemistry = model.chemistry_diagnostics
        print(
            f'Iteration {model.niterations}: '
            f'T residual={model.temperature_residual:.6g}, '
            f'log10P residual={model.log10pressure_residual:.6g} dex, '
            f'chemistry converged={model.chemistry_converged}, '
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
        f'chemistry converged={model.chemistry_converged}, '
        f'climate T_surf={model.c.T_surf:.3f} K, '
        f'photochemical surface pressure={model.pc.wrk.surface_pressure:.6g} bar',
        flush=True,
    )
    assert converged and not give_up, 'Coupled solver exhausted its iteration limit.'
    assert model.temperature_residual <= model.temperature_tol
    assert model.log10pressure_residual <= model.log10pressure_tol
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
        run_full()
    else:
        test_smoke()
        print('Coupled-model smoke test passed.')


if __name__ == '__main__':
    main()
