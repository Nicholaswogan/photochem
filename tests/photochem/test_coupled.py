import argparse
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


def main():
    parser = argparse.ArgumentParser(description='Run the coupled-model smoke test by default.')
    parser.add_argument('--full', action='store_true',
                        help='Run the full coupled test locally (not implemented yet).')
    args = parser.parse_args()
    if args.full:
        parser.error('The full coupled test has not been implemented yet.')
    test_smoke()
    print('Coupled-model smoke test passed.')


if __name__ == '__main__':
    main()
