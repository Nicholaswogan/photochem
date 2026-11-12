import pickle

import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from photochem.extensions import coupled

from photochem.utils import stars
from photochem.utils import zahnle_rx_and_thermo_files

TEST_DIR = Path(__file__).resolve().parent

def input_files(output_dir, test_dir):
    output_dir = Path(output_dir)
    # zahnle_rx_and_thermo_files(
    #     atoms_names=['H','N','O','C'],
    #     rxns_filename=str(test_dir / 'mechanism1.yaml'),
    #     thermo_filename=None,
    # )
    stars.solar_spectrum(
        outputfile=str(output_dir / 'flux.txt'),
        stellar_flux=1361.0,
        age=3.0
    )

def initialize():
    with TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        input_files(output_dir, TEST_DIR)
        model = coupled.CoupledPhotochemistryClimate(
            mechanism_file=str(TEST_DIR / 'mechanism1.yaml'),
            settings_file=str(TEST_DIR / 'settings1.yaml'),
            flux_file=str(output_dir / 'flux.txt'),
        )

    model.pc.var.atol = 1.0e-20
    model.pc.var.conv_longdy = 0.0
    model.pc.var.equilibrium_time = 1.0e17
    for i in range(len(model.pc.var.cond_params)):
        model.pc.var.cond_params[i].smooth_factor = 1.0
        model.pc.var.cond_params[i].k_cond = 100.0
        model.pc.var.cond_params[i].k_evap = 0.0
    
    guess = model.initial_guess(1.0e5)
    radii = model.pc.var.particle_radius
    for species in ('CO2aer', 'H2Oaer'):
        i = model.pc.dat.species_names.index(species)
        radii[i, :] = 1.0e-2
    model.pc.var.particle_radius = radii
    model.pc.var.verbose = 1
    model.pc.var.max_order = 1

    return model.pc

def run():
    pc = initialize()
    converged = pc.find_steady_state()
    print(converged)
    # with open(str(TEST_DIR / "state1.pkl"), 'wb') as f:
    #     pickle.dump(pc.model_state().to_dict(), f)
    
if __name__ == "__main__":
    run()
        
