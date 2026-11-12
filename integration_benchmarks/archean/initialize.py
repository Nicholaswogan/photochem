"""Fresh and saved-state initialization for the original Archean benchmark."""
import importlib.util
from pathlib import Path
import pickle
from photochem import EvoAtmosphere, EvoAtmosphereState
from photochem.utils import stars

HERE = Path(__file__).resolve().parent
MECHANISM = 'mechanism1.yaml'
SETTINGS = 'settings1.yaml'
CONTROLS = dict(atol=1e-20, rtol=1e-3, equilibrium_time=1e17,
                smooth_factor=1., k_cond=100., k_evap=0.,
                upwind_molec_diff=False, jacobian_method=3)


def initialize(destination, fresh=False):
    stars.solar_spectrum(outputfile=str(destination / 'flux.txt'), stellar_flux=1361., age=3.)
    if fresh:
        spec = importlib.util.spec_from_file_location('archean_fresh', HERE / 'benchmark1.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.initialize()
    pc = EvoAtmosphere(str(HERE / 'mechanism1.yaml'), str(HERE / 'settings1.yaml'),
                       str(destination / 'flux.txt'))
    with (HERE / 'state1.pkl').open('rb') as stream:
        pc.initialize_from_state(EvoAtmosphereState.from_dict(pickle.load(stream)))
    return pc
