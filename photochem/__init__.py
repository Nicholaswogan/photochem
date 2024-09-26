import os
os.environ['OMP_NUM_THREADS'] = '1'

from ._photochem import EvoAtmosphere, Atmosphere, PhotoException, __version__
from photochem_clima_data import DATA_DIR

zahnle_earth = DATA_DIR+'/reaction_mechanisms/zahnle_earth.yaml'