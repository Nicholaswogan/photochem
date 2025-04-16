from ._photochem import EvoAtmosphere, PhotoException, __version__
from photochem_clima_data import DATA_DIR

# Limits OpenMP threads to 1
from threadpoolctl import threadpool_limits
_ = threadpool_limits(limits=1, user_api='openmp')

zahnle_earth = DATA_DIR+'/reaction_mechanisms/zahnle_earth.yaml'