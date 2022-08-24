import os
os.environ['OMP_NUM_THREADS'] = '1'

from ._photochem import EvoAtmosphere, Atmosphere, __version__

zahnle_earth = os.path.dirname(os.path.realpath(__file__))+ \
                '/data/reaction_mechanisms/zahnle_earth.yaml'