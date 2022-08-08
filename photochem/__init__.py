from ._photochem import EvoAtmosphere, Atmosphere, sat_pressure_H2O, __version__

import os
zahnle_earth = os.path.dirname(os.path.realpath(__file__))+ \
                '/data/reaction_mechanisms/zahnle_earth.yaml'