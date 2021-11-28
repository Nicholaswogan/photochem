from ._Photochem import Atmosphere
from .FormatReactions import FormatReactions, FormatSettings

import os
zahnle_earth = os.path.dirname(os.path.realpath(__file__))+ \
                '/data/reaction_mechanisms/zahnle_earth.yaml'