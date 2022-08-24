from photochem import EvoAtmosphere, zahnle_earth
import numpy as np

from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)

pc = EvoAtmosphere(zahnle_earth,\
                    "testevo_settings.yaml",\
                    "../templates/ModernEarth/Sun_now.txt",\
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt")

t_eval = np.logspace(0,15,10)
success = pc.evolve('test.dat',0.0, pc.var.usol_init, t_eval, True)