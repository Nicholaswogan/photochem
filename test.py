import Photochem
Photochem.photochem_setup("data/reaction_mechanisms/zahnle_earth.yaml",\
                          "templates/ModernEarth/settings_ModernEarth.yaml",\
                          "templates/ModernEarth/Sun_now.txt",\
                          "templates/ModernEarth/atmosphere_ModernEarth.txt")          
rtol = 1.0e-3
atol = 1.0e-25
maxsteps = 10000
success, err = Photochem.photochem.photo_equilibrium(maxsteps, rtol, atol)

usol_out = Photochem.photochem_vars.usol_out
surface_fluxs, err = Photochem.photochem.compute_surface_fluxes(usol_out)
print(surface_fluxs)