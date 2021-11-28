from Photochem import Atmosphere, zahnle_earth

pc = Atmosphere(zahnle_earth,\
                "../templates/ModernEarth/settings_ModernEarth.yaml",\
                "../templates/ModernEarth/Sun_now.txt",\
                "../templates/ModernEarth/atmosphere_ModernEarth.txt")
success = pc.photochemical_equilibrium()
