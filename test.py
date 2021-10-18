from Photochem import Atmosphere
pc = Atmosphere("data/reaction_mechanisms/zahnle_earth.yaml",\
                "templates/ModernEarth/settings_ModernEarth.yaml",\
                "templates/ModernEarth/Sun_now.txt",\
                "templates/ModernEarth/atmosphere_ModernEarth.txt")
success = pc.photochemical_equilibrium()