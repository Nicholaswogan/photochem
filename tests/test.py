from photochem import Atmosphere, EvoAtmosphere, zahnle_earth

def main():

    pc = Atmosphere(zahnle_earth,\
                    "../templates/ModernEarth/settings_ModernEarth.yaml",\
                    "../templates/ModernEarth/Sun_now.txt",\
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt")

    pc.initialize_stepper(pc.wrk.usol)
    tn = pc.step()
    pc.destroy_stepper()

    pc1 = EvoAtmosphere(zahnle_earth,\
                    "testevo_settings2.yaml",\
                    "../templates/ModernEarth/Sun_now.txt",\
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt")
    
    pc1.initialize_stepper(pc1.wrk.usol)
    tn = pc1.step()
    pc1.destroy_stepper()

if __name__ == "__main__":
    main()