from photochem import Atmosphere, zahnle_earth

def main():

    pc = Atmosphere(zahnle_earth,\
                    "../templates/ModernEarth/settings_ModernEarth.yaml",\
                    "../templates/ModernEarth/Sun_now.txt",\
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt")

    pc.initialize_stepper(pc.wrk.usol)
    tn = pc.step()
    pc.destroy_stepper()

if __name__ == "__main__":
    main()