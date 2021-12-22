from Photochem import Atmosphere, zahnle_earth

pc = Atmosphere(zahnle_earth,\
                "../templates/ModernEarth/settings_ModernEarth.yaml",\
                "../templates/ModernEarth/Sun_now.txt",\
                "../templates/ModernEarth/atmosphere_ModernEarth.txt")

pc.initialize_stepper(pc.wrk.usol)
tn = 0
while tn<1e17:
    tn = pc.step()
pc.destroy_stepper()

con = pc.redox_conservation()
print('redox conservation = '+'%.2e'%con)
assert con < 1e-4

con = pc.atom_conservation('C').factor
print('C conservation = '+'%.2e'%con)
assert con < 1e-6

con = pc.atom_conservation('S').factor
print('S conservation = '+'%.2e'%con)
assert con < 1e-6

con = pc.atom_conservation('H').factor
print('H conservation = '+'%.2e'%con)
assert con < 1e-6

con = pc.atom_conservation('O').factor
print('O conservation = '+'%.2e'%con)
assert con < 1e-6

con = pc.atom_conservation('Cl').factor
print('Cl conservation = '+'%.2e'%con)
assert con < 1e-6