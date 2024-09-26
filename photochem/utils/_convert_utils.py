import os
import h5py
from ._format import Loader, MyDumper, yaml, FormatReactions_main
from photochem_clima_data import DATA_DIR
    
def compare2reactions(rx1, rx2):
    rx1 = rx1.replace('<=>','=>').replace('(','').replace(')','')
    rx2 = rx2.replace('<=>','=>').replace('(','').replace(')','')
    react1, prod1 = [sorted(a.replace(' ','').split('+')) for a in rx1.split('=>')]
    react2, prod2 = [sorted(a.replace(' ','').split('+')) for a in rx2.split('=>')]
    return all([react1 == react2, prod1 == prod2])

def generate_photo_yaml_entries(species_list):
    species_set = set(species_list)
    
    filenames = [a for a in os.listdir(DATA_DIR+'/xsections') if '.h5' in a and a != 'bins.h5']
    all_photoreactions = []
    for filename in filenames:
        with h5py.File(DATA_DIR+'/xsections/'+filename,'r') as f:
            if 'photodissociation-qy' in f.keys():
                for key in f['photodissociation-qy'].keys():
                    if key == 'wavelengths':
                        continue
                    all_photoreactions.append(key)
    rx_list = []
    for rx in all_photoreactions:
        tmp1, tmp2 = rx.split('=>')
        tmp = set([tmp1.split('+')[0].strip()] + [a.strip() for a in tmp2.split('+')])
        if tmp.issubset(species_set):
            entry = {
                'equation': rx,
                'type': 'photolysis'
            }
            rx_list.append(entry)
    return rx_list

def sort_photos(data_photo, possible_photo):
    
    missing = []
    not_missing = []

    for rx in data_photo:
        found = False
        for ph in possible_photo:
            if compare2reactions(ph['equation'],rx['equation']):
                found = True
                not_missing.append(ph)
                break
        if not found:
            missing.append(rx)
            
    return missing, not_missing

def resave_mechanism_with_atoms(infile, outfile, atoms_names):
    
    with open(infile,'r') as f:
        dat = yaml.load(f,Loader=Loader)

    atoms = []
    exclude_atoms = []
    for i,at in enumerate(dat['atoms']):
        if at['name'] in atoms_names:
            atoms.append(at)
        else:
            exclude_atoms.append(at['name'])

    species = []
    comp = {}
    comp['hv'] = []
    comp['M'] = []
    for i,sp in enumerate(dat['species']):
        comp[sp['name']] = [key for key in sp['composition'] if sp['composition'][key] > 0]

        exclude = False
        for tmp in comp[sp['name']]:
            if tmp in exclude_atoms:
                exclude = True
                break
        if not exclude:
            species.append(sp)

    if "particles" in dat:
        particles = []
        for i,sp in enumerate(dat['particles']):
            comp_tmp = [key for key in sp['composition'] if sp['composition'][key] > 0]

            exclude = False
            for tmp in comp_tmp:
                if tmp in exclude_atoms:
                    exclude = True
                    break
            if not exclude:
                particles.append(sp)

    reactions = []
    for i,rx in enumerate(dat['reactions']):
        eq = rx['equation']
        eq = eq.replace('(','').replace(')','')
        if '<=>' in eq:
            split_str = '<=>'
        else:
            split_str = '=>'

        a,b = eq.split(split_str)
        a = a.split('+')
        b = b.split('+')
        a = [a1.strip() for a1 in a]
        b = [b1.strip() for b1 in b]
        sp = a + b

        exclude = False
        for s in sp:
            for tmp in comp[s]:
                if tmp in exclude_atoms:
                    exclude = True
                    break
            if exclude:
                break
        if not exclude:
            reactions.append(rx)
                
    out = dat
    out['atoms'] = atoms
    out['species'] = species
    if 'particles' in dat:
        out['particles'] = particles
    out['reactions'] = reactions

    out = FormatReactions_main(out)

    with open(outfile,'w') as f:
        yaml.dump(out,f,Dumper=MyDumper,sort_keys=False,width=70)

class AtomData():
    atoms = {}
    atoms['H'] = {"name": "H", "mass": 1.00797, "redox": -0.5}
    atoms['N'] = {"name": "N", "mass": 14.0067, "redox": 0}
    atoms['O'] = {"name": "O", "mass": 15.9994, "redox": 1}
    atoms['C'] = {"name": "C", "mass": 12.011, "redox": -2}
    atoms['S'] = {"name": "S", "mass": 32.06, "redox": -2}
    atoms['Cl'] = {"name": "Cl", "mass": 35.453, "redox": -1}
    atoms['He'] = {"name": "He", "mass": 4.0026, "redox": 0}
    atoms['P'] = {"name": "P", "mass": 30.9738, "redox": 0}
    atoms['Na'] = {"name": "Na", "mass": 22.9897, "redox": 0}
    atoms['K'] = {"name": "K", "mass": 39.0983, "redox": 0}
    atoms['Si'] = {"name": "Si", "mass": 28.0855, "redox": 0}
    atoms['Fe'] = {"name": "Fe", "mass": 55.845, "redox": 0}
    atoms['Ar'] = {"name": "Ar", "mass": 39.948, "redox": 0}
    atoms['Ti'] = {"name": "Ti", "mass": 47.867, "redox": 0}
    atoms['V'] = {"name": "V", "mass": 50.9415, "redox": 0}
    atoms['Mg'] = {"name": "Mg", "mass": 54.938, "redox": 0}
    atoms['Ca'] = {"name": "Ca", "mass": 40.078, "redox": 0}

