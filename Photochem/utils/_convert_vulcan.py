import numpy as np

from ._format import FormatReactions_main, MyDumper, Loader, yaml
from ._convert_utils import AtomData, compare2reactions, \
                            generate_photo_yaml_entries, sort_photos
                            
def vulcan2yaml(vulcan_rx_filename, all_compose_filename, outfile = None, \
                rename_species = True, photo_database = "Photochem", vulcan_nasa9_data_folder = None):
    
    # get composition of each species
    fil = open(all_compose_filename,'r')
    lines = fil.readlines()
    fil.close()
    
    if rename_species:
        to_replace = [['O_1','O1D'],
                      ['N_2D','N2D'],
                      ['CH2_1','1CH2']]
    else:
        to_replace = []
        
    if vulcan_nasa9_data_folder == None:
        rx_arrow = "=>"
    else:
        rx_arrow = "<=>"
        
    atoms = lines[0].split()[1:17]
    composition = {}

    for line in lines[1:]:
        sp = line.split()[0]
        for rep in to_replace:
            if sp == rep[0]:
                sp = rep[1]
                break
        composition[sp] = {}
        na = line.split()[1:17]
        for i,a in enumerate(na):
            if int(a) > 0:
                composition[sp][atoms[i]] = int(a)
    
    # parse reactions
    fil = open(vulcan_rx_filename,'r')
    lines = fil.readlines()
    fil.close()

    type_ind = {}
    for i, line in enumerate(lines):
        if line.startswith("# Two-body"):
            type_ind['elementary'] = i + 3
        elif line.startswith("# 3-body and Disscoiation Reactions"):
            type_ind['falloff'] = i + 3
        elif line.startswith("# 3-body reactions without high-pressure rates"):
            type_ind['three-body'] = i + 3
        elif line.startswith("# photo disscoiation "):
            type_ind['photolysis'] = i + 2

    reactions = []

    # elementary
    for i, line in enumerate(lines[type_ind['elementary']:]):
        if len(line.strip()) == 0:
            break
        elif line[0] == '#':
            break
        rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
        for rep in to_replace:
            if rep[0] in rx_str:
                rx_str = rx_str.replace(rep[0],rep[1])
        A, b, Ea = [float(a) for a in line.split(']')[1].split()[:3]]
        rx = {}
        rx['equation'] = rx_str
        rx['rate-constant'] = {"A": A, "b": b, "Ea": Ea}
        reactions.append(rx)
        
    # falloff
    for i, line in enumerate(lines[type_ind['falloff']:]):
        if len(line.strip()) == 0:
            break
        elif line[0] == '#':
            break
        rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
        for rep in to_replace:
            if rep[0] in rx_str:
                rx_str = rx_str.replace(rep[0],rep[1])
        if "M" not in rx_str.split(rx_arrow)[0]:
            rx_str = rx_str.split(rx_arrow)[0] + "+ M "+rx_arrow+rx_str.split(rx_arrow)[1]
        if "M" not in rx_str.split(rx_arrow)[1]:
            rx_str = rx_str.split(rx_arrow)[0]+rx_arrow+rx_str.split(rx_arrow)[1]+ "+ M "
            
        A0, b0, Ea0, Ainf, binf, Eainf = [float(a) for a in line.split(']')[1].split()[:6]]
        rx = {}
        rx['equation'] = rx_str
        rx['type'] = "falloff"
        rx['low-P-rate-constant'] = {"A": A0, "b": b0, "Ea": Ea0}
        rx['high-P-rate-constant'] = {"A": Ainf, "b": binf, "Ea": Eainf}
        reactions.append(rx)
        
    # three-body
    for i, line in enumerate(lines[type_ind['three-body']:]):
        if len(line.strip()) == 0:
            break
        elif line[0] == '#':
            break
        rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
        for rep in to_replace:
            if rep[0] in rx_str:
                rx_str = rx_str.replace(rep[0],rep[1])
        if "M" not in rx_str.split(rx_arrow)[0]:
            rx_str = rx_str.split(rx_arrow)[0] + "+ M "+rx_arrow+rx_str.split(rx_arrow)[1]
        if "M" not in rx_str.split(rx_arrow)[1]:
            rx_str = rx_str.split(rx_arrow)[0]+rx_arrow+rx_str.split(rx_arrow)[1]+ "+ M "
        A, b, Ea = [float(a) for a in line.split(']')[1].split()[:3]]
        rx = {}
        rx['equation'] = rx_str
        rx['type'] = "three-body"
        rx['rate-constant'] = {"A": A, "b": b, "Ea": Ea}
        reactions.append(rx)
        
    photolysis = []
    # photolysis
    for i, line in enumerate(lines[type_ind['photolysis']:]):
        if len(line.strip()) == 0:
            break
        elif line[0] == '#':
            break
        rx_str = line.split(']')[0].split('[')[1].strip().replace('->','+ hv =>')
        for rep in to_replace:
            if rep[0] in rx_str:
                rx_str = rx_str.replace(rep[0],rep[1])
        rx = {}
        rx['equation'] = rx_str
        rx['type'] = "photolysis"
        photolysis.append(rx)
        
    sp_list = []
    for i,rx in enumerate(reactions):
        react, prod = [a.split('+') for a in rx['equation'].replace(' ','').split(rx_arrow)]
        sp_list+=react+prod
        
    for i,rx in enumerate(photolysis):
        react, prod = [a.split('+') for a in rx['equation'].replace(' ','').split('=>')]
        sp_list+=react+prod
    sp_list = list(set(sp_list))
    if "M" in sp_list:
        sp_list.remove('M')
    if "hv" in sp_list:
        sp_list.remove('hv')
        
    species = []
    atoms_list = []
    for sp in sp_list:
        entry = {}
        entry['name'] = sp
        entry['composition'] = composition[sp]
        atoms_list += [key for key in entry['composition'].keys()]
        species.append(entry)
        
    if vulcan_nasa9_data_folder != None:
        # We try to get thermodynamic data
        vulcan_sp_list = sp_list.copy()
        for rep in to_replace:
            if rep[1] in vulcan_sp_list:
                vulcan_sp_list[vulcan_sp_list.index(rep[1])] = rep[0]
        thermo = {} 
        for i,sp in enumerate(vulcan_sp_list):
            polys = np.loadtxt(vulcan_nasa9_data_folder+ '/'+ sp + '.txt')
            polys = polys.flatten()
            polys = [[float(a) for i,a in enumerate(polys[:10]) if i != 7], \
                     [float(a) for i,a in enumerate(polys[10:20]) if i != 7]]
            sp1 = sp_list[i]
            thermo[sp1] = {}
            thermo[sp1]['model'] = 'NASA9'
            thermo[sp1]['temperature-ranges'] = [0.0, 1000.0, 6000.0]
            thermo[sp1]['data'] = polys
            
        for i in range(len(species)):
            species[i]['thermo'] = thermo[species[i]['name']]    
        
    atoms_list = list(set(atoms_list))
    atoms = []
    for a in atoms_list:
        atoms.append(AtomData.atoms[a])
    
    # photolysis
    possible_photos = generate_photo_yaml_entries(sp_list)
    missing, not_missing = sort_photos(photolysis, possible_photos)
    
    # output
    out = {}
    if vulcan_nasa9_data_folder == None:
        out['reverse-reactions'] = False
    out['atoms'] = atoms
    out['species'] = species
    if photo_database == "Photochem":
        out['reactions'] = reactions + possible_photos
    elif photo_database == "Vulcan":
        out['reactions'] = reactions + not_missing
        out['missing'] = missing
    else:
        raise Exception('"Photochem" and "Vulcan" are the only options for photo_database')
    out = FormatReactions_main(out)
    
    if outfile == None:
        outfile = vulcan_rx_filename.replace('.txt','.yaml')
    
    # save
    fil = open(outfile,'w')
    yaml.dump(out,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
