import numpy as np
import os
import copy
import h5py
import shutil

from photochem_clima_data import DATA_DIR
from ._format import FormatReactions_main, MyDumper, Loader, yaml
from ._convert_utils import AtomData, compare2reactions, \
                            generate_photo_yaml_entries, sort_photos

def read_all_compose(all_compose_filename):
    "Retrieves the composition of all species in VULCAN"

    with open(all_compose_filename,'r') as f:
        lines = f.readlines()
    
    atoms = lines[0].split()[1:17]
    composition = {}

    for line in lines[1:]:
        sp = line.split()[0]
        composition[sp] = {}
        na = line.split()[1:17]
        for i,a in enumerate(na):
            if int(a) > 0:
                composition[sp][atoms[i]] = int(a)

    return composition

def get_zahnle_info():
    """Gets condensates and reactions from Zahnle, in case they
    are needed"""

    zahnle_earth = DATA_DIR+'/reaction_mechanisms/zahnle_earth.yaml'
    with open(zahnle_earth,'r') as f:
        dat = yaml.load(f,yaml.Loader)
    zahnle_condensates = {}
    for i,p in enumerate(dat['particles']):
        if p['formation'] == 'saturation':
            zahnle_condensates[p['gas-phase']] = p
    zahnle_reactions = {}
    for i,rx in enumerate(dat['reactions']):
        zahnle_reactions[rx['equation']] = rx

    return zahnle_condensates, zahnle_reactions

def check_xs_and_qy(out):
    """Checks cross sections and quantum yields to make sure they
    are valid."""

    # Must be sorted
    wv = out['wv'].astype(np.float32)
    assert np.all(wv[:-1] < wv[1:])

    if 'ratios' in out:
        wv = out['ratios']['wv'].astype(np.float32)
        assert np.all(wv[:-1] < wv[1:])

def make_h5_from_dict(species, out, outdir):
    "Creates an hdf5 photolysis cross section file"

    check_xs_and_qy(out)
    
    with h5py.File(outdir+species+'.h5','w') as f:
        dset = f.create_dataset("wavelengths", out['wv'].shape, 'f')
        dset[:] = out['wv']

        dset = f.create_dataset("photoabsorption", out['xsa'].shape, 'f')
        dset[:] = out['xsa']

        dset = f.create_dataset("photodissociation", out['xsp'].shape, 'f')
        dset[:] = out['xsp']

        dset = f.create_dataset("photoionisation", out['xsi'].shape, 'f')
        dset[:] = out['xsi']

        if 'ratios' in out:
            grp = f.create_group("photodissociation-qy")

            dset = grp.create_dataset("wavelengths", out['ratios']['wv'].shape, 'f')
            dset[:] = out['ratios']['wv']
            
            for ratio in out['ratios']:
                if ratio != 'wv':
                    dset = grp.create_dataset(ratio, out['ratios'][ratio].shape, 'f')
                    dset[:] = out['ratios'][ratio]  

def create_supporting_data(xs_info, vulcan_xs_folder, data_dir):
    """Creates a data folder with new VULCAN photolysis cross
    sections."""

    if data_dir is None:
        data_dir = DATA_DIR

    if os.path.isdir('vulcandata'):
        shutil.rmtree('vulcandata')
    _ = shutil.copytree(data_dir, 'vulcandata')
    _ = shutil.move('vulcandata/xsections/bins.h5','vulcandata/')
    shutil.rmtree('vulcandata/xsections')
    os.mkdir('vulcandata/xsections')
    _ = shutil.move('vulcandata/bins.h5','vulcandata/xsections/')

    photospecies = [a['sp'] for a in xs_info]
    photospecies = list(set(photospecies))
    branches = {}
    for sp in photospecies:
        tmp = []
        for a in xs_info:
            if sp == a['sp']:
                tmp.append({'equation': a['equation'],'ind': a['ind']})
        branches[sp] = tmp

    for sp in branches:
        xsfile = vulcan_xs_folder+'/'+sp+'/'+sp+'_cross.csv'
        wv, xsa, xsp, xsi = np.loadtxt(xsfile,skiprows=1,delimiter=',').T
        _, inds = np.unique(wv,return_index=True)
        wv = wv[inds]
        xsa = xsa[inds]
        xsp = xsp[inds]
        xsi = xsi[inds]

        qyfile = vulcan_xs_folder+'/'+sp+'/'+sp+'_branch.csv'
        with open(qyfile,'r') as f:
            lines = f.readlines()

        labels = lines[1][1:].strip().split(',')
        labels = [a.strip() for a in labels]
        for i in range(len(labels)-1):
            assert labels[i+1] == 'br_ratio_'+str(i+1)

        # Get the data
        dat = np.array([[float(b) for b in a.strip().split(',')] for a in lines[2:]])

        _, inds = np.unique(dat[:,0],return_index=True)
        ratios = {'wv': dat[inds,0]}
        for i,rx in enumerate(branches[sp]):
            assert rx['ind'] == i + 1
            ratios[rx['equation']] = dat[inds,i+1]

        out = {
            'wv': wv,
            'xsa': xsa,
            'xsp': xsp,
            'xsi': xsi,
            'ratios': ratios
        }
        # print(sp)
        make_h5_from_dict(sp, out, 'vulcandata/xsections/')

def vulcan2yaml(vulcan_rx_filename, thermo_folder, data_dir=None):
    """Converts Vulcan reactions and cross sections to a format that
    works with Photochem. Upon return, the routine will have saved a
    yaml file with a similar name to the input `vulcan_rx_filename`, and
    also a directory "vulcandata", which is a copy of the Photochem data, 
    except with the Vulcan photolysis cross sections.

    Parameters
    ----------
    vulcan_rx_filename : str
        Path to Vulcan input reactions file.
    thermo_folder : str
        Path to Vulcan "thermo" folder
    data_dir : str
        Path to the Photochem data folder. If `None`, then the data shipped
        with Photochem is used.
    """
    
    # Path to the "thermo" folder
    folder = thermo_folder+'/'
    # Composition file
    all_compose_filename = folder+'all_compose.txt'
    # NASA9 folder
    vulcan_nasa9_data_folder = folder+'NASA9'
    # Cross section folder
    vulcan_xs_folder = folder+'photo_cross'

    # Composition of all species
    composition = read_all_compose(all_compose_filename)

    # Get Zahnle condensates and reactions
    zahnle_condensates, zahnle_reactions = get_zahnle_info()
    
    rx_arrow = "<=>"
        
    # parse reactions
    with open(vulcan_rx_filename,'r') as f:
        lines = f.readlines()

    type_ind = {}
    for i, line in enumerate(lines):
        if line.startswith("# Two-body"):
            type_ind['elementary'] = i + 2
        elif line.startswith("# 3-body and Diss"):
            type_ind['falloff'] = i + 2
        elif line.startswith("# 3-body reactions without"):
            type_ind['three-body'] = i + 2
        elif line.startswith("# photo"):
            type_ind['photolysis'] = i + 2
        elif line.startswith("# condensation"):
            type_ind['condensation'] = i + 1
        elif line.startswith("# special cases"):
            type_ind['special'] = i + 2

    reactions = []

    # elementary
    if 'elementary' in type_ind:
        for i, line in enumerate(lines[type_ind['elementary']:]):
            if len(line.strip()) == 0:
                continue
            elif line[0] == '#':
                break
            rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
            A, b, Ea = [float(a) for a in line.split(']')[1].split()[:3]]
            rx = {}
            rx['equation'] = rx_str
            rx['rate-constant'] = {"A": A, "b": b, "Ea": Ea}
            reactions.append(rx)
        
    # falloff
    if 'falloff' in type_ind:
        for i, line in enumerate(lines[type_ind['falloff']:]):
            if len(line.strip()) == 0:
                continue
            elif line[0] == '#':
                break
            rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
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
    if 'three-body' in type_ind:
        for i, line in enumerate(lines[type_ind['three-body']:]):
            if len(line.strip()) == 0:
                continue
            elif line[0] == '#':
                break
            rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
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

    # Special
    if 'special' in type_ind:
        for i, line in enumerate(lines[type_ind['special']:]):
            if len(line.strip()) == 0:
                continue
            elif line[0] == '#':
                break
            rx_str = line.split(']')[0].split('[')[1].strip().replace('->',rx_arrow)
            found = False
            for key in zahnle_reactions:
                if compare2reactions(rx_str, key):
                    reactions.append(zahnle_reactions[key])
                    found = True
                    break
            if found:
                print('Using zahnle_earth.yaml rate for special reaction '+rx_str)
            else:
                print('Ignoring special reaction '+rx_str)
                
    # Condensation
    particles = []
    if 'condensation' in type_ind:
        for i, line in enumerate(lines[type_ind['condensation']:]):
            if len(line.strip()) == 0:
                continue
            elif line[0] == '#':
                break
            sp = line.split(']')[0].split('[')[1].split('->')[0].strip()
            if sp in zahnle_condensates:
                particles.append(zahnle_condensates[sp])
            else:
                print('Skipping condensate '+sp+' because we do not have data for it.')
        
    # photolysis
    photolysis = []
    xs_info = []
    for i, line in enumerate(lines[type_ind['photolysis']:]):
        if len(line.strip()) == 0:
            continue
        elif line[0] == '#':
            break
        rx_str = line.split(']')[0].split('[')[1].strip().replace('->','+ hv =>')
        react, prod = [sorted(a.replace(' ','').split('+')) for a in rx_str.split('=>')]
        rx_str = ' + '.join(react)+' => '+' + '.join(prod)
        rx = {}
        rx['equation'] = rx_str
        rx['type'] = "photolysis"
        photolysis.append(rx)
        rx = copy.deepcopy(rx)
        rx['sp'] = line.split(']')[1].split()[0].strip()
        rx['ind'] = int(line.split(']')[1].split()[1].strip())
        xs_info.append(rx)

    # List of species
    sp_list = []
    for i,rx in enumerate(reactions):
        react, prod = [a.split('+') for a in rx['equation'].replace(' ','').replace('(','').replace(')','').split(rx_arrow)]
        sp_list+=react+prod

    for i,rx in enumerate(photolysis):
        react, prod = [a.split('+') for a in rx['equation'].replace(' ','').split('=>')]
        sp_list+=react+prod
    sp_list = list(set(sp_list))
    if "M" in sp_list:
        sp_list.remove('M')
    if "hv" in sp_list:
        sp_list.remove('hv')

    # Get composition for species
    species = []
    atoms_list = []
    for sp in sp_list:
        entry = {}
        entry['name'] = sp
        entry['composition'] = composition[sp]
        atoms_list += [key for key in entry['composition'].keys()]
        species.append(entry)

    # Thermodynamics for species
    thermo = {} 
    for i,sp in enumerate(sp_list):
        polys = np.loadtxt(vulcan_nasa9_data_folder+ '/'+ sp + '.txt')
        polys = polys.flatten()
        polys = [[float(a) for i,a in enumerate(polys[:10]) if i != 7], \
                [float(a) for i,a in enumerate(polys[10:20]) if i != 7]]
        thermo[sp] = {}
        thermo[sp]['model'] = 'NASA9'
        thermo[sp]['temperature-ranges'] = [0.0, 1000.0, 6000.0]
        thermo[sp]['data'] = polys
    for i in range(len(species)):
        species[i]['thermo'] = thermo[species[i]['name']]    
        
    atoms_list = list(set(atoms_list))
    atoms = []
    for a in atoms_list:
        atoms.append(AtomData.atoms[a])

    out = {}
    out['atoms'] = atoms
    out['species'] = species
    if len(particles) > 0:
        out['particles'] = particles
    out['reactions'] = reactions + photolysis

    out = FormatReactions_main(out)

    outfile = os.path.basename(vulcan_rx_filename).replace('.txt','.yaml')

    with open(outfile,'w') as f:
        yaml.dump(out,f,Dumper=MyDumper,sort_keys=False,width=70)

    # Create the supporting data folder
    create_supporting_data(xs_info, vulcan_xs_folder, data_dir)