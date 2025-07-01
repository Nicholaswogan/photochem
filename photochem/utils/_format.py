import yaml
from photochem_clima_data import DATA_DIR
import copy

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

class flowmap( dict ): pass
def flowmap_rep(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=True)

class blockseqtrue( list ): pass
def blockseqtrue_rep(dumper, data):
    return dumper.represent_sequence( u'tag:yaml.org,2002:seq', data, flow_style=True )

yaml.add_representer(blockseqtrue, blockseqtrue_rep)
yaml.add_representer(flowmap, flowmap_rep)

class MyDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()


def FormatReactions(filename, outfile):
    """Formats .yaml chemical reaction network

    Parameters
    ----------
    filename : str
        Path of input reaction network file.
    outfile : str
        Path of output formated reaction network file.

    """
    fil = open(filename,'r')
    data = yaml.load(fil,Loader=Loader)
    fil.close()
    out = FormatReactions_main(data)
    fil = open(outfile,'w')
    yaml.dump(out,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
    
def FormatReactions_main(data):
    
    order = ['reverse-reactions','atoms','species','particles','reactions','missing']
    copy = data.copy()
    data.clear()
    for key in order:
        if key in copy.keys():
            data[key] = copy[key]
            
    # Atmos
    if 'atoms' in data:
        for i in range(len(data['atoms'])):
            data['atoms'][i] = flowmap(data['atoms'][i])
    
    # Species
    if 'species' in data:
        for i in range(len(data['species'])):
            
            if data['species'][i]['name'] == False:
                data['species'][i]['name'] = "NO"
            
            order = ['name', 'composition', 'condensate', 'thermo', 'saturation','note']
            copy = data['species'][i].copy()
            data['species'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['species'][i][key] = copy[key]
                        
            data['species'][i]['composition'] = flowmap(data['species'][i]['composition'])
            if 'thermo' in data['species'][i].keys():
                
                order = ['model', 'reference-pressure','temperature-ranges','data']
                copy = data['species'][i]['thermo'].copy()
                data['species'][i]['thermo'].clear()
                for key in order:
                    if key in copy.keys():
                        data['species'][i]['thermo'][key] = copy[key]
                    
                data['species'][i]['thermo']['temperature-ranges'] = blockseqtrue(data['species'][i]['thermo']['temperature-ranges'])
                
                data['species'][i]['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(data['species'][i]['thermo']['data'])]

            if 'saturation' in data['species'][i].keys():
                flowstyle = ['parameters','vaporization','sublimation','super-critical']
                for key in flowstyle:
                    data['species'][i]['saturation'][key] = flowmap(data['species'][i]['saturation'][key])

    # Particles
    if 'particles' in data:
        for i in range(len(data['particles'])):
            data['particles'][i]['composition'] = flowmap(data['particles'][i]['composition'])
            if 'formation' not in data['particles'][i]:
                continue
            if data['particles'][i]['formation'] == 'reaction':
                flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
                for key in flowstyle:
                    if key in data['particles'][i].keys():
                        data['particles'][i][key] = flowmap(data['particles'][i][key])
            elif data['particles'][i]['formation'] == 'saturation':
                flowstyle = ['parameters','vaporization','sublimation','super-critical']
                for key in flowstyle:
                    data['particles'][i]['saturation'][key] = flowmap(data['particles'][i]['saturation'][key])
            
    # Reactions
    if 'reactions' in data:
        for i in range(len(data['reactions'])):
            order = ['equation','type','rate-constant','rate-constants','low-P-rate-constant',
                     'high-P-rate-constant','duplicate','efficiencies','JPL','citation']
            copy = data['reactions'][i].copy()
            data['reactions'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['reactions'][i][key] = copy[key]
                    
            flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
            for key in flowstyle:
                if key in data['reactions'][i].keys():
                    data['reactions'][i][key] = flowmap(data['reactions'][i][key])

            if 'rate-constants' in data['reactions'][i]:
                for j in range(len(data['reactions'][i]['rate-constants'])):
                    data['reactions'][i]['rate-constants'][j] = flowmap(data['reactions'][i]['rate-constants'][j])
                    
    return data
    
def FormatSettings(infile, outfile):
    """Formats a photochem settings file (e.g. settings_ModernEarth.yaml)

    Parameters
    ----------
    infile : str
        Path to input settings file.
    outfile : str
        Path of output formatted settings file.

    """
    fil = open(infile,'r')
    data = yaml.load(fil,Loader=Loader)
    fil.close()
    data = FormatSettings_main(data)
    fil = open(outfile,'w')
    yaml.dump(data,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
    
def FormatSettings_main(data):
    
    if 'planet' in data:
        if "rainout-species" in data['planet']['water'].keys():
            data['planet']['water']['rainout-species'] = blockseqtrue(data['planet']['water']['rainout-species'])

        if "condensation-rate" in data['planet']['water'].keys():
            data['planet']['water']['condensation-rate'] = flowmap(data['planet']['water']['condensation-rate'])
    
    if 'particles' in data:
        for i in range(len(data['particles'])):
            if "condensation-rate" in data['particles'][i]:
                data['particles'][i]["condensation-rate"] = \
                flowmap(data['particles'][i]["condensation-rate"])

    if 'boundary-conditions' in data:
        for i in range(len(data['boundary-conditions'])):
            if "lower-boundary" in data['boundary-conditions'][i]:
                order = ['type','vdep','mix','press','den','flux','height']
                copy = data['boundary-conditions'][i]['lower-boundary'].copy()
                data['boundary-conditions'][i]['lower-boundary'].clear()
                for key in order:
                    if key in copy.keys():
                        data['boundary-conditions'][i]['lower-boundary'][key] = copy[key]

                data['boundary-conditions'][i]['lower-boundary'] = flowmap(data['boundary-conditions'][i]['lower-boundary'])

                order = ['type','veff','flux']
                copy = data['boundary-conditions'][i]['upper-boundary'].copy()
                data['boundary-conditions'][i]['upper-boundary'].clear()
                for key in order:
                    if key in copy.keys():
                        data['boundary-conditions'][i]['upper-boundary'][key] = copy[key]

                data['boundary-conditions'][i]['upper-boundary'] = flowmap(data['boundary-conditions'][i]['upper-boundary'])
        
    return data

###
### Several routines for altering reaction files.
###

def species_in_reaction(rx):
    rx1 = rx.replace('<=>','=>').replace('(','').replace(')','')
    react, prod = [a.replace(' ','').split('+') for a in rx1.split('=>')]
    return react + prod

def mechanism_dict_with_atoms(dat_orig, atoms_names, exclude_species=[], remove_particles=False, remove_reaction_particles=False):

    dat = copy.deepcopy(dat_orig)

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

        if sp['name'] in exclude_species:
            exclude = True

        if not exclude:
            species.append(sp)

    if "particles" in dat and remove_particles:
        del dat['particles']
    if "particles" in dat:
        particles = []
        for i,sp in enumerate(dat['particles']):
            comp_tmp = [key for key in sp['composition'] if sp['composition'][key] > 0]

            exclude = False
            for tmp in comp_tmp:
                if tmp in exclude_atoms:
                    exclude = True
                    break

            if remove_reaction_particles and sp['formation'] == "reaction":
                exclude = True

            if sp['name'] in exclude_species:
                exclude = True
            if sp['formation'] == "saturation":
                if sp['gas-phase'] in exclude_species:
                    exclude = True
            elif not remove_reaction_particles and sp['formation'] == "reaction":
                rx_sp = species_in_reaction(sp['equation'])
                if any(a in exclude_species for a in rx_sp):
                    exclude = True
            
            if not exclude:
                particles.append(sp)

    if "reactions" in dat:
        reactions = []
        for i,rx in enumerate(dat['reactions']):
            sp = species_in_reaction(rx['equation'])
    
            exclude = False
            for s in sp:
                for tmp in comp[s]:
                    if tmp in exclude_atoms:
                        exclude = True
                        break
                if exclude:
                    break

            if any(a in exclude_species for a in sp):
                exclude = True

            if not exclude:
                reactions.append(rx)
                
    out = dat
    out['atoms'] = atoms
    out['species'] = species
    if 'particles' in dat:
        out['particles'] = particles
    if 'reactions' in dat:
        out['reactions'] = reactions

    return out

def resave_mechanism_with_atoms(
        infile, 
        outfile, 
        atoms_names, 
        exclude_species=[], 
        remove_particles=False, 
        remove_reaction_particles=False
    ):
    """Alters a reaction mechanism file with altered atoms.

    Parameters
    ----------
    infile : str
        Path to input yaml file
    outfile : str
        Path to output yaml file
    atoms_names : list
        List of atoms to include in the output file
    exclude_species : list, optional
        List of species to exclude, by default []
    remove_particles : bool, optional
        If True, then all particles are removed, by default False
    remove_reaction_particles : bool, optional
        If True, then partcles forming from reactions are removed, by default False
    """    

    with open(infile,'r') as f:
        dat = yaml.load(f,Loader=Loader)
    
    out = mechanism_dict_with_atoms(
        dat, 
        atoms_names, 
        exclude_species,
        remove_particles,
        remove_reaction_particles
    )

    out = FormatReactions_main(out)
    with open(outfile,'w') as f:
        yaml.dump(out,f,Dumper=MyDumper,sort_keys=False,width=70)

def generate_zahnle_earth_thermo(outfile='zahnle_earth_thermo.yaml', atoms_names=None, exclude_species=[], remove_particles=False):
    """Generates a thermodynamic file for equilibrium solving that includes
    condensible species (e.g., H2O condensate).

    Parameters
    ----------
    outfile : str, optional
        Name of the output file, by default 'zahnle_earth_thermo.yaml'
    atoms_names : list, optional
        List of atoms to keep. By default all atoms in the mechanism are kept
    exclude_species : list, optional
        List of species to exclude.
    remove_particles : bool, optional
        If True, then particles (i.e. condensates) will be removed, by default False.
    """    

    rx_folder = DATA_DIR+'/reaction_mechanisms/'

    with open(rx_folder+'zahnle_earth.yaml','r') as f:
        dat = yaml.load(f, Loader=Loader)

    with open(rx_folder+'condensate_thermo.yaml','r') as f:
        dat1 = yaml.load(f, Loader=Loader)

    # Delete information that is not needed
    for i,atom in enumerate(dat['atoms']):
        del dat['atoms'][i]['redox'] 
    del dat['particles']
    del dat['reactions']

    if not remove_particles:
        for i,sp in enumerate(dat1['species']):
            dat['species'].append(sp)

    if atoms_names is None:
        atoms_names = [a['name'] for a in dat['atoms']]
        
    dat = mechanism_dict_with_atoms(dat, atoms_names, exclude_species)

    dat = FormatReactions_main(dat)

    with open(outfile, 'w') as f:
        yaml.dump(dat,f,Dumper=MyDumper,sort_keys=False,width=70)

def zahnle_rx_and_thermo_files(
        atoms_names=['H','He','N','O','C','S'], 
        rxns_filename='photochem_rxns.yaml', 
        thermo_filename='photochem_thermo.yaml',
        exclude_species=[],
        remove_particles=False, 
        remove_reaction_particles=False
    ):
    """Generates input reactions and thermodynamic files for photochem by altering the main
    reaction network (zahnle_earth.yaml).

    Parameters
    ----------
    atoms_names : list, optional
        Atoms to include in the thermodynamics, by default ['H','He','N','O','C','S']
    rxns_filename : str, optional
        Name of output reactions file, by default 'photochem_rxns.yaml'
    thermo_filename : str, optional
        Name of output thermodynamic file, by default 'photochem_thermo.yaml'
    exclude_species : list, optional
        List of species to exclude.
    remove_particles : bool, optional
        If True, then particles will be removed, by default False.
    remove_reaction_particles : bool, optional
        If True, then reactions particles are removed, by default False.
    """

    zahnle_earth = DATA_DIR+'/reaction_mechanisms/zahnle_earth.yaml'
    # Kinetics
    if rxns_filename is not None:
        resave_mechanism_with_atoms(zahnle_earth, rxns_filename, atoms_names, exclude_species, remove_particles, remove_reaction_particles)

    # Thermodynamics
    if thermo_filename is not None:
        generate_zahnle_earth_thermo(thermo_filename, atoms_names, exclude_species, remove_particles)
