from ._equilibrate import ChemEquiAnalysis, EquilibrateException
from ._equilibrate import __version__
from .utils._format import FormatReactions_main, yaml, Loader, MyDumper, mechanism_dict_with_atoms
from photochem_clima_data import DATA_DIR

def generate_zahnle_earth_thermo(outfile='zahnle_earth_thermo.yaml', atoms_names=None, exclude_species=[]):
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

    for i,sp in enumerate(dat1['species']):
        dat['species'].append(sp)

    if atoms_names is None:
        atoms_names = [a['name'] for a in dat['atoms']]
        
    dat = mechanism_dict_with_atoms(dat, atoms_names, exclude_species)

    dat = FormatReactions_main(dat)

    with open(outfile, 'w') as f:
        yaml.dump(dat,f,Dumper=MyDumper,sort_keys=False,width=70)
