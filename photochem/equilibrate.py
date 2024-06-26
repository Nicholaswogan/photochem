from ._equilibrate import ChemEquiAnalysis, EquilibrateException
from ._equilibrate import __version__
from .utils._format import FormatReactions_main, yaml, Loader, MyDumper
import os

def generate_zahnle_earth_thermo(outfile='zahnle_earth_thermo.yaml'):
    """Generates a thermodynamic file for equilibrium solving that includes
    condensible species (e.g., H2O condensate).

    Parameters
    ----------
    outfile : str, optional
        Name of the output file, by default 'zahnle_earth_thermo.yaml'
    """    

    rx_folder = os.path.dirname(os.path.realpath(__file__))+'/data/reaction_mechanisms/'

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

    dat = FormatReactions_main(dat)

    with open(outfile, 'w') as f:
        yaml.dump(dat,f,Dumper=MyDumper,sort_keys=False,width=70)
