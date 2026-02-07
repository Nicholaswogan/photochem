from photochem_clima_data import DATA_DIR
import os
from ._format import yaml, Loader, MyDumper, FormatReactions_main, flowmap

def species_dict_for_climate(species, condensates, particles=None):
    """Generates a species dictionary for climate simulations with AdiabatClimate.

    Parameters
    ----------
    species : list
        List of species to include in the file.
    condensates : list
        List of species that should be condensates.
    particles : dict, optional
        List of particles with their composition. For example, 
        `particles = [{'name': 'HCaer', 'composition': {'C': 6, 'H': 2}}`, 
        by default {}.
    
    Returns
    -------
    species_file : dict
        Dictionary of the species file.
    """

    if not set(condensates).issubset(species):
        raise ValueError("`condensates` must be a subset of `species`.")

    zahnle_earth = os.path.join(DATA_DIR,'reaction_mechanisms/zahnle_earth.yaml')

    with open(zahnle_earth,'r') as f:
        dat = yaml.load(f, Loader=Loader)

    # Get saturation information
    saturation = {}
    for sp in dat['particles']:
        if 'gas-phase' not in sp:
            continue
        saturation[sp['gas-phase']] = sp['saturation']

    # Get species
    species_new = []
    for sp in dat['species']:
        if sp['name'] in species:
            if sp['name'] in condensates:
                sp['saturation'] = saturation[sp['name']]
            species_new.append(sp)

    # Get atoms
    atoms = []
    for sp in species_new:
        atoms += [key for key in sp['composition'] if sp['composition'][key] > 0]
    atoms = list(set(atoms))

    atoms_new = []
    for a in dat['atoms']:
        if a['name'] in atoms:
            a.pop('redox')
            atoms_new.append(a)

    species_file = {
        'atoms': atoms_new,
        'species': species_new,
    }
    if particles is not None:
        species_file['particles'] = particles

    return species_file
    
def species_file_for_climate(filename, species, condensates, particles=None):
    """Generates a species dict for climate simulations with AdiabatClimate.

    Parameters
    ----------
    filename : str
        Path of the output yaml file.
    species : list
        List of species to include in the file.
    condensates : list
        List of species that should be condensates.
    particles : dict, optional
        List of particles with their composition. For example, 
        `particles = [{'name': 'HCaer', 'composition': {'C': 6, 'H': 2}}`, 
        by default {}.
    
    Returns
    -------
    species_file : dict
        Dictionary of the species file.
    """

    species_file = species_dict_for_climate(species, condensates, particles)

    species_file = FormatReactions_main(species_file)
    with open(filename, 'w') as f:
        yaml.dump(species_file,f,Dumper=MyDumper,sort_keys=False,width=70)

def settings_dict_for_climate(planet_mass, planet_radius, surface_albedo, 
                              number_of_layers=50, number_of_zenith_angles=4, 
                              photon_scale_factor=1.0, opacities=None):
    """Generates a settings dictionary for the AdiabatClimate model.

    Parameters
    ----------
    planet_mass : float
        Planet mass in grams.
    planet_radius : float
        Planet radius in cm.
    surface_albedo : float
        Surface albedo of the planet.
    number_of_layers : int, optional
        Number of atmospheric layers, by default 50.
    number_of_zenith_angles : int, optional
        Number of solar zenith angles to do solar radiative transfer at, 
        by default 4.
    photon_scale_factor : float, optional
        A value to multiply the stellar flux by, by default 1.0.
    opacities: dict, optional
        A dictionary describing the opacities, by default None.
    
    Returns
    -------
    settings : dict
        Dictionary of the settings file.
    """    
    
    if opacities is None:
        opacities = {
            'k-distributions': True, 
            'CIA': True, 
            'rayleigh': True, 
            'photolysis-xs': True, 
            'water-continuum': 'MT_CKD'
        }

    settings = {
        'atmosphere-grid': {
            'number-of-layers': number_of_layers
        },
        'planet': {
            'planet-mass': planet_mass, 
            'planet-radius': planet_radius, 
            'surface-albedo': surface_albedo,
            'number-of-zenith-angles': number_of_zenith_angles,
            'photon-scale-factor': photon_scale_factor
        },
        'optical-properties': {
            'k-method': 'RandomOverlapResortRebin',
            'opacities': opacities
        }
    }
    
    return settings

def settings_file_for_climate(filename, planet_mass, planet_radius, surface_albedo, 
                              number_of_layers=50, number_of_zenith_angles=4, 
                              photon_scale_factor=1.0, opacities=None):
    """Generates a settings file for the AdiabatClimate model.

    Parameters
    ----------
    filename : str
        Path to output file.
    planet_mass : float
        Planet mass in grams.
    planet_radius : float
        Planet radius in cm.
    surface_albedo : float
        Surface albedo of the planet.
    number_of_layers : int, optional
        Number of atmospheric layers, by default 50.
    number_of_zenith_angles : int, optional
        Number of solar zenith angles to do solar radiative transfer at, 
        by default 4.
    photon_scale_factor : float, optional
        A value to multiply the stellar flux by, by default 1.0.
    opacities: dict, optional
        A dictionary describing the opacities, by default None.
    """
    settings = settings_dict_for_climate(
        planet_mass, 
        planet_radius, 
        surface_albedo, 
        number_of_layers, 
        number_of_zenith_angles, 
        photon_scale_factor,
        opacities
    )
    
    settings['optical-properties']['opacities'] = flowmap(settings['optical-properties']['opacities'])
    with open(filename, 'w') as f:
        yaml.dump(settings,f,Dumper=MyDumper,sort_keys=False,width=70)