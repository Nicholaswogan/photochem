from ._photochem import photochem, photochem_vars, photochem_data, photochem_setup, photochem_const
import numpy as np
import os

from .utils import get_character_arr

string_format = "{:1024}"
rootdir = os.path.dirname(os.path.realpath(__file__))+'/'
photochem_vars.data_dir = string_format.format(rootdir+"../data")

class PhotoException(Exception):
    pass

class Photochem():
    
    def __init__(self, mechanism_file, settings_file, flux_file, atmosphere_txt):
        self.data = photochem_data
        self.vars = photochem_vars
        self.const = photochem_const
        err = photochem_setup.setup(mechanism_file,\
                                    settings_file,\
                                    flux_file,\
                                    atmosphere_txt)
        if len(err.strip()) > 0:
            raise PhotoException(err.decode("utf-8").strip())
            
    def photochemical_equilibrium(self, maxsteps = 10000, rtol = 1e-3, atol = 1e-25):
        success, err = photochem.photo_equilibrium(maxsteps, rtol, atol)
        if len(err.strip()) > 0:
            raise PhotoException(err.decode("utf-8").strip())
        return success
    
    def out2in(self):
        err = photochem_setup.out2in()
        if len(err.strip()) > 0:
            raise PhotoException(err.decode("utf-8").strip())
    
    def out2atmosphere_txt(self, filename, overwrite = False, clip = True):
        err = photochem_setup.out2atmosphere_txt(filename, overwrite, clip)
        if len(err.strip()) > 0:
            raise PhotoException(err.decode("utf-8").strip())
    
    def out_dict(self):
        if not self.vars.at_photo_equilibrium:
            raise PhotoException("Must integrate to equilibrium before making a dictionary")
        out = {}
        out['alt'] = self.vars.z/1e5
        names = self.species_names
        for i in range(self.data.nq):
            out[names[i]] = self.vars.usol_out[i,:]
        return out
        
    @property
    def species_names(self):
        names = get_character_arr(photochem_setup.get_species_names, \
                                  self.data.nsp + 2, 15)
        return names
        
    @property
    def atoms_names(self):
        names = get_character_arr(photochem_setup.get_atoms_names, \
                                  self.data.natoms, 8)
        return names
        