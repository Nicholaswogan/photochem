from ._photochem import photochem, photochem_vars, photochem_setup, photochem_const
import os
string_format = "{:1024}"

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'
photochem_vars.data_dir = string_format.format(rootdir+"../data")