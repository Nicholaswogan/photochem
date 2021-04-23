import os
from ruamel.yaml import YAML
yaml = YAML(typ='safe')

def check_consistency(verbose=True):
    '''Checks for consistency between xs_metadata.yaml, and the data folders/files
    '''
    fil = open('xs_metadata.yaml','r')
    meta_data = yaml.load(fil)
    fil.close()
    for species in meta_data.keys():
        # check directory exists
        if not os.path.isdir(species):
            raise Exception(species+' in xsections meta-data is not a directory.')

        if not os.path.isfile(species+'/'+species+'_xs.txt'):
            raise Exception(species+' in xsections meta-data does not have corresponding cross section data.')

        for rx in meta_data[species]['reactions'].keys():
            if not os.path.isfile(species+'/'+rx.replace(' ','_')+'.txt'):
                raise Exception('Reaction '+rx+' in xsections meta-data does not have a corresponding data file.')
    if verbose:
        print('xs_metadata.yaml is consistent with the data files.')

def generate_photo_yaml_entries(species_list):
    '''Generates list of photolysis reactions given a list of species.
    '''
    check_consistency(verbose=False)
    species_set = set(species_list)
    fil = open('xs_metadata.yaml','r')
    meta_data = yaml.load(fil)
    fil.close()
    all_photo_species = [key for key in meta_data.keys()]
    photo_species = list(species_set.intersection(set(all_photo_species)))

    # grab reactions
    rx_list = []
    for species in photo_species:
        for rx in meta_data[species]['reactions'].keys():
            tmp = set([a.strip() for a in rx.split('=>')[1].split('+')])
            if tmp.issubset(species_set):        
                rx_list.append({})
                rx_list[-1]['reaction'] = rx
                rx_list[-1]['type'] = 'photolysis'
    return rx_list
    
if __name__ == "__main__":
    check_consistency()