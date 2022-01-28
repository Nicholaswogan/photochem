import numpy as np
import os
from ruamel.yaml import YAML
yaml = YAML(typ='safe')

def check_consistency(verbose=True):
    '''Checks for consistency between metadata.yaml, and the data folders and files
    '''
    fil = open('metadata.yaml','r')
    meta_data = yaml.load(fil)
    fil.close()
    all_species = []
    for species in meta_data.keys():
        if species!='overall-notes':
            all_species.append(species)
            # check directory exists
            if not os.path.isdir(species):
                raise Exception(species+' in xsections meta-data is not a directory.')

            if not os.path.isfile(species+'/'+species+'_xs.txt'):
                raise Exception(species+' in xsections meta-data does not have corresponding cross section data.')

            for rx in meta_data[species]['reactions'].keys():
                if not os.path.isfile(species+'/'+rx.replace(' ','_')+'.txt'):
                    raise Exception('Reaction '+rx+' in xsections meta-data does not have a corresponding data file.')

    directories = [a for a in os.listdir() if os.path.isdir(a) and a[0]!='.']
    directories2remove = ['__pycache__', "Zahnle_Kevin_data"]
    for directory in directories2remove:
        try:
            directories.remove(directory)
        except:
            pass
    
    # Make sure meta data file has all the reactions in the folder
    for directory in directories:
        reactions = []
        for key in meta_data[directory]['reactions'].keys():
            reactions.append(key)
        files = os.listdir(directory)
        for file in files:
            if 'xs' not in file and file[0]!='.':
                reaction = ' '.join(file.replace('.txt','').split('_'))
                if reaction not in reactions:
                    raise Exception("metadata.yaml does not have an entry for reaction "+ \
                                    reaction+" but data for this reaction exists")

    if set(all_species) != set(directories):
        bad = list(set(all_species).symmetric_difference(set(directories)))
        bad_str = ''
        for b in bad:
            bad_str += ' "'+b+'"'
        raise Exception('The folder(s)'+bad_str+' are not in the meta-data file.')

    # now read in data, and make sure it is sensible
    for spec in all_species:
        # get the xs
        fil = open(spec+'/'+spec+'_xs.txt','r')
        lines = fil.readlines()
        fil.close()

        ntemps = len(lines[1].split())-1
        wv_xs = []
        xs = []

        for i,line in enumerate(lines[2:]):
            wv_xs.append(float(line.split()[0]))
            tmp = []
            for j in range(1,ntemps+1):
                tmp.append(float(line.split()[j]))
            xs.append(np.array(tmp))

        xs_temps = [float(a.strip('K')) for a in lines[1].strip().split()[1:]]
        for i in range(0,len(xs_temps)-1):
            if not xs_temps[i] < xs_temps[i+1]:
                raise Exception('Temperature columns must be sorted for species '+sp)


        wv_qys = []
        qys = []
        qy_temps = []
        for rx in meta_data[spec]['reactions']:
            rxname = '_'.join(rx.split())
            fil = open(spec+'/'+rxname+'.txt','r')
            lines = fil.readlines()
            fil.close()

            ntemps = len(lines[1].split())-1
            wv_qy = []
            qy = []
            for line in lines[2:]:
                tmp = []
                wv_qy.append(float(line.split()[0]))
                for j in range(1,ntemps+1):
                    tmp.append(float(line.split()[j]))
                qy.append(tmp)
            wv_qys.append(wv_qy)
            qys.append(qy)

            qy_temps.append([float(a.strip('K')) for a in lines[1].strip().split()[1:]])
            for i in range(0,len(qy_temps[-1])-1):
                if not qy_temps[-1][i] < qy_temps[-1][i+1]:
                    raise Exception('Temperature columns must be sorted for reaction '+rx+' quantum yields')

        nr = len(meta_data[spec]['reactions'])
        for i in range(1,nr):
            if len(qys[0][0]) != len(qys[i][0]):
                raise Exception('The quantum yields for '+spec+' must all have the same number of temperature columns, with the same corresponding temperature')

        for i in range(1,nr):
            for j in range(len(qy_temps[0])):
                if qy_temps[0][j] != qy_temps[i][j]:
                    raise Exception('The quantum yields for '+spec+' have data at different temperatures.\n'+\
                                    'This is not allowed.')

        for i in range(1,nr):
            if not np.isclose(np.min(wv_qys[0]),np.min(wv_qys[i])):
                raise Exception('Wavelengths for the quantum yields for ',spec,' have different lower bounds. They must be the same.')
            if not np.isclose(np.max(wv_qys[0]),np.max(wv_qys[i])):
                raise Exception('Wavelengths for the quantum yields for ',spec,' have different upper bounds. They must be the same.')
        if not np.isclose(np.min(wv_qys[0]),np.min(wv_xs)):
            raise Exception('Wavelengths for the quantum yields and cross sections for ',spec,' have different lower bounds. They must be the same.')
        if not np.isclose(np.max(wv_qys[0]),np.max(wv_xs)):
            raise Exception('Wavelengths for the quantum yields and cross sections for ',spec,' have different upper bounds. They must be the same.')

        # check qys add to less than 1
        wv_new = np.linspace(np.min(wv_xs),np.max(wv_xs),1000)
        qy_tot = np.zeros(1000)

        for i in range(len(qy_temps[0])):
            qy_tot = np.zeros(1000)
            for j in range(len(qy_temps)):
                qy_tot += np.interp(wv_new, wv_qys[j], np.array(qys[j])[:,i])

        if not all(qy_tot<1.01):
            raise Exception('Quantum yields sum for ',spec,' sum to >1. They should sum to <= 1.')

        # check meta-data wv citations align
        for rx in meta_data[spec]['reactions']:
            low = []
            high = []
            for citation in meta_data[spec]['reactions'][rx]['citations']:
                tmp = [float(a) for a in citation['nm-range'].split('-')]
                low.append(tmp[0])
                high.append(tmp[1])
            if not np.isclose(np.min(low),np.min(wv_qys[0])) or not np.isclose(np.max(high),np.max(wv_qys[0])):
                raise Exception('The citation range in the meta-data does not align with the data for reaction '+rx)
        low = []
        high = []
        for citation in meta_data[spec]['xsections']['citations']:
            tmp = [float(a) for a in citation['nm-range'].split('-')]
            low.append(tmp[0])
            high.append(tmp[1])
        if not np.isclose(np.min(low),np.min(wv_xs)) or not np.isclose(np.max(high),np.max(wv_xs)):
                raise Exception('The citation range in the meta-data does not align with the xs data for '+spec)

    if verbose:
        print('metadata.yaml is consistent with the data files.')

def generate_photo_yaml_entries(species_list):
    '''Generates list of photolysis reactions given a list of species.
    '''
    species_set = set(species_list)
    fil = open('metadata.yaml','r')
    meta_data = yaml.load(fil)
    fil.close()
    all_photo_species = [key for key in meta_data.keys() if key != 'overall-notes']
    photo_species = list(species_set.intersection(set(all_photo_species)))

    # grab reactions
    rx_list = []
    for species in photo_species:
        for rx in meta_data[species]['reactions'].keys():
            tmp = set([a.strip() for a in rx.split('=>')[1].split('+')])
            if tmp.issubset(species_set):
                rx_list.append({})
                rx_list[-1]['equation'] = rx
                rx_list[-1]['type'] = 'photolysis'
    return rx_list

if __name__ == "__main__":
    check_consistency()
