from ._format import FormatReactions_main, FormatSettings_main, MyDumper, yaml
from ._convert_utils import generate_photo_yaml_entries, sort_photos
    
def atmos2yaml(rx_file, species_file, outfile, photo_database = "Photochem", with_citations = False):
    """Converts Atmos reactions to .yaml format compatable with Photochem

    Parameters
    ----------
    rx_file : str
        Path to Atmos reactions file (e.g. reactions.rx)
    species_file : str
        Path to Atmos species file (e.g. species.dat)
    outfile : str
        Name of output .yaml file in photochem format
    photo_database : str
        Options are "Photochem" or "Atmos". If "Photochem", then will use all
        possible photolysis reactions given the cross section data contained
        within Photochem. If "Atmos", then will photolysis reactions in
        rx_file that are possible, given Photochem cross section data.
        

    """
    data = rx2dict(rx_file, with_citations = with_citations)
    species, particles = get_species(species_file)
    out = make_rx_yaml(species, particles, data, photo_database = photo_database)
    out = FormatReactions_main(out)
    fil = open(outfile,'w')
    yaml.dump(out,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()

def atmosbc2yaml(species_file, outfile, short_lived = False):
    """Converts Atmos boundary conditions into .yaml format compatible with Photochem.

    Parameters
    ----------
    species_file : str
        Path to Atmos species file which contains boundary conditions (e.g. species.dat)
    outfile : _type_
        Name of output .yaml file which will contain boundary conditions
    """
    out = atmosbc2yaml_main(species_file, short_lived = short_lived)
    with open(outfile,'w') as f:
        f.write(out)

def atmosbc2yaml_main(species_file, short_lived = False):

    with open(species_file) as f:
        lines = f.readlines()

    bc_list = []
    for i,line in enumerate(lines):
        if not line.startswith("*"):
            tmp = line.split()
            spname = tmp[0]

            if tmp[1] == "LL":
                lbound = int(tmp[8])
                vdep = float(tmp[9])
                mix = float(tmp[10])
                sgflux = float(tmp[11])
                disth = float(tmp[12])
                mbound = int(tmp[13])
                smflux = float(tmp[14])
                veff = float(tmp[15])

                if lbound == 0 and vdep == 0.0 \
                    and mbound == 0 and veff == 0.0:
                    # skip the boundary condition
                    # because it is default.
                    continue
                
                if lbound == 0:
                    # vdep
                    lb = {"type": "vdep", "vdep": vdep}
                elif lbound == 1:
                    lb = {"type": "mix", "mix": mix}
                elif lbound == 2:
                    lb = {"type": "flux", "flux": sgflux}
                elif lbound == 3:
                    lb = {"type": "vdep + dist flux", "vdep": vdep, "flux": sgflux, "height": disth}
                else:
                    raise Exception('"'+species_file+'" contains an invalid lower boundary condition')

                if mbound == 0:
                    ub = {"type": "veff", "veff": veff}
                elif mbound == 1:
                    raise Exception('"'+species_file+\
                        '" contains a fixed mixing ratio upper boundary '+\
                        'condition which is not compatible with Photochem')
                elif mbound == 2:
                    ub = {"type": "flux", "flux": smflux}
                
                entry = {}
                entry['name'] = spname
                entry['lower-boundary'] = lb
                entry['upper-boundary'] = ub
                bc_list.append(entry)

            elif tmp[1] == "SL" and short_lived:

                entry = {}
                entry['name'] = spname
                entry['type'] = 'short lived'
                bc_list.append(entry)

    out = {}
    out['boundary-conditions'] = bc_list
    out = FormatSettings_main(out)
    return yaml.dump(out, Dumper=MyDumper ,sort_keys=False, width=70)
  
def rx2dict(filename, with_citations = False):

    fil = open(filename,'r')
    lines = fil.readlines()
    fil.close()

    data = {}
    data['reactions'] = []
    data['photolysis'] = []
    data['weird'] = []

    for line in lines:
        
        if line.startswith("REACTANTS"):
            continue
        
        if "HCAER" in line or "S8AER" in line:
            continue
            
        reaction = {}
        weird = False
        
        if '!' in line:
            ind = line.index('!')
            citation = line[ind+1:].strip() 
        else:
            citation = None
            
        rx_type = line[48:53]
        if " " in rx_type:
            # UGH!
            rx_type = line[50:55]
            rx = line[:50]
            rates = line[55:95]
        else:
            rx = line[:48]
            rates = line[53:93]

        rx = rx.replace('HV','hv').replace('L','l')
        react = rx[:20]
        prod = rx[20:].replace('hv','')
        rx_list = react.split() + prod.split()
        rx_name = ' + '.join(react.split())+" => "+' + '.join(prod.split())

        if 'M' in rx and rx_type == '2BODY':
            rx_type = 'three-body'

        if rx_type == "ELEMT":
            A, b, Ea = [float(a) for a in rates.split()]
            reaction['equation'] = rx_name
            reaction['rate-constant'] = {'A':A, 'b':b, 'Ea':Ea}
            if with_citations:
                reaction['citation'] = citation
        elif rx_type == "2BODY":
            
            react = react.replace('hv','')
            rx_name = ' + '.join(react.split())+" => "+' + '.join(prod.split())
            if len([float(a) for a in rates[:28].split()]) == 1:
                A = [float(a) for a in rates[:28].split()][0]
                Ea = 0
            else:
                A, Ea = [float(a) for a in rates[:28].split()]
            Ea = -Ea
            reaction['equation'] = rx_name
            reaction['rate-constant'] = {'A':A, 'b':0.0, 'Ea':Ea}
            if with_citations:
                reaction['citation'] = citation
                
        elif rx_type == "3BODY":
            a0,a1,cn,cm = [float(a) for a in rates.replace('D','E').split()]

            rx_name = ' + '.join(react.split())+" + M => "+' + '.join(prod.split())+' + M'

            A_low = a0*(300)**(cn)
            b_low = -cn
            A_high = a1*(300)**(cm)
            b_high = -cm
            reaction['equation'] = rx_name
            reaction['type'] = 'falloff'
            reaction['low-P-rate-constant'] = {'A':A_low, 'b':b_low, 'Ea':0.0}
            reaction['high-P-rate-constant'] = {'A':A_high, 'b':b_high, 'Ea':0.0}
            reaction['JPL'] = True
            if with_citations:
                reaction['citation'] = citation
        elif rx_type == "three-body":
            
            A, Ea = [float(a) for a in rates[:28].split()]
            Ea = -Ea
            reaction['equation'] = rx_name
            reaction['type'] = 'three-body'
            reaction['rate-constant'] = {'A':A, 'b':0.0, 'Ea':Ea}
            if with_citations:
                reaction['citation'] = citation

        elif rx_type == "2BD3P":
            
            A, Ea, b = [float(a) for a in rates[:34].split()]
            Ea = -Ea
            A = A*(1/298)**b
            
            reaction['equation'] = rx_name
            reaction['rate-constant'] = {'A':A, 'b':b, 'Ea':Ea}
            if with_citations:
                reaction['citation'] = citation
            
        elif rx_type == "3BD3P":
            
            A, Ea, b = [float(a) for a in rates[:34].split()]
            Ea = -Ea
            A = A*(1/298)**b
            rx_name = ' + '.join(react.split())+" + M => "+' + '.join(prod.split())+' + M'

            reaction['equation'] = rx_name
            reaction['type'] = 'three-body'
            reaction['rate-constant'] = {'A':A, 'b':b, 'Ea':Ea}
            if with_citations:
                reaction['citation'] = citation
            
        elif rx_type == "WEIRD":
            weird = True
            reaction['equation'] = rx_name
            if with_citations:
                reaction['citation'] = citation
            
        elif rx_type == "PHOTO" or rx_type[:-1] == "PHOT":
            reaction['equation'] = rx_name
            reaction['type'] = 'photolysis'
            if with_citations:
                reaction['citation'] = citation
        else:
            raise Exception(rx_type+' is an unknown reaction type')
                
        if weird:
            data['weird'].append(reaction)
        else:
            if "type" in reaction:
                if reaction['type'] == "photolysis":
                    data['photolysis'].append(reaction)
                else:
                    data['reactions'].append(reaction)
            else:
                data['reactions'].append(reaction)
    
    return data
    
def get_species(filename):

    fil = open(filename,'r')
    lines = fil.readlines()
    fil.close()

    species = []
    particles = []

    for line in lines[:-2]:
        if line[0]=='*':
            pass
        else:
            tmp = line.split()
            name = tmp[0].replace('CL','Cl')
            if 'AER' in name:
                particles.append(name)
            else:
                O = int(tmp[2])
                H = int(tmp[3])
                C = int(tmp[4])
                S = int(tmp[5])
                N = int(tmp[6])
                Cl = int(tmp[7])

                comp = {}

                atm = ['O','H','C','S','N','Cl']
                lst = [O,H,C,S,N,Cl]
                for i in range(len(lst)):
                    if lst[i] != 0:
                        comp[atm[i]] = lst[i]

                species.append({})
                species[-1]['name'] = name
                species[-1]['composition'] = comp
    return species, particles

def make_rx_yaml(species, particles, rx, photo_database = "Photochem"):
    out = {}
    out['reverse-reactions'] = False
    out['atoms'] = []
    out['atoms'].append({"name": "H", "mass": 1.00797, "redox": -0.5})
    out['atoms'].append({"name": "N", "mass": 14.0067, "redox": 0})
    out['atoms'].append({"name": "O", "mass": 15.9994, "redox": 1})
    out['atoms'].append({"name": "C", "mass": 12.011, "redox": -2})
    out['atoms'].append({"name": "S", "mass": 32.06, "redox": -2})
    out['atoms'].append({"name": "Cl", "mass": 35.453, "redox": -1})
    
    out['species'] = species
    
    possible_photos = generate_photo_yaml_entries([sp['name'] for sp in species])
    missing, not_missing = sort_photos(rx['photolysis'], possible_photos)
    
    if photo_database == "Photochem":
        out['reactions'] = rx['reactions'] + possible_photos
        out['missing'] = rx['weird']
    elif photo_database == "Atmos":
        out['reactions'] = rx['reactions'] + not_missing
        out['missing'] = rx['weird'] + missing
    else:
        raise Exception('"Photochem" and "Atmos" are the only options for photo_database')
        
    return out