
from ._format import yaml, MyDumper, Loader, blockseqtrue, flowmap, FormatReactions_main

def photochem2cantera(infile, outfile):
    """Generates .yaml file compatible with Cantera from a Photochem reactions file.

    Parameters
    ----------
    infile : str
        Photochem .yaml reactions file
    outfile : str
        Output .yaml file in Cantera format.

    """
    fil = open(infile,'r')
    data = yaml.load(fil,Loader=Loader)
    fil.close()

    data = photochem2cantera_main(data)
            
    fil = open(outfile,'w')
    yaml.dump(data,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
    
def photochem2cantera_main(data):
    # remove particles
    if "particles" in data:
        del data['particles']

    # remove photolysis
    reactions_new = []
    for i in range(len(data['reactions'])):
        if "type" not in data['reactions'][i]:
            rxt = "elementary"
        else:
            rxt = data['reactions'][i]['type']
        if rxt != 'photolysis':
            reactions_new.append(data['reactions'][i])
            
        # check for JPL falloff. 
        if rxt == "falloff" and 'JPL' in data['reactions'][i]:
            raise Exception("Reaction "+data['reactions'][i]['equation']\
            +" has a JPL falloff function, which does not work with Cantera")
            
    data['reactions'] = reactions_new

    # get elements and delete
    elements = [a['name'] for a in data['atoms']]
    del data['atoms']

    # species
    species = [a['name'] for a in data['species']]

    # delete reverse-reactions
    if "reverse-reactions" in data:
        if data['reverse-reactions'] == False:
            raise Exception("Can only convert to Cantera if reactions are reversable")
        del data['reverse-reactions']  

    # Cantera can't handle > 2 temperature ranges for thermodynamic data
    for i in range(len(data['species'])):
        if len(data['species'][i]['thermo']['temperature-ranges']) > 3:
            data['species'][i]['thermo']['temperature-ranges'] = data['species'][i]['thermo']['temperature-ranges'][:3]
            data['species'][i]['thermo']['data'] = data['species'][i]['thermo']['data'][:2] 

    units = flowmap({"length": "cm", "quantity": "molec", "activation-energy": "K"})
    phases = [{}]
    phases[0]["name"] = "gas name"
    phases[0]["thermo"] = 'ideal-gas'
    phases[0]["elements"] = blockseqtrue(elements)
    phases[0]["species"] = blockseqtrue(species)
    phases[0]["kinetics"] = "gas"
    phases[0]["reactions"] = "all"
    phases[0]["state"] = flowmap({"T": 300.0,"P":1.0e5})

    data = FormatReactions_main(data)
    data['units'] = units
    data['phases'] = phases
    order = ['units','phases','species','reactions']
    copy = data.copy()
    data.clear()
    for key in order:
        if key in copy.keys():
            data[key] = copy[key]
    
    return data