import yaml

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
            
            order = ['name', 'composition', 'condensate', 'thermo','note']
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

    # Particles
    if 'particles' in data:
        for i in range(len(data['particles'])):
            data['particles'][i]['composition'] = flowmap(data['particles'][i]['composition'])
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
                order = ['type','vdep','mix','flux','height']
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

def mechanism_dict_with_atoms(dat, atoms_names, remove_particles=False, remove_reaction_particles=False):

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
        if not exclude:
            species.append(sp)

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
            
            if not exclude:
                particles.append(sp)

    if "reactions" in dat:
        reactions = []
        for i,rx in enumerate(dat['reactions']):
            eq = rx['equation']
            eq = eq.replace('(','').replace(')','')
            if '<=>' in eq:
                split_str = '<=>'
            else:
                split_str = '=>'
    
            a,b = eq.split(split_str)
            a = a.split('+')
            b = b.split('+')
            a = [a1.strip() for a1 in a]
            b = [b1.strip() for b1 in b]
            sp = a + b
    
            exclude = False
            for s in sp:
                for tmp in comp[s]:
                    if tmp in exclude_atoms:
                        exclude = True
                        break
                if exclude:
                    break
            if not exclude:
                reactions.append(rx)
                
    out = dat
    out['atoms'] = atoms
    out['species'] = species
    if 'particles' in dat and not remove_particles:
        out['particles'] = particles
    if 'reactions' in dat:
        out['reactions'] = reactions

    return out

def resave_mechanism_with_atoms(infile, outfile, atoms_names, remove_particles=False, remove_reaction_particles=False):
    
    with open(infile,'r') as f:
        dat = yaml.load(f,Loader=Loader)
    
    out = mechanism_dict_with_atoms(
        dat, 
        atoms_names, 
        remove_particles,
        remove_reaction_particles
    )

    out = FormatReactions_main(out)
    with open(outfile,'w') as f:
        yaml.dump(out,f,Dumper=MyDumper,sort_keys=False,width=70)

        