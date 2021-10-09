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

def FormatReactions(filename, outfilename):
    fil = open(filename,'r')
    data = yaml.load(fil,Loader=Loader)
    fil.close()
    
    order = ['reverse-reactions','atoms','species','particles','reactions']
    copy = data.copy()
    data.clear()
    for key in order:
        if key in copy.keys():
            data[key] = copy[key]

    # Species
    for i in range(len(data['species'])):
        
        if data['species'][i]['name'] == False:
            data['species'][i]['name'] = "NO"
        
        order = ['name', 'composition', 'thermo','note']
        copy = data['species'][i].copy()
        data['species'][i].clear()
        for key in order:
            if key in copy.keys():
                data['species'][i][key] = copy[key]
                    
        if 'thermo' in data['species'][i].keys():
            data['species'][i]['composition'] = flowmap(data['species'][i]['composition'])
            
            order = ['model','temperature-ranges','data','note']
            copy = data['species'][i]['thermo'].copy()
            data['species'][i]['thermo'].clear()
            for key in copy:
                data['species'][i]['thermo'][key] = copy[key]
                
            data['species'][i]['thermo']['temperature-ranges'] = blockseqtrue(data['species'][i]['thermo']['temperature-ranges'])
            
            data['species'][i]['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(data['species'][i]['thermo']['data'])]
        else:
            order = ['name', 'composition', 'thermo','note']
            copy = data['species'][i].copy()
            data['species'][i].clear()
            for key in order:
                data['species'][i][key] = copy[key]
            
            data['species'][i]['composition'] = flowmap(data['species'][i]['composition'])

    # Reactions
    for i in range(len(data['reactions'])):
        order = ['equation','type','rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
        copy = data['reactions'][i].copy()
        data['reactions'][i].clear()
        for key in order:
            if key in copy.keys():
                data['reactions'][i][key] = copy[key]
                
        flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
        for key in flowstyle:
            if key in data['reactions'][i].keys():
                data['reactions'][i][key] = flowmap(data['reactions'][i][key])
                
                
    fil = open(outfilename,'w')
    yaml.dump(data,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
    
if __name__ == "__main__":
    FormatReactions('data/reaction_mechanisms/zahnle_earth.yaml', 'zahnle_earth.yaml')
                
    