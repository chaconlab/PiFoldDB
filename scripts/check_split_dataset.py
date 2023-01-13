import os, json, random
from collections import defaultdict

# Get SIFTS annotations from
# pdb_chain_cath_uniprot.tsv.gz

def write_fasta(splits, namekey):
    nfile = namekey + ".fa"
    with open(nfile, 'w') as outfile:
    	for name in splits[namekey] :
    		fname = 'data/seq/' + name + ".fa"
    		with open(fname) as infile:
    			outfile.write(infile.read())
    return 0

def write_tab(splits, cath_nodes, namekey):
    nfile = namekey + ".tab"
    with open(nfile, 'w') as outfile:
    	for name in splits[namekey] :
    		line = name
    		for t in cath_nodes[name] :
    			line = line + ' ' + t
    		outfile.write ('%s\n' % line)
    return 0       		
    		

def load_cath_topologies(cath_domain_file = 'data/cath-domain-list-4.3.0.txt'):
    print('Loading CATH domain nodes from', cath_domain_file)
    cath_nodes = defaultdict(list)
    with open(cath_domain_file,'r') as f:
        lines = [line.strip() for line in f if not line.startswith('#')]
        for line in lines:
            entries = line.split()
            cath_id, cath_node = entries[0], '.'.join(entries[1:4])
            chain_name = cath_id[:4].lower() + '.' + cath_id[4]
            cath_nodes[chain_name].append(cath_node)
    # Uniquify the list
    cath_nodes = {key:list(set(val)) for key,val in iter(cath_nodes.items())}
    return cath_nodes


if __name__ == '__main__':
    valid_fraction = 0.1

 
    # What topologies are in the Ollikainen set?
    print('Loading topologies from splits.json set')
    cath_nodes = load_cath_topologies()

    # Opening JSON file
    f = open('splits.json')
    
    # returns JSON object as a dictionary
    splits = json.load(f)
    dataset_names = splits.keys()

    print("Test:", len(splits['test']), "Train:", len(splits['train']), "Validation:",  len(splits['validation']))
    
    topo_test  = list(set([t for name in splits['test'] for t in cath_nodes[name]]))
    topo_train = list(set([t for name in splits['train'] for t in cath_nodes[name]]))
    topo_validation = list(set([t for name in splits['validation'] for t in cath_nodes[name]]))
    

    for name in splits['test'] :
        if not set(cath_nodes[name]).isdisjoint(topo_train) :
               print (name," present", set(cath_nodes[name]) )
    
    
    
    write_fasta(splits, 'validation')
    write_fasta(splits, 'train')
    write_fasta(splits, 'test')
    
    write_tab(splits, cath_nodes, 'test')
    write_tab(splits, cath_nodes, 'train')
    write_tab(splits, cath_nodes, 'validation')
    
    

    
    
    
    	
