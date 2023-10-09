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
    		
def load_chain_set(chain_file, verbose=True):
    chain_set = {}
    with open(chain_file,'r') as f:
        for i, line in enumerate(f):
            #print (line)
            #chain_name, jsons = line.split('\t')
            x = json.loads(line)
            #print(x)
            chain_set[x['name']] = x
            #print(x['name'])
            if verbose and (i + 1) % 1000 == 0:
                print('Loaded {} chains'.format(i+1))
    return chain_set

def load_cath_topologies(cath_domain_file = './cath-domain-list-4.3.0.txt'):
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

    cath_domain_file = 'data/cath-domain-list-4.3.0.txt'
    cath_json_file = 'chain_set800.jsonl'	
    cath_testset_file = 'testset_ingraham.txt'
    cath_json_out = 'splits.json'
        
    # We will make splits
    splits = {key:set() for key in ['test', 'train', 'validation']}
    
    # What topologies are in CATH?
    print('Loading topologies from CATH domains dict file:',cath_domain_file)
    cath_nodes = load_cath_topologies(cath_domain_file)
    

    # What topologies are in the test set?    
    print('Loading test set topologies from ', cath_testset_file, ' file')
    ingra = []
    with open(cath_testset_file,'r') as f:
       for line in f: 
        line = line.strip() #or some other preprocessing
        ingra.append(line) #storing everything in memory!

    topo_test = sorted(list(set(ingra)))
    test_topology_set = set(topo_test)
    #print(topo_test)


    # Build the test set from all topologies that appeared in the test set
    print('Loading CATH json file',cath_json_file )
    dataset = load_chain_set(cath_json_file)
    dataset_names = dataset.keys()
    #print(dataset_names)
    dataset_values = dataset.values()
    #print(dataset_values)

    
#    for name in dataset_names :
#        if ((not set(cath_nodes[name]).isdisjoint(test_topology_set)) and ( len(set(cath_nodes[name])) == 1)) :
#        	print (name," test", set(cath_nodes[name]), " len",  len(set(cath_nodes[name])) )
#        else :
#        	print (name," train", set(cath_nodes[name]), " len ", len(set(cath_nodes[name])))
         
                
  
    initial = [
        name for name in dataset_names
        if (not set(cath_nodes[name]).isdisjoint(test_topology_set) and ( len(set(cath_nodes[name])) <= 3)) 
    ]
                   
    topo_test2 = sorted(list(set([t for name in initial for t in cath_nodes[name]])))
    diferentes = sorted(set(topo_test2).difference(set(topo_test)))
            
    splits['test'] = [
        name for name in dataset_names
        if (not set(cath_nodes[name]).isdisjoint(test_topology_set) and (set(cath_nodes[name]).isdisjoint(diferentes)) and ( len(set(cath_nodes[name])) <= 3)) 
    ]
    
    topo_test2 = sorted(list(set([t for name in splits['test']  for t in cath_nodes[name]])))
    test_topology_set2 = set(topo_test2)
    print("Intial top", len(topo_test), " different", len(diferentes), " final ", len(topo_test2))
    
    
            
    splits['train'] = [
        name for name in dataset_names
        if (set(cath_nodes[name]).isdisjoint(test_topology_set2)) 
    ]

    valid_fraction = 0.05
    topo_train = list(set([t for name in splits['train'] for t in cath_nodes[name]]))
    num_topo_validation = int(valid_fraction * len(topo_train))
 
    random.seed(42)
    random.shuffle(topo_train)
    topo_validation = set(topo_train[:num_topo_validation])
   
    splits['validation'] = [
        name for name in splits['train']
        if not set(cath_nodes[name]).isdisjoint(topo_validation)
    ]

	
    	  
#    for name in splits['test'] :
#    	  print ("Test", name, set([t for t in cath_nodes[name]]))
    	  

    # Everything else is train
    splits['train'] = set(splits['train']).difference(set(splits['validation']))
    splits['train'] = set(splits['train']).difference(set(splits['test']))
    splits['train'] = list(splits['train'])
    
    ### just single chain to prevent overlap
    #splits['test'] = [
    #    name for name in dataset_names
    #    if ( (not set(cath_nodes[name]).isdisjoint(test_topology_set)) and ( len(set(cath_nodes[name])) == 1))
    #]



    print("Test:", len(splits['test']), "Train:", len(splits['train']), "Validation:",  len(splits['validation']))
    print("Total:", len(dataset_names))
    
    write_fasta(splits, 'validation')
    write_fasta(splits, 'train')
    write_fasta(splits, 'test')
    
    write_tab(splits, cath_nodes, 'test')
    write_tab(splits, cath_nodes, 'train')
    write_tab(splits, cath_nodes, 'validation')

    with open('splits.json', 'w') as f:
        json.dump(splits, f)
