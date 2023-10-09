import os, time, gzip, urllib, json
import mmtf
from collections import defaultdict
from urllib.error import URLError, HTTPError
import gemmi 
from gemmi import cif
import numpy as np


subs_one= {
'2AS':'D', 'CMH':'C', 'CS3':'C', 'CSA':'C', 'CSK':'C', 'CSR':'C', 'CSU':'C', 'CSZ':'C', 'CCS':'C', 'CSB':'C', 'CSP':'C', 'MCS':'C',
'DMG':'G', 'DAH':'A', 'DDZ':'A', 'KYN':'A', '4GJ':'A', 'LYX':'K', 'M2L':'M', 'MSO':'M', 'OMT':'M', 'PXU':'P', 'ELY':'K', 'FHL':'K',
'HSO':'H', 'TRO':'W', 'TRN':'W', 'TYI':'Y', 'BYR':'Y', 'LP6':'L', 'ACE':'S', 'SEB':'S', 'ALS':'S'
}

def remove_char(s, i):
    b = bytearray(s, 'utf-8')
    del b[i]
    return b.decode()
    
def write_fasta(seq, name, wrap=80):
    fasta_file = "data/seq/" +name + ".fa" 
    with open(fasta_file, 'w') as f:
    	f.write('>%s\n' % (name))
    	for i in range(0, len(seq), wrap):
    		f.write('%s\n' % (seq[i:i + wrap]))

def download_cached(url, target_location):
    """ Download with caching """
    target_dir = os.path.dirname(target_location)
    if not os.path.isfile(target_location):
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # Use MMTF for speed

        try:
           response = urllib.request.urlopen(url)
           size = int(float(response.headers['Content-Length']) / 1e3)
           print('Downloading {}, {} KB'.format(target_location, size))
           with open(target_location, 'wb') as f:
              f.write(response.read())
        except HTTPError as e:
            # do something
            print('Error code: ', e.code)
           
        
    return target_location


def mcif_fetch(pdb, cache_dir='data/cifs/'):

    """ Retrieve cif record from PDBREDO or PDB  with local caching """
    cif_file = cache_dir + pdb + '.cif'
    cif_filegz = cache_dir + pdb + '.cif.gz'
    cif_file_pdbredo = cache_dir + pdb + '_final.cif'
    if not os.path.exists(cache_dir):
    	os.makedirs(cache_dir)
    	
    if not os.path.isfile(cif_filegz):
    	url = "https://pdb-redo.eu/db/"+ pdb+ "/"+ pdb +"_final.cif"	
    	print (url)
    	fileA = download_cached(url, cif_file_pdbredo)

    	if not os.path.isfile(cif_filegz):
    		url = "https://files.rcsb.org/download/"+  pdb +".cif.gz"		
    		print (url)
    		fileA = download_cached(url, cif_filegz)
    	else :
    		cmd = "mv " + cif_file_pdbredo + " " + cif_file +"; gzip " +  cif_file
    		print (cmd)
    		os.system(cmd)

	
    if not os.path.isfile(cif_filegz):
    	print("ERROR downloading",pdb)
    else :
    	st = gemmi.read_structure(cif_filegz)
    	st.remove_alternative_conformations()
    	st.remove_hydrogens()
    	st.remove_ligands_and_waters()
    	st.remove_empty_chains()


    return st


def mcif_parse(pdb_id, sel_chain, target_atoms = ['N', 'CA', 'C', 'O']):

	st = mcif_fetch(pdb_id)
# Build a dictionary
	mmtf_dict = {}
	mmtf_dict['seq'] = []
	mmtf_dict['coords'] = {code:[] for code in target_atoms}



	firstRes=1
	oldResId=1
	resId = 0
	oldResId = 0
	gap = 0



	detectC = 0
	chainsD =" "
	for ent in st.entities:
		#if ((ent.polymer_type == gemmi.PolymerType.PeptideL) or (ent.polymer_type == gemmi.PolymerType.Unknown))  :  # not working PDB is wrongly annotated
		if  (detectC == 0 ) :
			for subchain in ent.subchains :
				#print("Subchain", subchain," ", ent.polymer_type )
				chainsD = chainsD + " " + subchain
				if subchain == sel_chain:
					selA = ent.name
					seqF = gemmi.one_letter_code(ent.full_sequence)
					mmtf_dict['num_chains'] = len(ent.subchains)
					detectC = 1
					break
		if (detectC == 1 ) :
			break
	
	if (detectC == 0 ) :
		mmtf_dict['num_chains'] = len(ent.subchains) 
		#print("warning chain not detected in", pdb_id," ",sel_chain," present", chainsD)
	
	# print(seqF,"len", len(seqF))
					
	

        ### first round to get the sequence with gaps
	seq0 = ''
	nres = 0
	for chain in st['1']:
	      if chain.name == sel_chain:
		      for residue in chain:
		      		detectC = 1		      		
		      		resId0 = int(residue.seqid.num)
		      		if (firstRes==1):
		      			oldResId = int(residue.seqid.num)
		      			firstRes = 0

		      		gap = np.int_(resId0*1 - oldResId*1)
		      		if (gap > 1) :
		      			print ("Gap ", gap, " ", pdb_id, " ", oldResId+1,"-", oldResId+gap)
		      			for x in range(gap):
		      				seq0=seq0+"X";
		      				nres+=1;		
		      		else :
		      			aa = gemmi.find_tabulated_residue(residue.name).one_letter_code
		      			if (aa == " ") :
		      				if residue.name in subs_one:
		      					aa = subs_one[residue.name]
		      					#print (" subs", aa, " ", residue.name) 
		      				else :
		      					aa = "X"
		      										
		      			seq0=seq0+aa
		      			nres+=1;
		      						
		      		oldResId = resId0
		      break

	if (detectC == 0 ) :
		print("warning chain not detected in", pdb_id," ",sel_chain," present", chainsD)
		return mmtf_dict


	chain_len = len(seq0)
	#print(seq0,"len", len(seq0), nres)

	seq = ''
	nres = 0
	coords_null0 = [float('nan')] * 3
	detectC = 0
	for chain in st['1']:
	      if chain.name == sel_chain:
		     	      
		      coords_null = [[float('nan')] * 3] * (chain_len)
		      mmtf_dict['coords'] = {code : list(coords_null) for code in target_atoms}
		      #print ("mierda ", chain.name, " ", chain[1].seqid.num,"-", chain[-1].seqid.num)
		      for residue in chain:
		      		
		      		resId0 = int(residue.seqid.num)
		      		if (firstRes==1):
		      			oldResId = int(residue.seqid.num)
		      			firstRes = 0

		      		gap = np.int_(resId0*1 - oldResId*1)
		      		if (gap > 1) :
		      			print ("Gap ", gap, " ", pdb_id, " ", oldResId+1,"-", oldResId+gap)
		      			for x in range(gap):
		      				# print ("X  ", oldResId+x+1)
		      				seq=seq+"X";
		      				for code in target_atoms:
		      					mmtf_dict['coords'][code][nres] = coords_null0
		      				nres+=1;		
		      		else :
		      			has_check=0;
		      			for atom in residue:
		      				for code in target_atoms:
		      					if atom.name == code:
		      						#print(atom.pos[0]," ",atom.pos[1]," ",atom.pos[2]) 
		      						xyz =  [atom.pos[0],atom.pos[1],atom.pos[2]]
		      						mmtf_dict['coords'][code][nres] = xyz 
		      						has_check += 1
		      			aa = gemmi.find_tabulated_residue(residue.name).one_letter_code
		      			if (aa == " ") :
		      				if residue.name in subs_one:
		      					aa = subs_one[residue.name]
		      					#print (" subs", aa, " ", residue.name) 
		      				else :
		      					has_check += 100
		      										
		      			if ( has_check == 4 ) :
		      				#print (gemmi.find_tabulated_residue(residue.name).one_letter_code, " ",residue.seqid.num, "type-", residue.entity_type,"-")

		      				
		      				seq=seq+aa

		      				  
		      				nres+=1;
		      				
		      			else :
		      				print ("Missing atoms ", pdb_id, " ", has_check, " res_name", residue.name, "(",gemmi.find_tabulated_residue(residue.name).one_letter_code,") res_num ",residue.seqid.num)
		      				seq=seq+"X";
		      				for code in target_atoms:
		      					mmtf_dict['coords'][code][nres] = coords_null0
		      				nres+=1;

		      			        		
		      			
		      		oldResId = resId0
		      break
		      

   
	# print("se0 ",seq)
	###### remove gaps from Nt and Ct ######
	for runs in range(6):	
		jj=5
		contX =0;
		while jj < 10:
			if (seq[jj] == 'X') :
				contX +=1;
				# print("seq ",seq[jj]," ",jj)
			jj += 1
		
		if (contX == 5) :
			for index in range(5):
				seq = seq[:index] + "X" + seq[index + 1:]
		

		jj = nres -10
		contX =0;

		while jj < nres - 5 :
			#print("seq-> ",seq[jj]," ",jj)
			if (seq[jj] == 'X') :
				contX +=1;
			jj += 1
		
		if (contX == 5) :
			seq = seq[:(nres - 5)] 
			for index in range(5):
				seq = seq + "X" 
			
		

		# print("se0 ",seq)
		     
	    ### Remove at Nt missing res
		jj=0
		while jj < 1:
			detect =0  
			if (seq[jj] == 'X') :
				for code in target_atoms :
					del mmtf_dict['coords'][code][jj]
				seq = seq[1:] 
				nres -= 1
			else:
				jj=1
				

		
	   ### Remove at Ct missing res
		jj=nres -1
		while jj == nres -1:
			detect =0
			if (seq[jj] == 'X') :
				for code in target_atoms :
					del mmtf_dict['coords'][code][jj]
				seq = seq[:-1]
				nres -= 1
				jj -=1
			else:
				jj=nres+1

	 # print("seqF ",seq)   	
	

			
	
	if (len(seq)<=2) :
		print("chain not detected in", pdb_id," ",sel_chain," present", chainsD)
	
	contX=0
	MaxGAP=0
	countGAP =0
	for x in range(len(seq)) :
		if (seq[x] == 'X') :
			contX +=1
			countGAP +=1
		else :
			if (countGAP>MaxGAP) :
				MaxGAP = countGAP
			countGAP=0
				
	if (contX/len(seq) > 0.3)  :
		print ("removed sequence with gaps >1/3 of the length", pdb_id," ",sel_chain," ", seq)
		mmtf_dict['seq'] = []
		mmtf_dict['coords'] = []
		return mmtf_dict
		
	if (MaxGAP > 100 )  :
		print ("removed GAP > 100 aa", pdb_id," ",sel_chain,"MaxGAP",MaxGAP," ", seq)
		mmtf_dict['seq'] = []
		mmtf_dict['coords'] = []
		return mmtf_dict
		
	if (len(seq)-contX< 20)  :
		print ("removed sequence less than 20 aa", pdb_id," ",sel_chain," ", seq)
		mmtf_dict['seq'] = []
		mmtf_dict['coords'] = []
	else :

		x=0
		limit =len(seq)-1
		while x < limit:
			if (seq[x] == 'X') :
			   for code in target_atoms:
			     mmtf_dict['coords'][code][x] = coords_null0
			   jj=x+1	
			   while jj < len(seq)-1:
			     if (seq[jj] == 'X') :
			        seq = remove_char(seq,jj)
			        for code in target_atoms :
			        	del mmtf_dict['coords'][code][jj]
			        limit-=1
			     else: 
			        jj=1000000000	
			x+=1	     

		fasta = pdb_id+"."+sel_chain
		fasta_file = "data/seq/" + fasta + ".fa"
		if not os.path.isfile(fasta):
			write_fasta(seq.upper(), fasta, 80)
		mmtf_dict['seq'] = seq.upper()
		if (len(mmtf_dict['coords']['CA']) != len(seq)):
			print ("\n\n Error", len(mmtf_dict['coords']['CA']), " seq", len(seq))

	
	return mmtf_dict
