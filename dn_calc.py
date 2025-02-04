from Bio import AlignIO
import re
import pandas as pd
import os
import statistics as st


codon_table = pd.read_csv('codon_usage.csv')
print(codon_table)
#print(codon_table.iloc[0][0])

# dictionary of Amino acid: #codons coding the amino acid

codon_syn_numbers = {'Glycine' : 4, 'Glutamate' : 2, 'Aspartate' : 2, 'Valine' : 4, 'Alanine': 4, 'Arginine' : 6, 'Lysine' : 2, 'Asparagine': 2 , 'Methionine': 1, 'Isoleucine': 3, 'Threonine' :4, 'Tryptophan':1, 'Cysteine':2, 'Stop':3, 'Tyrosine': 2, 'Phenylalanine': 2, 'Serine':6, 'Glutamine' : 2, 'Histadine': 2, 'Leucine': 6, 'Proline' : 4, 'None': 1}

#print(codon_syn_numbers['Glycine'])


# function to identify amino acid by given a as codon input 
def codon_to_aminoacid(codon):

    for i in range(codon_table.shape[0]):
        if codon == codon_table.iloc[i][1]:
            return     codon_table.iloc[i][0]
    return "None"
    
    
# function to identify the mismatches with respect to the reference sequence used
def find_codon_at_mismatch(codon_genome):
    
    checked_codon =[]
    checked_aminoacid =[]
    variant_codon = []
    variant_aa = []
    variant_syn = [] # append 1 if syn, append 2 if non-syn
    syn_check = 0
    nonsyn_check = 0
    both_check = 0
    all_codon =[]
    
    checked_aminoacid =[codon_to_aminoacid(codon_genome[0])]
    checked_codon = [codon_genome[0]]

    # number of alternative codons present relative to reference codon
    for i in codon_genome:
        #print(i)
        all_codon.append(i)
        if i not in checked_codon:
            checked_codon.append(i)
            aminoacid = codon_to_aminoacid(i)
            #print(aminoacid)
            if  aminoacid not in checked_aminoacid:
                checked_aminoacid.append(aminoacid)
                variant_codon.append(i)
                variant_aa.append(aminoacid)
                variant_syn.append(2)
                nonsyn_check =1
            elif aminoacid in checked_aminoacid:
                variant_codon.append(i)
                variant_aa.append(aminoacid)
                variant_syn.append(1)
                syn_check =1
        
        if syn_check and nonsyn_check:
            both_check =1
            
    return checked_codon, checked_aminoacid, variant_codon, variant_aa, variant_syn, syn_check, nonsyn_check, both_check, all_codon
   

def dnds_calc(codon_all,uniqcodon_list,aa_list):
	

	ref = st.mode(codon_all)
	den_n = sum(list(codon_syn_numbers.values())) - int(codon_syn_numbers[codon_to_aminoacid(ref)])-1
	den_s = int(codon_syn_numbers[codon_to_aminoacid(ref)])-1
		
	diffn = len(aa_list)-1
	diffs = 0
	
	for i in uniqcodon_list:
		if i!=ref:
			if codon_to_aminoacid(i) == codon_to_aminoacid(ref):
				diffs+=1
			
	return den_n, diffn, den_s, diffs
		
		
		
			
 
    
######################################################################################################################################################################

#'''
main_fold = os.listdir('cyto_align_tsv/')

results = []

for file_name in main_fold:
    
    print(file_name)
    align = pd.read_csv('cyto_align_tsv/'+str(file_name),header=None)  
    shape = align.shape

    total_loci_mut = 0 ## Total loci mismatches are found 
    syn_mismatch =0
    non_syn_mismatch =0
    both_mismatch =0
    all_codon_pos =[]
    aval_n = 0
    aval_s = 0
    diffn = 0
    diffs = 0
     
      
    for i in range(2,shape[1]):

        codon_list_ = list(align.iloc[:][i])

        check_syn = find_codon_at_mismatch(codon_list_)
        #all_codon_pos.append(check_syn[-1])
        #dnds_res = dnds_calc(check_syn[-1], check_syn[0], check_syn[1])
        #aval_n += dnds_res[0]
        #aval_s += dnds_res[2]
        diffn += len(check_syn[1])-1
        print(diffn)
        #diffs += dnds_res[3]
        
        '''if check_syn[5]:
            syn_mismatch+=1
        if check_syn[6]:
            non_syn_mismatch+=1
        if check_syn[7]:
            both_mismatch+=1'''
    
	   
    print((diffn/(shape[1]-2)))
    print(diffn)  

    results.append([(file_name.split('.')[0]),(diffn/(shape[1]-2))])
   

df = pd.DataFrame(results, columns=['Gene', 'Avg-non_syn_mut_per_codon'])
df.to_csv('Results_data/'+'dN_cyto.csv', index=False)


############################################################################################################################################################

    
#'''
main_fold = os.listdir('membrane_align_tsv/')

results = []

for file_name in main_fold:
    
    print(file_name)
    align = pd.read_csv('membrane_align_tsv/'+str(file_name),header=None)  
    shape = align.shape

    total_loci_mut = 0 ## Total loci mismatches are found 
    syn_mismatch =0
    non_syn_mismatch =0
    both_mismatch =0
    aval_n = 0
    aval_s = 0
    diffn = 0
    diffs = 0
     
      
    for i in range(2,shape[1]):

        codon_list_ = list(align.iloc[:][i])
        
        check_syn = find_codon_at_mismatch(codon_list_)
        #all_codon_pos.append(check_syn[-1])
        #dnds_res = dnds_calc(check_syn[-1], check_syn[0], check_syn[1])
        #aval_n += dnds_res[0]
        #aval_s += dnds_res[2]
        diffn += len(check_syn[1])-1
        print(diffn)
        #diffs += dnds_res[3]
        
        '''if check_syn[5]:
            syn_mismatch+=1
        if check_syn[6]:
            non_syn_mismatch+=1
        if check_syn[7]:
            both_mismatch+=1'''
        
    print((diffn/(shape[1]-2)))  

    results.append([(file_name.split('.')[0]),(diffn/(shape[1]-2))])
   

df = pd.DataFrame(results, columns=['Gene', 'Avg-non_syn_mut_per_codon'])
df.to_csv('Results_data/'+'dN_membrane.csv', index=False)
    
#'''
