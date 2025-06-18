from Bio import AlignIO
import re
import pandas as pd
import os
import statistics as st


codon_table = pd.read_csv('codon_usage.csv')
print(codon_table)

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
    
    checked_codon =[]  ###### codon parsed - already seen
    checked_aminoacid =[] ########## amino_acid aparsed- already seen
    variant_codon = [] ######### codon mutation
    variant_aa = [] ############ amino acid mutation
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
	degeneracy = int(codon_syn_numbers[codon_to_aminoacid(ref)])
	den_n = sum(list(codon_syn_numbers.values())) - degeneracy -1
	den_s = degeneracy -1
		
	diffn = len(aa_list)-1
	diffs = 1
	
	for i in uniqcodon_list:
		if i!=ref:
			if codon_to_aminoacid(i) == codon_to_aminoacid(ref):
				diffs+=1
	
	if codon_to_aminoacid(ref) ==  'None' or degeneracy==1:
		 degeneracy = 100		
	return den_n, diffn, den_s, diffs, degeneracy 
		
		
			
 
    
######################################################################################################################################################################

'''
main_fold = os.listdir('cyto_align_tsv/')
print(main_fold)
results = []

for file_name in main_fold:
    
    print(file_name)
    align = pd.read_csv('cyto_align_tsv/'+str(file_name),header=None)  
    shape = align.shape

    all_codon_pos =[]
    aval_s = 0
    diffs = 0
    
    ### 100 = non-synonymous mutation
    degn_count = {1:0,2:0,3:0,4:0,6:0,100:0}  ########## count total number of positions having codons with particular degeneracy (particular amino aicd coded by mutiple codons) ###################
    degn_pos_1 = {1:0, 100:0} ########## count total number of positions with degenaracy 'Deg-1': number of codons out of 'Deg-1' observed in the total protein length ############ 
    degn_pos_2 = {1:0, 2:0,100:0} ########## count total number of positions with degenaracy 'Deg-2': number of codons out of 'Deg-2' observed in the total protein length ############ 
    degn_pos_3 = {1:0, 2:0,3:0,100:0} ########## count total number of positions with degenaracy 'Deg-3': number of codons out of 'Deg-3' observed in the total protein length ############ 
    degn_pos_4 = {1:0, 2:0,3:0,4:0,100:0}########## count total number of positions with degenaracy 'Deg-4': number of codons out of 'Deg-4' observed in the total protein length ############ 
    degn_pos_6 = {1:0, 2:0,3:0,4:0,5:0,6:0,100:0} ########## count total number of positions with degenaracy 'Deg-6': number of codons out of 'Deg-6' observed in the total protein length ############ 
    dS_M = 1 
    
      
    for i in range(2,shape[1]):

        codon_list_ = list(align.iloc[:][i])
        check_syn = find_codon_at_mismatch(codon_list_)
        dnds_res = dnds_calc(check_syn[-1], check_syn[0], check_syn[1])
        #print(dnds_res)
        deg_ref = dnds_res[-1]
        degn_count[deg_ref]+=1
        if deg_ref == 1:
            if dnds_res[3] ==1:
                degn_pos_1[1]+=1
            if dnds_res[1] >=1:
                degn_pos_1[100]+=1

        elif deg_ref == 2:
            if dnds_res[3] ==1:
                degn_pos_2[1]+=1
            elif dnds_res[3] ==2:
                degn_pos_2[2]+=1
            if dnds_res[1] >=1:
                degn_pos_2[100]+=1
           
        elif deg_ref == 3:
            if dnds_res[3] ==1:
                degn_pos_3[1]+=1
            elif dnds_res[3] ==2:
                degn_pos_3[2]+=1
            elif dnds_res[3] ==3:
                degn_pos_3[3]+=1
            if dnds_res[1] >=1:
                degn_pos_3[100]+=1
            
        
        elif deg_ref == 4:
            if dnds_res[3] ==1:
                degn_pos_4[1]+=1
            elif dnds_res[3] ==2:
                degn_pos_4[2]+=1
            elif dnds_res[3] ==3:
                degn_pos_4[3]+=1
            elif dnds_res[3] ==4:
                degn_pos_4[4]+=1
            if dnds_res[1] >=1:
                degn_pos_4[100]+=1
           
        
        elif deg_ref == 6:
            if dnds_res[3] ==1:
                degn_pos_6[1]+=1
            elif dnds_res[3] ==2:
                degn_pos_6[2]+=1
            elif dnds_res[3] ==3:
                degn_pos_6[3]+=1
            elif dnds_res[3] ==4:
                degn_pos_6[4]+=1
            elif dnds_res[3] ==5:
                degn_pos_6[5]+=1
            elif dnds_res[3] ==6:
                degn_pos_6[6]+=1
            if dnds_res[1] >=1:
                degn_pos_6[100]+=1

        
    #print([degn_num[2]/degn_count[2], degn_num[3]/degn_count[3], degn_num[4]/degn_count[4], degn_num[6]/degn_count[6]])
    results.append([(file_name.split('.')[0]) ,shape[1]-1 , degn_pos_2[1]/degn_count[2] , degn_pos_2[2]/degn_count[2] ,degn_pos_2[100]/degn_count[2] ,degn_pos_3[1]/degn_count[3] ,degn_pos_3[2]/degn_count[3] ,degn_pos_3[3]/degn_count[3] ,degn_pos_3[100]/degn_count[3] ,degn_pos_4[1]/degn_count[4] ,degn_pos_4[2]/degn_count[4] ,degn_pos_4[3]/degn_count[4] ,degn_pos_4[4]/degn_count[4] ,degn_pos_4[100]/degn_count[4] ,degn_pos_6[1]/degn_count[6] ,degn_pos_6[2]/degn_count[6] ,degn_pos_6[3]/degn_count[6] ,degn_pos_6[4]/degn_count[6] ,degn_pos_6[5]/degn_count[6] , degn_pos_6[6]/degn_count[6],degn_pos_6[100]/degn_count[6]])
   

df = pd.DataFrame(results, columns=['Gene','genome_size_by#codon','Deg-21','Deg-22','Deg-20','Deg-31','Deg-32','Deg-33','Deg-30','Deg-41','Deg-42','Deg-43','Deg-44','Deg-40','Deg-61','Deg-62','Deg-63','Deg-64','Deg-65','Deg-66','Deg-60'])
df.to_csv('Results_data/'+'dSfracdist_cyto.csv', index=False)


############################################################################################################################################################
'''
    

main_fold = os.listdir('membrane_align_tsv/')

results = []

for file_name in main_fold:
    
    print(file_name)
    align = pd.read_csv('membrane_align_tsv/'+str(file_name),header=None)  
    shape = align.shape

    all_codon_pos =[]
    aval_s = 0
    diffs = 0
    
    ### 100 = non-synonymous mutation
    degn_count = {1:0,2:0,3:0,4:0,6:0,100:0}  ########## count total number of positions having codons with particular degeneracy (particular amino aicd coded by mutiple codons) ###################
    degn_pos_1 = {1:0, 100:0} ########## count total number of positions with degenaracy 'Deg-1': number of codons out of 'Deg-1' observed in the total protein length ############ 
    degn_pos_2 = {1:0, 2:0,100:0} ########## count total number of positions with degenaracy 'Deg-2': number of codons out of 'Deg-2' observed in the total protein length ############ 
    degn_pos_3 = {1:0, 2:0,3:0,100:0} ########## count total number of positions with degenaracy 'Deg-3': number of codons out of 'Deg-3' observed in the total protein length ############ 
    degn_pos_4 = {1:0, 2:0,3:0,4:0,100:0}########## count total number of positions with degenaracy 'Deg-4': number of codons out of 'Deg-4' observed in the total protein length ############ 
    degn_pos_6 = {1:0, 2:0,3:0,4:0,5:0,6:0,100:0} ########## count total number of positions with degenaracy 'Deg-6': number of codons out of 'Deg-6' observed in the total protein length ############ 
    dS_M = 1 
    
      
    for i in range(2,shape[1]):

        codon_list_ = list(align.iloc[:][i])
        check_syn = find_codon_at_mismatch(codon_list_)
        dnds_res = dnds_calc(check_syn[-1], check_syn[0], check_syn[1])
        #print(dnds_res)
        deg_ref = dnds_res[-1]
        degn_count[deg_ref]+=1
        if deg_ref == 1:
            if dnds_res[3] ==1:
                degn_pos_1[1]+=1
            if dnds_res[2] >=1:
                degn_pos_1[100]+=1

        elif deg_ref == 2:
            if dnds_res[3] ==1:
                degn_pos_2[1]+=1
            if dnds_res[3] ==2:
                degn_pos_2[2]+=1
            if dnds_res[2] >=1:
                degn_pos_2[100]+=1
           
        elif deg_ref == 3:
            if dnds_res[3] ==1:
                degn_pos_3[1]+=1
            if dnds_res[3] ==2:
                degn_pos_3[2]+=1
            if dnds_res[3] ==3:
                degn_pos_3[3]+=1
            if dnds_res[3] >=1:
                degn_pos_3[100]+=1
            
        
        elif deg_ref == 4:
            if dnds_res[3] ==1:
                degn_pos_4[1]+=1
            if dnds_res[3] ==2:
                degn_pos_4[2]+=1
            if dnds_res[3] ==3:
                degn_pos_4[3]+=1
            if dnds_res[3] ==4:
                degn_pos_4[4]+=1
            if dnds_res[3] >=1:
                degn_pos_4[100]+=1
           
        
        elif deg_ref == 6:
            if dnds_res[3] ==1:
                degn_pos_6[1]+=1
            if dnds_res[3] ==2:
                degn_pos_6[2]+=1
            if dnds_res[3] ==3:
                degn_pos_6[3]+=1
            if dnds_res[3] ==4:
                degn_pos_6[4]+=1
            if dnds_res[3] ==5:
                degn_pos_6[5]+=1
            if dnds_res[3] ==6:
                degn_pos_6[6]+=1
            if dnds_res[3] >=1:
                degn_pos_6[100]+=1
        
       
    results.append([(file_name.split('.')[0]) ,shape[1]-1 , degn_pos_2[1]/degn_count[2] , degn_pos_2[2]/degn_count[2] ,degn_pos_2[100]/degn_count[2] ,degn_pos_3[1]/degn_count[3] ,degn_pos_3[2]/degn_count[3] ,degn_pos_3[3]/degn_count[3] ,degn_pos_3[100]/degn_count[3] ,degn_pos_4[1]/degn_count[4] ,degn_pos_4[2]/degn_count[4] ,degn_pos_4[3]/degn_count[4] ,degn_pos_4[4]/degn_count[4] ,degn_pos_4[100]/degn_count[4] ,degn_pos_6[1]/degn_count[6] ,degn_pos_6[2]/degn_count[6] ,degn_pos_6[3]/degn_count[6] ,degn_pos_6[4]/degn_count[6] ,degn_pos_6[5]/degn_count[6] , degn_pos_6[6]/degn_count[6],degn_pos_6[100]/degn_count[6]])
   

df = pd.DataFrame(results, columns=['Gene','genome_size_by#codon','Deg-21','Deg-22','Deg-20','Deg-31','Deg-32','Deg-33','Deg-30','Deg-41','Deg-42','Deg-43','Deg-44','Deg-40','Deg-61','Deg-62','Deg-63','Deg-64','Deg-65','Deg-66','Deg-60'])
df.to_csv('Results_data/'+'dSfracdist_membrane.csv', index=False)
    
#'''
