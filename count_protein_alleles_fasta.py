from Bio import AlignIO
from Bio import SeqIO
import re
import pandas as pd
import os

#Count the number of protein sequence alleles

main_fold = os.listdir('protein_analysis/proteins_cytoplasm')

results =[]
for file_name in main_fold: 
	align = SeqIO.parse('protein_analysis/proteins_cytoplasm/'+str(file_name), 'fasta')
	to_df =[]
	allele_count = 0
	total_count = 0

	read_seq_list = []
	for record in align:
		total_count+=1
		if record.seq not in read_seq_list:
			read_seq_list.append(record.seq)
			allele_count+=1
			#print(record.seq)
		else:
			continue
            
        
	print(allele_count/total_count) #noemalizing based on total filtered sequences considered
	results.append([file_name.split('.')[0],allele_count/total_count,allele_count,total_count])

df = pd.DataFrame(results, columns=['Gene', 'Alleles_count','raw_allele_count','tot_genome_count'])
df.to_csv('Results_data/'+'count_protein_allele_cyto_fasta.csv', index=False)



results =[]
main_fold = os.listdir('protein_analysis/proteins_membrane/')

for file_name in main_fold: 
	align = SeqIO.parse('protein_analysis/proteins_membrane/'+str(file_name), 'fasta')

	to_df =[]
	allele_count = 0
	total_count = 0

	read_seq_list = []
	for record in align:
		total_count+=1
		if record.seq not in read_seq_list:
			read_seq_list.append(record.seq)
			allele_count+=1
			#print(record.seq)
		else:
			continue
            
        
	print(allele_count/total_count) #noemalizing based on total filtered sequences considered
	results.append([file_name.split('.')[0],allele_count/total_count,allele_count,total_count])

df = pd.DataFrame(results, columns=['Gene', 'Alleles_count','raw_allele_count','tot_genome_count'])
df.to_csv('Results_data/'+'count_protein_allele_membrane_fasta.csv', index=False)

