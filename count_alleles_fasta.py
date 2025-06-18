from Bio import AlignIO
from Bio import SeqIO
import re
import pandas as pd
import os


main_fold = os.listdir('SIMGenes_alleles/cyto_seqs')

results =[]
for file_name in main_fold: 
	align = SeqIO.parse('/home/namratha/Namratha/NewDrive_link/CUB_codon_usage_bias/Simulate-null-hypothesis/Ecoli/SIMGenes_alleles/cyto_seqs/'+str(file_name), 'fasta')
	print(file_name)
	#print ("Alignment length %i" % align.get_alignment_length())

	to_df =[]
	#ref = align[0].seq
	allele_count = 0
	total_count = 0

	read_seq_list = []
	for record in align:
		total_count+=1
		if record.seq not in read_seq_list:
			read_seq_list.append(record.seq)
			allele_count+=1
			print(record.seq)
		else:
			continue
            
        
	print(allele_count/total_count) #normalizing based on total filtered sequences considered
	results.append([file_name.split('.')[0],allele_count])

df = pd.DataFrame(results, columns=['Gene', 'Alleles_count'])
df.to_csv('Results_data/'+'count_allele_cyto_fasta.csv', index=False)



results =[]

main_fold = os.listdir('/home/namratha/Namratha/NewDrive_link/CUB_codon_usage_bias/Simulate-null-hypothesis/Ecoli/SIMGenes_alleles/membrane_seqs')

for file_name in main_fold: 

	align = SeqIO.parse('/home/namratha/Namratha/NewDrive_link/CUB_codon_usage_bias/Simulate-null-hypothesis/Ecoli/SIMGenes_alleles/membrane_seqs/'+str(file_name), 'fasta')

	print(file_name)
	#print ("Alignment length %i" % align.get_alignment_length())


	to_df =[]
	#ref = align[0].seq
	allele_count = 0
	total_count = 0

	read_seq_list = []
	for record in align:
		total_count+=1
		if record.seq not in read_seq_list:
			read_seq_list.append(record.seq)
			allele_count+=1
			print(record.seq)
		else:
			continue
            
        
	print(allele_count/total_count) #noemalizing based on total filtered sequences considered
	results.append([file_name.split('.')[0],allele_count])

df = pd.DataFrame(results, columns=['Gene', 'Alleles_count'])
df.to_csv('Results_data/'+'count_allele_membrane_fasta.csv', index=False)

