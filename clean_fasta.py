import numpy as np
import matplotlib as mplt
import pandas as pd
import sys
import csv
from Bio import SeqIO
import os


## In every fasta file remove the gene sequences which has characters other than ['A','T','G','C']


def verify(seq):

	DNAchar = ['A','T','G','C']
	
	flag =1
	for i in seq:
		if i not in DNAchar:
			flag =0
			break
	return flag
	
		

main_fold = os.listdir('Genes_cytoplasm_raw')

for file_name in main_fold: 
	align = SeqIO.parse('Genes_cytoplasm_raw/'+str(file_name), 'fasta')

	for record in align:
    
		if verify(str(record.seq)):
			with open("Genes_cytoplasm/"+str(file_name)[:-4]+".fna", "a") as genefile:
				to_write = '>'+record.id+'\n'+str(record.seq)+'\n'
				print(to_write)
				genefile.write(str(to_write))

        

main_fold = os.listdir('Genes_membrane_raw')

for file_name in main_fold: 
	align = SeqIO.parse('Genes_membrane_raw/'+str(file_name), 'fasta')

	for record in align:
		if verify(str(record.seq)):
			with open("Genes_membrane/"+str(file_name)[:-4]+".fna", "a") as genefile:
				to_write = '>'+record.id+'\n'+str(record.seq)+'\n'
				print(to_write)
				genefile.write(str(to_write))
