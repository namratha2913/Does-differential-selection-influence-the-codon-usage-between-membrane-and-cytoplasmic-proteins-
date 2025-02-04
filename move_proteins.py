import shutil
import os
import csv


main_fold = os.listdir('Genes_all_raw')

def count_sequences(f):
	
	c = 0
	file = open("Genes_all_raw/"+str(f)) 
	for l in file:
		#print(l)
		if l.startswith(">"):
			c+=1
	return c

'''
gene_list = []
with open("cytoplasm_proteins_annotated.tsv") as file:

	tsv_file = csv.reader(file, delimiter="\t")
    
	for line in tsv_file:
		gene_list.append(line[6].split())
    
	print(gene_list)


for f in main_fold:
    
	#print("Genes_all_raw/"+str(f)[:-4])
        
	if [str(f)[:-4]] in gene_list:
		print("yes")
		c = count_sequences(f)
		print(c)
		
		if c>900:
        		shutil.copy("Genes_all/"+str(f), "Genes_cytoplasm_raw")

'''

gene_list = []
with open("membrane_proteins_annotated.tsv") as file:

	tsv_file = csv.reader(file, delimiter="\t")
    
	for line in tsv_file:
		gene_list.append(line[6].split())
    
	print(gene_list)


for f in main_fold:
    
	#print("Genes_all/"+str(f)[:-4])
        
	if [str(f)[:-4]] in gene_list:
		print("yes")
		c = count_sequences(f)
		print(c)
		
		if c>900:
        		shutil.copy("Genes_all_raw/"+str(f), "Genes_membrane_raw")
  	
       	
        
