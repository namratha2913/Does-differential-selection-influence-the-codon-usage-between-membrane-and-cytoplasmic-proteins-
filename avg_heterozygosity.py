import os
from Bio import SeqIO
import itertools
from collections import Counter
import pandas as pd



main_fold = os.listdir('Genes_cytoplasm')
results =[]

for file_name in main_fold: 
    align = SeqIO.parse('Genes_cytoplasm/'+str(file_name), 'fasta')
    sequences = []

    for record in align:
        sequences.append(str(record.seq).upper())
    print("1")   
    n = len(sequences)
    if n < 2:
        continue
    print("11")
    allele_counts = Counter(sequences)
    H = 1.0 - sum((count / n) ** 2 for count in allele_counts.values())
    print(n)
    print(len(allele_counts))
    print(H) #noemalizing based on total filtered sequences considered
    results.append([file_name.split('.')[0],H])

df = pd.DataFrame(results, columns=['Gene', 'Heterozygosity'])
df.to_csv('Results_data/'+'Hetero_cyto.csv', index=False)

 
##################################################################################################
   
main_fold = os.listdir('Genes_membrane')
results =[]

for file_name in main_fold: 
    align = SeqIO.parse('Genes_membrane/'+str(file_name), 'fasta')
    sequences = []

    for record in align:
        sequences.append(str(record.seq).upper())
        
    n = len(sequences)
    if n < 2:
        continue
    
    allele_counts = Counter(sequences)
    H = 1.0 - sum((count / n) ** 2 for count in allele_counts.values())
    
    print(H) #noemalizing based on total filtered sequences considered
    results.append([file_name.split('.')[0],H])

df = pd.DataFrame(results, columns=['Gene', 'Heterozygosity'])
df.to_csv('Results_data/'+'Hetero_membrane.csv', index=False)
