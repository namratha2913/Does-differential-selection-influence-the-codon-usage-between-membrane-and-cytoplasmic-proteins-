import numpy as np
import matplotlib as mplt
import pandas as pd
import sys
import csv
from Bio import SeqIO
import os



df = pd.read_csv('codon_usage.csv')

aa = list((df['AminoAcid']))
codon = list(df['Codon'])
codon_to_aa = dict(zip(codon,aa))


df = pd.read_csv('amino-acids.csv')

codon = list((df['full name']))
symbol = list(df['single letter code'])

aa_symbol = dict(zip(codon,symbol))

print(codon_to_aa)
print(aa_symbol)

def translate(seq):
    
    protein =""
    if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                #print(codon)
                s = codon_to_aa[codon]
                protein+= aa_symbol[s]

    return protein


main_fold = os.listdir('Genes_cytoplasm')

for file_name in main_fold: 
    align = SeqIO.parse('Genes_cytoplasm/'+str(file_name), 'fasta')

    for record in align:
        with open("protein_analysis/proteins_cytoplasm/"+str(file_name)[:-4]+".fna", "a") as genefile:
            to_write = '>'+record.id+'\n'+translate(record.seq)+'\n'
            print(to_write)
            genefile.write(to_write)
            
'''          

main_fold = os.listdir('Genes_membrane')

for file_name in main_fold: 
    align = SeqIO.parse('Genes_membrane/'+str(file_name), 'fasta')

    for record in align:
        with open("protein_analysis/proteins_membrane/"+str(file_name)[:-4]+".fna", "a") as genefile:
            print(file_name)
            to_write = '>'+record.id+'\n'+translate(record.seq)+'\n'
            print(to_write)
            genefile.write(to_write)
            
'''
