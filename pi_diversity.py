from Bio import AlignIO
from Bio import SeqIO
import re
import pandas as pd
import os
from itertools import combinations




def pairwise_differences(seq1, seq2):
    """Count nucleotide differences ignoring gaps and Ns"""
    k = 0
    L = 0
    for a, b in zip(seq1, seq2):
        if a in "ACGT" and b in "ACGT":
            L += 1
            if a != b:
                k += 1
    return k, L

    
main_fold = os.listdir('Genes_cyto_alleles')

results =[]
for file_name in main_fold: 
    align = SeqIO.parse('Genes_cyto_alleles/'+str(file_name), 'fasta')

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
            #print(record.seq)
        else:
            continue
            
    total = 0.0
    pair_count = 0

    for seq1, seq2 in combinations(read_seq_list, 2):
        k, L = pairwise_differences(seq1, seq2)
        if L > 0:
            total += k / L
            pair_count += 1

    pi = total / pair_count
           
    print(pi) #noemalizing based on total filtered sequences considered
    results.append([file_name.split('.')[0],pi])

df = pd.DataFrame(results, columns=['Gene', 'pi-nucleotide_diversity'])
df.to_csv('Results_data/'+'pi_cyto.csv', index=False)



results =[]
main_fold = os.listdir('Genes_membrane_alleles/')

for file_name in main_fold: 
    align = SeqIO.parse('Genes_membrane_alleles/'+str(file_name), 'fasta')

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
            #print(record.seq)
        else:
            continue
            
    total = 0.0
    pair_count = 0

    for seq1, seq2 in combinations(read_seq_list, 2):
        k, L = pairwise_differences(seq1, seq2)
        if L > 0:
            total += k / L
            pair_count += 1

    pi = total / pair_count
           
    print(pi) #noemalizing based on total filtered sequences considered
    results.append([file_name.split('.')[0],pi])

df = pd.DataFrame(results, columns=['Gene', 'pi-nucleotide_diversity'])
df.to_csv('Results_data/'+'pi_membrane.csv', index=False)


