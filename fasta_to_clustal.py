from Bio import AlignIO
import os

'''
main_fold = os.listdir('alignment_cyt/')


for file_name in main_fold:

    AlignIO.convert("alignment_cyt/"+str(file_name), "fasta", "cyto_clustal_align/"+str(file_name)[:-5]+"clu", "clustal")
    
'''    
main_fold = os.listdir('alignment_membrane/')


for file_name in main_fold:

    AlignIO.convert("alignment_membrane/"+str(file_name), "fasta", "membrane_clustal_align/"+str(file_name)[:-5]+"clu", "clustal")
    

