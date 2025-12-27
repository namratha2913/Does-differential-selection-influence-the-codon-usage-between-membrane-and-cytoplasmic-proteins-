# This project includes the code used in the research study titled: Genomic signatures of asymmetric selection on membrane and cytoplasmic proteins.
The following code are used to analyze genome data and identify genomic signatures of the evolution of cytoplasmic and membrane proteins.

1. multiple_files_macse.sh - To perform multiple sequence alignment
2. count_alleles_fasta.py - To count the number of Coding DNA sequence alleles
3. count_protein_alleles_fasta.py - To count the number of protein sequence alleles. We use data generated from files 2,3 to calculate the fraction of synonymous alleles for the genes.
4. dn_calc.py - To quantify average non-synonymous substitutions per codon
5. ds_fraction_dist.py - To quantify the number of alternate synonymous codons used at all the amino acid sites with degeneracy> 1.

List of _E. coli_ genes used in the study:

- Cytoplasmic_prot_Ecoli.csv
- Membrane_prot_Ecoli.csv

List of _S. cerevisiae_ genes used in the study:

- Cyto-prot_Scer.csv
- Membrane-prot_Scer.csv
