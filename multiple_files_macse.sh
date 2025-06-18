########## Shell script to perform multiple sequence analysis #############

for file in Genes_cyto_alleles/*; do
	output=$"alignment_cyt/"$(basename -a $file)-"output.fasta"
	output_AA=$"protein_analysis/protein_alignment_cyt/"$(basename -a $file)-"AAalign.fasta"
	#pseudo_file=$"pseudo_cyto/"$(basename -a $file)-"pseudo.fasta"

echo $file
echo $(basename -a $file)
#echo $pseudo_file
java -jar macse_v2.07.jar -prog alignSequences -seq $file -out_NT $output -out_AA $output_AA -max_refine_iter 0

done
