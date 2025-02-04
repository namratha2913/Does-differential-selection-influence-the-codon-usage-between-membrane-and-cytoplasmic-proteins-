for file in /home/namratha/Namratha/CUB_codon_usage_bias/1-3-24Analysis/Results_all/alignment_cyt/*; do

	echo $file
	echo $(basename -a $file)
	#echo $pseudo_file
	pseqsid -i $file

done
