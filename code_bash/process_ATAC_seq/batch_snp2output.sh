#!/bin/bash

# Siwei 27 Apr 2023

# batch creating WASP calib directories from multiple vcf files
# integrate code from snp2output.sh

# 1hr_CD_03_NEFM_neg_glut_hg38only_scATAC_no_VQSR_uncalib.vcf
# Ast-31_S9_trimmed.unfiltered.vcf


vcf_list=( $(ls *.vcf ) )

cell_line_string=( $(ls *.vcf \
        | sed 's/_trimmed/\t/g' \
        | cut -f 1 ) )

echo ${cell_line_string[@]}

for ((i=0; i<${#vcf_list[@]}; i++ ))
do

	INPUT_VCF=${vcf_list[$i]}
	OUTPUT_DIR="WASP_calib"/${cell_line_string[$i]}

	k=1
	char_chr='chr'
	rm -r $OUTPUT_DIR
	
	mkdir -p $OUTPUT_DIR
	
	while [ $k -lt 23 ]
	do
	        echo $char_chr$k
	        k_chr=$char_chr$k
	        echo '$char_chr='$char_chr
	        echo '$k_chr='$k_chr
	        echo '$k='$k
	        OUTPUT_FILE=$OUTPUT_DIR/$k_chr.snps.txt.gz
	        cat $INPUT_VCF | awk -F '\t' -v awk_k=$k_chr '$1 == awk_k {print $2"\t"$4"\t"$5}' | gzip > $OUTPUT_FILE
	        k=$[$k+1]
	done

done


