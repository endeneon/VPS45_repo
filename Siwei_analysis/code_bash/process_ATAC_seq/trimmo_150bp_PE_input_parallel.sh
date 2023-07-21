#!/bin/bash

# Siwei 10 May 2023
# Trim 2x150 bp bulk ATAC-seq input for their adaptor contents

output_dir="trimmed_fastqs"

rm -r $output_dir
mkdir -p $output_dir

input_fastq_list="fastq_files_list.txt"

ls *_R1_001.fastq.gz > fastq_files_list.txt

trim_adaptor () {

	java="/data/zhangs3/software/jre1.8.0_291/bin/java"
	trimmomatic="/data/zhangs3/software/Trimmomatic-0.39/trimmomatic-0.39.jar"

	EACHFILE_1P=$1
	# echo $EACHFILE_1P
	EACHFILE_2P=${EACHFILE_1P/%_R1_001.fastq.gz/_R2_001.fastq.gz} 	

	output_dir=$2

	$java -jar \
		$trimmomatic \
		PE \
		-threads 6 \
		-phred33 \
		$EACHFILE_1P \
		$EACHFILE_2P \
		$output_dir/${EACHFILE_1P/%_R1_001.fastq.gz/_trimmed_R1_001.fastq.gz} \
		$output_dir/${EACHFILE_1P/%_R1_001.fastq.gz/_trimmed_U1_001.fastq.gz} \
		$output_dir/${EACHFILE_2P/%_R2_001.fastq.gz/_trimmed_R2_001.fastq.gz} \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
		LEADING:20 \
		TRAILING:25 \
		SLIDINGWINDOW:4:29 \
		MINLEN:36

}

export -f trim_adaptor

parallel \
	-a $input_fastq_list \
	-j 10 \
	trim_adaptor {} $output_dir

