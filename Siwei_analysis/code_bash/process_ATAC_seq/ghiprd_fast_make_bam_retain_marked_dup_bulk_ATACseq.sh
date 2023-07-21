#!/bin/bash

# Siwei 3 Apr 2018

# Siwei 13 Sept 2018
# Siwei 12 Dec 2019
# Siwei 04 May 2021
# Siwei 21 Apr 2023

bowtie2="/data/zhangs3/software/bowtie2-2.5.1-linux-x86_64/bowtie2"
jre_8="/data/zhangs3/software/jre1.8.0_291/bin/java"

source /data/zhangs3/software/cellranger-arc-2.0.2/sourceme.bash

mkdir -p pre_bams

for EACHFILE_1P in *_R1_001.fastq.gz
do
	####init
	date
	rm *.bam
	rm *.list
	rm *.table
	rm *.bai

	####variables
        echo $EACHFILE_1P
	echo ${EACHFILE_1P/%_R1_001.fastq.gz/_R2_001.fastq.gz}  #variable for 2P

	####alignments


	####align for paired reads
	date >> stat.txt
#	echo "Paired"
	echo ${EACHFILE_1P/%_R1_001.fastq.gz/} >> stat.txt
	($bowtie2 -p 72 -X 2000 \
		--mm --qc-filter --met 1 \
		-t \
		--sensitive \
		--no-mixed \
		--no-discordant \
		-x /data/zhangs3/Databases/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa \
		-1 $EACHFILE_1P \
		-2 ${EACHFILE_1P/%_R1_001.fastq.gz/_R2_001.fastq.gz} \
		| samtools sort \
		-l 9 -m 5G -@ 24 \
		-o merged_sorted.bam) \
		2>>stat.txt
	ls -lh
	
	cat stat.txt

	####merge reads and sort
#	samtools merge -l 0 -@ 23 -f merged_unsorted.bam Paired.bam Unpaired.bam
#	samtools sort -l 9 -m 2G -@ 23 -o merged_sorted.bam merged_unsorted.bam
	samtools index -@ 30 merged_sorted.bam
	ls -lh

	####dedup
	$jre_8 -Xmx80g -jar \
		/data/zhangs3/software/picard291.jar \
		MarkDuplicates \
		INPUT=merged_sorted.bam \
		OUTPUT=merged_sorted_dedup.bam \
		METRICS_FILE=metrics.txt \
		REMOVE_DUPLICATES=False \
		ASSUME_SORTED=True
	samtools index -@ 30 merged_sorted_dedup.bam

	####add read group
	$jre_8 -Xmx80g -jar /data/zhangs3/software/picard291.jar \
		AddOrReplaceReadGroups \
		INPUT=merged_sorted_dedup.bam \
		OUTPUT=pre_bams/${EACHFILE_1P/%_R1_001.fastq.gz/}.bam \
		RGID=${EACHFILE_1P/%_R1_001.fastq.gz/} \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=${EACHFILE_1P/%_R1_001.fastq.gz/}

	samtools index -@ 30 pre_bams/${EACHFILE_1P/%_R1_001.fastq.gz/}.bam

done
