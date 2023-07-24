#!/bin/bash

# Siwei 03 Jun 2022

# !!! input SAF file is the first parameter

# @@@ output file name is the second parameter

# will count all bam files under the current directory

mkdir -p featurecount_output

saf_reference=$1
output_file_name=$2

# make index for all BAM files in case
#for eachfile in *.bam
#do
#	echo $eachfile
#	samtools index -@ 20 $eachfile
#done

/home/zhangs3/Data/Tools/subread-2.0.3-Linux-x86_64/bin/featureCounts \
	-a $saf_reference \
	-B -C \
	-d 1 -D 2000 \
	-F SAF \
	-o featurecount_output/$output_file_name \
	-p \
	-T 32 \
	--countReadPairs \
	--ignoreDup \
	*.bam


