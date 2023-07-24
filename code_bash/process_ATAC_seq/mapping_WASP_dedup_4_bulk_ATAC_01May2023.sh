#!/bin/bash

# Siwei 01 May 2023
# Adapt for bulk ATAC-seq input

# Siwei 27 Apr 2023
# rewrite for Duan_025 series

# Siwei 27 Sept 2022
# Solve the non-unique BAM read name issue
# It appears that CellRanger only make unique read names (L+R use the same name, column 1) within the same barcode
# Hence WASP might discarded large amount of reads due to duplicated read names (across different barcodes)
# Need to split and WASP-calibrate BAM files at the barcode level and merge all

# Siwei 26 Aug 2022
# UnmarkDuplicates and AddOrRemoveReadGroups have to be re-written using samtools
# instead of GATK+picard due to the aforementioned JAVA VM int max issue.
# Apparently splitting bams by chr did not solve the problem
# use samtools view --remove-flags and 
#     samtools addreplacerg

# Need to use samtools v1.14 (anaconda env aligners) by define $samtools in variables
# ~/Data/Anaconda3-envs/aligners/bin/samtools


# Siwei 25 Aug 2022
# Pre-split large BAM file by chromosome and process sequentially
# Since Java VM int max value is 2^31-1 (32-bit, regardless of 32/64 bit system), the one used in Picard and GATK,
# which corresponding to approximately ~2G BAM file.
# Integer larger than this should use long format

# Siwei 16 Aug 2022
# Randomly remove duplication use WASP rmdup at cell-barcode level
# by de-convoluting reads by CR:Z with 10x subset-bam;
# Use barcode list from the original barcode lists generated from RNA-Seq + genotyping that contains
# cell type + time + line information.

# Run deconvolution and duplication removal steps in parallel
# using GNU parallel, 18 jobs in total;
# Remember to pass variables as arguments into functions (or they cannot be found).

# Use system java openjdk version "1.8.0_332", probably faster than 
# the previously used jre1.8.0_291 for Picard AddOrRemoveReadGroups on Picard 2.25.4

# update Picard to v2.25.4 (the same one comes with GATK 4.1.8.1)

###############

# Siwei 16 Jul 2022

# Siwei 10 May 2021
# customise code for genres01

# Siwei 20 Dec 2019
# updated WASP to 0.3.4 (30 Apr 2019)


# Set these environment vars to point to
# your local installation of WASP

# source /home/zhangs3/Data/Tools/cellranger-arc-2.0.2/sourceme.bash

export OMP_NUM_THREADS=8

WASP="/home/zhangs3/Data/Tools/WASP"
# jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin/java"
picard_291="/home/zhangs3/Data/Tools/picard291.jar"
picard_2254="/home/zhangs3/Data/Tools/gatk-4.1.8.1/picard2_25_4.jar"
gatk4="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"
samtools="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/samtools"


DATA_DIR="/home/zhangs3/NVME/package_temp/WASP_temp_1"


OUTPUT_FINAL="bowtie2_WASPed_BAMs"

# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2
INDEX="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
# use BWA reference from Cell Ranger_ref
# bwa_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

# set error trap
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

mkdir -p $OUTPUT_FINAL

# 1hr_CD_03_NEFM_neg_glut_hg38only.bam

# Start the loop
for eachfile_input in *.bam
do
	echo $eachfile_input

	# pre-cleanup, only remove if exists
	[[ -d $DATA_DIR ]] && rm -r $DATA_DIR/

	# initialise
	mkdir -p $DATA_DIR/

	WASP_run_dir=$DATA_DIR/'WASP_run'
	mkdir -p $WASP_run_dir


	# set WASP calibration SNP path 	
#        SNP_DIR=$(echo $eachfile_input | \
#                sed 's/.\+\(CD_[0-9][0-9]\).*/\1/g')
#        SNP_DIR=WASP_calib/$SNP_DIR"_NEFM_neg_glut_unfiltered"

	SNP_DIR=$( echo $eachfile_input \
		| sed 's/\_trimmed/\t/g' \
		| cut -f 1 )
	SNP_DIR=../WASP_calib/$SNP_DIR

        echo $SNP_DIR

	final_output_file_name_prefix=${eachfile_input/%.bam/}

	$samtools index -@ 16 $eachfile_input

	# Pull out reads that need to be remapped to check for bias
	# Use the -p option for paired-end reads.
	#
	python $WASP/mapping/find_intersecting_snps.py \
		--is_paired_end \
		--is_sorted \
		--output_dir $WASP_run_dir/ \
		--snp_dir $SNP_DIR \
		$eachfile_input

	# Remap the reads, using same the program and options as before.
	# NOTE: If you use an option in the first mapping step that modifies the
	# reads (e.g. the -5 read trimming option to bowtie2) you should omit this
	# option during the second mapping step here (otherwise the reads will be modified
	# twice)!         
        echo "map for paired ends"
        bowtie2 -p 24 -X 2000 --mm --qc-filter --met 1 -t \
                --sensitive --no-mixed --no-discordant \
                -x $INDEX \
                -1 $WASP_run_dir/${eachfile_input/%.bam/.remap.fq1.gz} \
                -2 $WASP_run_dir/${eachfile_input/%.bam/.remap.fq2.gz} \
                | $samtools sort -l 9 -@ 20 -m 10G \
                -o $WASP_run_dir/${eachfile_input/%.bam/.remap_paired.bam}
	
        $samtools index $WASP_run_dir/${eachfile_input/%.bam/.remap_paired.bam}
	
	## remap to same position
        python $WASP/mapping/filter_remapped_reads.py \
                $WASP_run_dir/${eachfile_input/%.bam/.to.remap.bam} \
                $WASP_run_dir/${eachfile_input/%.bam/.remap_paired.bam} \
                $WASP_run_dir/${eachfile_input/%.bam/.remap_paired.keep.bam}

	#
	################################### merge, sort, and dedup
	#
	## Create a merged BAM containing [1] reads (PE+SE) that did
	## not need remapping [2] filtered remapped reads
	#	echo "merge bam"

       $samtools merge \
	       -f \
	       -@ 24 -l 5 \
                $WASP_run_dir/${eachfile_input/%.bam/.sim_pe_reads.merged.bam} \
                $WASP_run_dir/${eachfile_input/%.bam/.remap_paired.keep.bam} \
                $WASP_run_dir/${eachfile_input/%.bam/.keep.bam}

	#################### Sort and index the bam file
	####### !! Retain Uniquely mapped reads only (-v XA:Z:, SA:Z:)before send
	echo "filter and sort bam"
        $samtools view \
		-h \
		$WASP_run_dir/${eachfile_input/%.bam/.sim_pe_reads.merged.bam} \
	       | $samtools sort \
	       -l 9 -m 10G -@ 16 \
               -o $WASP_run_dir/${eachfile_input/%.bam/.sim_pe_reads.merged.sorted.bam}
	
	$samtools index -@ 20 \
		$WASP_run_dir/${eachfile_input/%.bam/.sim_pe_reads.merged.sorted.bam}

	### Dedup by WASP rmdup

	python $WASP/mapping/rmdup_pe.py \
		$WASP_run_dir/${eachfile_input/%.bam/.sim_pe_reads.merged.sorted.bam} \
		$DATA_DIR/merged_dedup.bam

	### !! remove the subsetted files in $DATA_DIR !!
	rm -r $WASP_run_dir

        $samtools sort -l 9 -m 10G -@ 16 \
                -o $DATA_DIR/merged_sorted_dedup.bam \
                $DATA_DIR/merged_dedup.bam

	$samtools index -@ 20 $DATA_DIR/merged_sorted_dedup.bam
	#
	############### Add SO:tag to header (as required by AddOrReplaceReadGroups
		# need to extract sam header first
		$samtools view -H \
			$DATA_DIR/merged_sorted_dedup.bam \
			| sed 's/^@HD/@HD\tVN:1\.5/g' \
			> $DATA_DIR/new_header.sam
	
		# replace the bam header
		$samtools reheader \
			$DATA_DIR/new_header.sam \
			$DATA_DIR/merged_sorted_dedup.bam \
			| $samtools sort \
			-@ 16 -m 5G -l 9 \
			-o $DATA_DIR/merged_sorted_dedup_header.bam
	
		$samtools index -@ 20 $DATA_DIR/merged_sorted_dedup_header.bam
	
	################### Strip all duplicate flags (0x400)

		$samtools view -h -S -@ 10 \
			--remove-flags 0x400 \
			$DATA_DIR/merged_sorted_dedup_header.bam \
			| $samtools sort \
			-@ 16 -m 5G -l 9 \
                        -o $DATA_DIR/merged_sorted_dedup_header_nodup.bam

	################### Re-add readgroup     
        rg_string=$final_output_file_name_prefix

        $samtools addreplacerg \
                -w \
                -r "ID:$rg_string\tLB:lib1\tPL:illumina\tSM:$rg_string\tPU:unit1" \
                -@ 16 \
                $DATA_DIR/merged_sorted_dedup_header_nodup.bam \
                | $samtools sort \
                -l 9 -m 10G -@ 16 \
                -o $OUTPUT_FINAL/$final_output_file_name_prefix'_bowtie2_barcode_WASPed.bam'


######## remove all intermediate bam files
#        	rm -r $DATA_DIR/*.bam
#        	rm -r $DATA_DIR/*.bai

#		rm -r $DATA_DIR/barcode_bams/$1*
#		rm -r $DATA_DIR/WASP_run/$1*
		
#		rm -r $DATA_DIR
################################## clean up the temp folder

	echo '############################'

done

[[ -d $DATA_DIR ]] && rm -r $DATA_DIR/

unset OMP_NUM_THREADS
