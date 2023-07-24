#!/bin/bash

# 03 May 2023
# Siwei

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

# 13 Jul 2023
# use gatk 4261
# investigate the missing rs2027349 issue

vcf_suffix="_995.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

temp_dir="/home/zhangs3/NVME/package_temp/vcf_995_R21_2"

mkdir -p vcf_output_995 

# date > log.log

rm -r $temp_dir

for EACHFILE in *.bam
do

#	rm *.vcf
#	rm *.vcf.idx
#	rm *.recal
#	rm *.tranches
#	rm *.R
	mkdir -p $temp_dir

	echo $EACHFILE
	echo $EACHFILE >> $temp_dir"/log.log"
	date >> $temp_dir"/log.log"
############## call variants ##############
	
	samtools index -@ 20 $EACHFILE
	
	$gatk4 --java-options "-Xmx150g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 20 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $EACHFILE \
		-O $temp_dir"/raw_variants.vcf"

############variants filtering##########
############SNP#########################
############Remove previous calibration files######

############ SNP Recalibration ###########

        $gatk4 --java-options "-Xmx150g" \
	 	VariantRecalibrator \
		-R $ref_genome \
		-V $temp_dir"/raw_variants.vcf" \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ref_path/hapmap_3.3.hg38.vcf \
		--resource:omni,known=false,training=true,truth=true,prior=12.0 $ref_path/1000G_omni2.5.hg38.vcf \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 $ref_path/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
		--max-gaussians 4 \
		-an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
		-mode SNP \
		-AS \
		-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 90.0 \
		-O $temp_dir/raw_recalibrate_SNP.recal \
		--tranches-file $temp_dir/raw_recalibrate_SNP.tranches \
		--rscript-file $temp_dir/raw_recalibrate_SNP_plots.R  

########### apply SNP recalib###########

        $gatk4 --java-options "-Xmx150g" \
		ApplyVQSR \
		-R $ref_genome \
		-V $temp_dir/raw_variants.vcf \
		-mode SNP \
		-AS \
		--truth-sensitivity-filter-level 99.5 \
		--recal-file $temp_dir/raw_recalibrate_SNP.recal \
		--tranches-file $temp_dir/raw_recalibrate_SNP.tranches \
		-O $temp_dir/variants_calibrated.vcf

########## filter SNPs located in the hg38 blacklisted region #######

	$gatk4 --java-options "-Xmx150g" \
		VariantFiltration \
		-R $ref_genome \
		-V $temp_dir/variants_calibrated.vcf \
		--mask $ref_path/hg38_blacklisted_regions.bed --mask-name "blacklisted_regions" \
		--filter-name "PASS" \
		--filter-expression "AC > 0" \
		-O $temp_dir/variants_calibrated_blacklisted.vcf

########## select SNP for output from autosomes only #######

	$gatk4 --java-options "-Xmx150g" \
		SelectVariants \
		-R $ref_genome \
		-V $temp_dir/variants_calibrated_blacklisted.vcf \
		--select-type-to-include SNP \
		--restrict-alleles-to BIALLELIC \
		-L $ref_path/autosomes.list \
		--exclude-filtered true \
		-O $temp_dir/pre_output.vcf

######## extract het only

        cat $temp_dir/pre_output.vcf \
                | grep "^#" \
                > $temp_dir/header.txt

        cat $temp_dir/pre_output.vcf \
                | grep -v "^#" \
                | grep "0/1" \
                > $temp_dir/body.txt

        cat $temp_dir/header.txt \
		$temp_dir/body.txt \
                > vcf_output_995/${EACHFILE/%.bam/}$vcf_suffix

###### cleanup

#	rm -r $temp_dir

done

