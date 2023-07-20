#!/bin/bash

# 03 May 2023
# Siwei

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

vcf_suffix="_het.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

mkdir -p het_vcf 

date > log.log

rm *.recal
rm *.tranches
rm *.R


for EACHFILE in *.vcf
do


########## filter SNPs located in the hg38 blacklisted region #######

	$gatk4 --java-options "-Xmx150g" \
		VariantFiltration \
		-R $ref_genome \
		-V $EACHFILE \
		--mask $ref_path/hg38_blacklisted_regions.bed --mask-name "blacklisted_regions" \
		--filter-expression "DP > 0" --filter-name "PASS" \
		-O variants_calibrated_blacklisted.vcf

########## select SNP for output from autosomes only #######

	$gatk4 --java-options "-Xmx150g" \
		SelectVariants \
		-R $ref_genome \
		-V variants_calibrated_blacklisted.vcf \
		--select-type-to-include SNP \
		--restrict-alleles-to BIALLELIC \
		-L $ref_path/autosomes.list \
		--exclude-filtered true \
		-O pre_output.vcf

######## extract het only

        cat pre_output.vcf \
                | grep "^#" \
                > header.txt

        cat pre_output.vcf \
                | grep -v "^#" \
                | grep "0/1" \
                > body.txt

        cat header.txt body.txt \
                > het_vcf/${EACHFILE/%.vcf/}$vcf_suffix

###### cleanup

        rm *.recal
        rm *.tranches
        rm *.R
	rm *.txt

done

