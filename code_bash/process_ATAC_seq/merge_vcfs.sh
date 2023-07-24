#!/bin/bash


# 24 May 2022
# 07 May 2021 Siwei

## use conda aligners environment

## Set these environment vars to point to
## your local installation

gatk37_path="/home/zhangs3/Data/Tools/GATK37"
gatk36_path="/home/zhangs3/Data/Tools/gatk_36/opt/gatk-3.6"
gatk4_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1"

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin"
openjdk_16="/home/zhangs3/Data/Tools/jdk-16.0.1/bin/java"

hg38_index="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

rm -r output
mkdir -p output

## send all .vcf files to a .list file
ls *.vcf > vcf_list.list

$jre_8/java -Xmx200G -jar $gatk36_path/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R $hg38_index \
        --variant vcf_list.list \
        -genotypeMergeOptions UNIQUIFY \
        -o output/raw_merged_vcfs.vcf

## select variants that only includes SNP
$gatk4_path/gatk SelectVariants \
        -R $hg38_index \
        -V output/raw_merged_vcfs.vcf \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        -XL /home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed \
	-O output/"trimmed_merged_SNP.vcf"

## select variants with DP >= 2
$gatk4_path/gatk VariantFiltration \
        -R $hg38_index \
        -V output/"trimmed_merged_SNP.vcf" \
        --filter-name "DP_filter" -filter "DP < 1" \
        -O output/"NGN2_22Jun2023_trimmed_merged_SNP_VF.vcf"

## export vcfs to table
$gatk4_path/gatk VariantsToTable \
        -V output/"NGN2_22Jun2023_trimmed_merged_SNP_VF.vcf" \
        -F CHROM -F POS -F ID -F REF -F ALT -F NCALLED \
        -GF AD \
        -O output/"trimmed_merged_SNP_VTT.txt"

## convert the txt title line
cat output/"trimmed_merged_SNP_VTT.txt" \
        | sed 's/NA/0,0/g' \
        | sed 's/,/\t/g' \
        | sed 's/variant\.AD/variant\tAD/g' \
        | sed 's/variant.\.AD/variant\tAD/g' \
        | sed 's/variant..\.AD/variant\tAD/g' \
        > output/"NGN2_22JunMay2023_trimmed_scATAC_noBQSR_4_R.txt"

