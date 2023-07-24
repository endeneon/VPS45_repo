#!bin/bash

# Siwei 9 Aug 2021

# Siwei 3 Sept 2018
# Siwei 12 Feb 2021
# make adaptation to genres01 and STAR 2.7.7

rm -r ~/NVME/temp_cache_STAR

date > stat.txt

for EACHFILE_L1 in *_1.fq.gz
do
	echo $EACHFILE_L1
	echo ${EACHFILE_L1/%_1.fq.gz/_2.fq.gz}
	echo $EACHFILE_L1 >> stat.txt
	####init
        #date
        rm -rf ~/NVME/temp_cache_STAR
	#mkdir ~/4TB_2/RNASeq_22May2019_C202SC19010639/raw_data/STAR_output/${EACHFILE_L1/%_1.fq.gz/}

	STAR --runThreadN 30 \
        --genomeDir ~/Data/Databases/Genomes/hg38/STAR_db \
        --readFilesIn $EACHFILE_L1 ${EACHFILE_L1/%_1.fq.gz/_2.fq.gz} \
        --readFilesCommand zcat \
        --outFileNamePrefix /home/zhangs3/Data/FASTQ/20_rapid_neuron_raw_data/STAR_output/${EACHFILE_L1/%_1.fq.gz/} \
	--bamRemoveDuplicatesType UniqueIdentical \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile ~/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf \
        --quantMode GeneCounts \
        --outSAMattrIHstart 0 \
        --outSAMstrandField intronMotif \
        --outSAMmultNmax 1 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outBAMcompression 10 \
        --outBAMsortingThreadN 40 \
        --outBAMsortingBinsN 40 \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 2 \
        --outSJfilterReads Unique \
        --limitBAMsortRAM 40000000000 \
        --alignSoftClipAtReferenceEnds No \
        --quantTranscriptomeBAMcompression 10 10 \
        --outFilterScoreMinOverLread 0.30 \
        --outFilterMatchNminOverLread 0.30 \
        --outTmpDir ~/NVME/temp_cache_STAR
	#--twopassMode Basic

done

