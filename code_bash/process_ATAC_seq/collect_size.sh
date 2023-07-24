#!/bin/bash

# Siwei 29 Jun 2021

# Use Picard from GATK 4.1.8.1

picard_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1/picard2_25_4.jar"

for eachfile in *.bam
do
	java -jar $picard_path CollectInsertSizeMetrics \
		-I $eachfile \
		-O size_dist/${eachfile/%_new_WASPed\.bam}.txt \
		-H size_dist/${eachfile/%_new_WASPed\.bam}.pdf \
		-M 0.5
done


