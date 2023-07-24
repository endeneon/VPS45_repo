#!/bin/bash

# Siwei 13 Sept 2021

# count both total and unique reads of each .bam file

rm read_sum.txt

for eachfile in *.bam
do
	echo $eachfile
	printf "$eachfile\t" >> read_sum.txt
	samtools view -@ 40 $eachfile | wc -l >> read_sum.txt
	# move back to the end of the previous line
#	echo -e '\e[1A'
	# print a second tab
#	printf "\t" >> read_sum.txt
#	samtools view \
#		-@ 40 \
#		-F 4 \
#		$eachfile \
#		| wc -l >> read_sum.txt
done

