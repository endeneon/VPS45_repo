#!/bin/bash

# Siwei 10 May 2022
# use the improved method as in 
# https://www.biorxiv.org/content/10.1101/496521v1

macs2 callpeak \
	-t *.bam \
	-f BAMPE \
	-g 2.7e9 \
	-q 0.05 \
	--keep-dup all \
	--nolambda \
	--min-length 100 \
	--max-gap 50 \
	--buffer-size 1000000 \
	-n output/$1 \
	--seed 42
