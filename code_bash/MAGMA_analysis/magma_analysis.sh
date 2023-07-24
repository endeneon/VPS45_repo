#!/bin/bash

# Siwei 14 Sept 2020

for eachfile in *.raw
do
	magma --gene-results $eachfile \
		--set-annot VPS45_gene_set_4_magma.txt col=2,1 \
		--out output/${eachfile/%hg19_SETD1A_organoid_line18_/}
done

