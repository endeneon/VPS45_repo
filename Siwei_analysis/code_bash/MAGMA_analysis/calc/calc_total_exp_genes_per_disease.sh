#!/bin/bash
# Siwei 15 Sept 2020

rm record.txt

for eachfile in ../*.raw
do	
	echo $eachfile
	echo -e "$eachfile\t\c" >> record.txt
	cat $eachfile | grep -v "^#" | cut -f 1 -d ' ' | sort > disease_genes.txt
	join -1 2 -2 1 VPS45_all_exp_genes_name_sorted.txt disease_genes.txt | wc -l >> record.txt
done

