#!/bin/sh

infile=$1 #.vcf.gz to edit
outfile=$2 #.vcf.gz to save to

zcat ${infile} | awk -v OFS="\t" '{if($3 == "."){$3 = $1"_"$2"_"$4"_"$5}{print}}' | bgzip > ${outfile}

tabix -p vcf ${outfile}

echo New vcf file called ${outfile} and index file called ${outfile}.tbi created.\n

