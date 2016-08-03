#!/bin/sh

infile=$1 # to edit (just the prefix)
outfile=$2 # to save (just the prefix)

zcat ${infile}.vcf.gz | cut -f1-8 | grep -v "^##" | tr -s ' ' | tr ' ' '\t' > ${outfile}.temp

echo CHROM POS ID REF ALT QUAL FILTER INFO_RefPanelAF INFO_AN INFO_AC INFO_Info | tr ' ' '\t' > ${outfile}.info

tail -n +2 ${outfile}.temp | tr '=' '\t' | tr ';' '\t' | cut -f1-7,9,11,13,15 >> ${outfile}.info

rm ${outfile}.temp

