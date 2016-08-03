#!/bin/sh

infile=$1 # to edit (just the prefix)
outfile=$2 # to save (just the prefix)

## added compensation for ARIC impute lines with "TYPED" in INFO set
zcat ${infile}.vcf.gz | cut -f1-8 | grep -v "^##" | tr -s ' ' | tr ' ' '\t' > ${outfile}.temp

echo CHROM POS ID REF ALT QUAL FILTER INFO_RefPanelAF INFO_AN INFO_AC INFO_Info | tr ' ' '\t' > ${outfile}.info

tail -n +2 ${outfile}.temp | grep -v TYPED | tr '=' '\t' | tr ';' '\t' | cut -f1-7,9,11,13,15 >> ${outfile}.info
tail -n +2 ${outfile}.temp | grep TYPED | tr '=' '\t' | tr ';' '\t' | cut -f1-7,9,11,13,16 >> ${outfile}.info

rm ${outfile}.temp


