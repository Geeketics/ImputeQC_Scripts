#!/bin/sh

vcf=$1 #.vcf.gz to work with
subset=$2 #list of samples to subset by
out=$3 #outfile prefix

## calculate HWE and Allele Freq for all SNPs - subsetted by subset list
plink2 --vcf ${vcf} --keep ${subset} --hardy --freq --out ${out}

## save a list of SNPs with extreme HWE p-values
cat ${out}.hwe | tr -s ' ' | tr ' ' '\t' | sed 's/^\t//g' | sed 's/\t$//g' > ${out}.hwe.tab
awk -v OFS="\t" '{if($9 < 1e-10){print}}' ${out}.hwe.tab > ${out}.hwe.out

## count how many SNPs had extreme HWE p-values
wc -l ${out}.hwe.out

echo Files called ${out}.frq, ${out}.hwe, ${out}.hwe.tab, and ${out}.hwe.out have been created.
echo ${out}.frq - PLINK output from the --freq flag
echo ${out}.hwe - PLINK output from the --hardy flag
echo ${out}.hwe.tab - a tab-separated version of ${out}.hwe
echo ${out}.hwe.out - a list of SNPs with a HWE p-value less than 1e-10

