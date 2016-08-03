Strand fixing file: new_v4_CompareRS100GWithBimFiles.r

The script requires a Plink bim file and 1000Genomes legend files. The script is run per chromosome. The script generates:
- lists of SNPs to be kept (filename.chr)
- SNPs that need to be updated by position (filename.pos)
- SNPs to be deleted (filename.delete) 
- SNPs to be flipped. (filename.flip)

The files (per chromosome) are then merged and used in the first step of imputation to generated strand-fixed Plink files. The first step of the imputation pipeline is in the file: 1_strand_align.sh
