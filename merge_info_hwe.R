### Usage Rscript merge_info_hwe.R {file1}.info {file2}.hwe {file3}.frq {outfile}

args <- commandArgs(TRUE)

info_file = read.delim(args[1], header = TRUE, stringsAsFactors = FALSE)
hwe_file = read.delim(args[2], header = TRUE, stringsAsFactors = FALSE)
freq_file = read.delim(args[3], header = TRUE, stringsAsFactors = FALSE)

out_merged = merge(info_file, hwe_file, by.x = "ID", by.y = "SNP", all = TRUE)
out_merged = merge(out_merged, freq_file, by.x = "ID", by.y = "SNP", all = TRUE)
out_merged = out_merged[,c("CHROM","POS","ID","REF","ALT","A1.x","A2.x","NCHROBS","GENO","MAF","INFO_Info","P")]
names(out_merged)[6:12] = c("A1","A2","N*2","GENO","MAF","INFO","HWE.P")

out_name = paste(args[4], ".info.hwe.frq", sep = "")

write.table(out_merged, file = out_name, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

cat(paste(out_name, " successfully generated.\n", sep = ""))

