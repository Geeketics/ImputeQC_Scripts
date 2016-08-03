### Rscript merge_info_hwe.R {file1}.info {file2}.hwe {outfile}

args <- commandArgs(TRUE)

info_file = read.delim(args[1], header = TRUE, stringsAsFactors = FALSE)
hwe_file = read.delim(args[2], header = TRUE, stringsAsFactors = FALSE)

out_merged = merge(info_file, hwe_file, by.x = "ID", by.y = "SNP", all = TRUE)
out_merged = out_merged[,c("CHROM","POS","ID","REF","ALT","A1","A2","GENO","INFO_Info","P")]
names(out_merged)[9:10] = c("INFO","HWE.P")

out_name = paste(args[3], ".info.hwe", sep = "")

write.table(out_merged, file = out_name, quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

cat(paste(out_name, " successfully generated.\n", sep = ""))