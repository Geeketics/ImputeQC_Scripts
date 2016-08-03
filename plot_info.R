### Usage: Rscript plot_info.R {file}.info.hwe.frq {outfile}

args <- commandArgs(TRUE)

file = read.delim(args[1], header = TRUE, stringsAsFactors = FALSE)

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], ".info.frq.jpeg", sep = ""))
plot(x = file$MAF, y = file$INFO, xlim = c(0,0.5), ylim = c(0,1), xlab = "Minor Allele Frequency", ylab = "Imputation Info Score", pch = 4)
dev.off()

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], "_hist.info.frq.jpeg", sep = ""))
hist(x = file$MAF, y = file$INFO, xlim = c(0,0.5), ylim = c(0,1), xlab = "Minor Allele Frequency", ylab = "Imputation Info Score", pch = 4)
dev.off()

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], ".info.hwe.jpeg", sep = ""))
plot(x = -log10(file$HWE.P), y = file$INFO, ylim = c(0,1), xlab = expression("Hardy Weinberg"~~-log[10](italic(P))), ylab = "Imputation Info Score", pch = 4)
dev.off()

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], "_hist.info.hwe.jpeg", sep = ""))
plot(x = -log10(file$HWE.P), y = file$INFO, ylim = c(0,1), xlab = expression("Hardy Weinberg"~~-log[10](italic(P))), ylab = "Imputation Info Score", pch = 4)
dev.off()

cat(paste(args[2], "info.frq.jpeg, ", args[2] "_hist.info.frq.jpeg, ", args[2], ".info.hwe.jpeg, and ", args[2], "_hist.info.hwe.jpeg were successfully generated\n"))

