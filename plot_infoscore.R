## Usage: Rscript plot_infoscore.R {infile}.info.hwe.frq {outfile}

## To plot the average info score across allele frequency bins and hardy weinberg bins - saves summary info too

args <- commandArgs(TRUE)

datafile = read.delim(args[1], header = TRUE, stringsAsFactors = FALSE)

## Minor Allele Frequency Plot
datafile = within(datafile, MAFbins <- cut(MAF, seq(0, 0.5, 0.01), include.lowest = FALSE))
datafile$MAFbins = do.call(rbind,strsplit(as.character(datafile$MAFbins),c(",")))[,2]
datafile$MAFbins = do.call(rbind,strsplit(as.character(datafile$MAFbins),c("]")))[,1]

datasummaryMAF = aggregate(datafile$INFO ~ datafile$MAFbins, FUN = function(x) c("Mean" = mean(x)))
names(datasummaryMAF) = c("MAFbins","Mean_Info")
datasummaryMAF$Count = aggregate(datafile$INFO ~ datafile$MAFbins, FUN = function(x) c("Count" = table(!is.na(x))))[,2]

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], ".info.maf.jpeg", sep = ""))
plot(x = datasummaryMAF$MAFbins, y = datasummaryMAF$Mean_Info, axes = FALSE, xlab = "Minor Allele Frequency", ylab = "Imputation Info Score (Mean per MAF bin)", pch = 4, type = "b", xlim = c(0,0.5), ylim = c(0,1))
axis(side = 1, labels = round(seq(0,0.5, 0.05),2), at = seq(0,0.5, 0.05))
axis(side = 2, labels = round(seq(0,1,0.1),1), at = seq(0,1,0.1))
dev.off()

write.table(datasummaryMAF, file = paste(args[2], "_summary.info.maf", sep = ""), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

## Hardy Weinberg Plot
datafile = within(datafile, HWEbins <- cut(HWE.P, seq(0, 1, 0.01), include.lowest = FALSE))
datafile$HWEbins = do.call(rbind,strsplit(as.character(datafile$HWEbins),c(",")))[,2]
datafile$HWEbins = do.call(rbind,strsplit(as.character(datafile$HWEbins),c("]")))[,1]

datasummaryHWE = aggregate(datafile$INFO ~ datafile$HWEbins, FUN = function(x) c("Mean" = mean(x)))
names(datasummaryHWE) = c("HWEbins","Mean_Info")
datasummaryHWE$Count = aggregate(datafile$INFO ~ datafile$HWEbins, FUN = function(x) c("Count" = table(!is.na(x))))[,2]

jpeg(width = 6, height = 6, units = "in", res = 300, file = paste(args[2], ".info.hwe.jpeg", sep = ""))
plot(x = datasummaryHWE$HWEbins, y = datasummaryHWE$Mean_Info, axes = FALSE, xlab = "Hardy-Weinberg Equilibrium", ylab = "Imputation Info Score (Mean per HWE bin)", pch = 4, type = "b", xlim = c(0,1), ylim = c(0,1))
axis(side = 1, labels = round(seq(0, 1, 0.05),2), at = seq(0, 1, 0.05))
axis(side = 2, labels = round(seq(0, 1, 0.1),1), at = seq(0, 1, 0.1))
dev.off()

write.table(datasummaryHWE, file = paste(args[2], "_summary.info.hwe", sep = ""), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
