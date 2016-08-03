### Usage: Rscript imputationinfo.R filename

args <- commandArgs(TRUE)

infoscore = read.delim(args[1], stringsAsFactors = FALSE, header = TRUE)
names(infoscore) = c("CHR","POS","SNP","INFO")

## Make manhattan.generic function
manhattan.generic = function(x, chr = NULL, bp = NULL, value = NULL, ylabel = NULL, col = c("gray10","gray60"), suggestiveline = NA, genomewideline = NA){
  CHR = BP = VALUE = index = NULL
  data = data.frame(CHR = x[[chr]], BP = x[[bp]], VALUE = x[[value]])
  data = na.omit(data)
  data = data[order(data$CHR, data$BP),]
  data$index = data$pos = NA
  
  indexing = 0
  for(i in unique(data$CHR)){
    indexing = indexing + 1
    data$index[data$CHR == i] = indexing
  }
  
  lastbase = 0
  tickpos = NULL
  data$BPnew = data$BP/1e+06
  for(i in unique(data$index)){
    if(i == 1){
      data$pos[data$index == i] = data$BPnew[data$index == i]
    }
    else{
      lastbase = lastbase + max(data$BPnew[data$index == (i - 1)])
      data$pos[data$index == i] = data$BPnew[data$index == i] + lastbase
    }
    tickpos = c(tickpos, (min(data$pos[data$index == i]) + max(data$pos[data$index == i]))/2 + 1)
  }
  xlabel = "Chromosome"
  label = unique(data$CHR)
  
  xmin = (max(data$pos) * -0.03)
  xmax = (max(data$pos) * 1.03)
  ymax = (max(data$VALUE) * 1.03)
  ymin = (max(data$VALUE) * -0.03)
  yinterval = ((ymax - ymin)/4)
  
  plot(0, bty = "n", pch = "", xaxt = "n", xlim = c(xmin,xmax), ylim = c(ymin,ymax), xlab = xlabel, ylab = ylabel)
  axis(side = 1, at = tickpos, labels = label, lwd = 1, lwd.ticks = 1, line = 0.25)
  for(i in unique(data$CHR)){
    if((i %% 2) == 0){
      points(x = data$pos[data$CHR == i], y = data$VALUE[data$CHR == i], col = col[1], pch = 20)
    }
    else{
      points(x = data$pos[data$CHR == i], y = data$VALUE[data$CHR == i], col = col[2], pch = 20)
    }
  }
  if(!is.na(suggestiveline)){
    abline(h = suggestiveline, col = "blue")
  }
  if(!is.na(genomewideline)){
    abline(h = genomewideline, col = "red")
  }
}

manhattanname = paste(args[1], "manhattan", "jpeg", sep = ".")
histogramname = paste(args[1], "histogram", "jpeg", sep = ".")

jpeg(file = manhattanname, width = 1600, height = 800)
manhattan.generic(x = infoscore, chr = "CHR", bp = "POS", value = "INFO", ylabel = "Info Score", suggestiveline = 0.3, genomewideline = 0.8)
dev.off()

jpeg(file = histogramname, width = 800, height = 800)
hist(infoscore$INFO, density = TRUE, main = "Imputation Info Score", xlab = "Info Score")
abline(v = c(0.3,0.8), col = c("blue","red"))
dev.off()

cat(paste(manhattanname, " and ", histogramname, " files succesfully generated.\n", sep = ""))
