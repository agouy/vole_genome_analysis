### RNDsp plots
source("./scripts/Resequencing-data-analysis-functions.R")
require(GenomicRanges)

load("./data/annotation.rda")
load("./data/pi_rndsp.rda")

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%m%d_%H%M%S")

span.loess <- 0.01

vole.col <- c(rep("#228B22", 4), RColorBrewer::brewer.pal(9, "Set1"))
samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")

rndsps <- stats[, grep("rndsp", colnames(stats))]

ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)

# 3 x 4 individuals
## arvalis
# pdf(file = paste0(output.dir, stamp, "_RNDsp.pdf"), paper = "a4", width = 10, height = 8)
png(file = paste0(output.dir, stamp, "_RNDsp.png"), width = 17, height = 22, res = 600, units = "cm")

par(mfrow=c(5,1), mar = c(1,5,1,0), omi = c(0.8,0,0,0))
for(i in c(1:4, 8)) {
  Y <- rndsps[ind.ord][[i]]

  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0, ifelse(i == 8, 0.6, 0.4)), scale = ifelse(i == 1, 0.3, 0.4), clust = 100, ylab = b)
  if(i==8) {
    segments(1900, 0.4, 2000, 0.4, lwd = 2)
    text(1950, 0.4, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)
dev.off()

png(file = paste0(output.dir, stamp, "_RNDsp_supp1.png"), width = 10, height = 8, res = 600, units = "in")

## close
par(mfrow=c(3,1), mar = c(1,5,1,0), omi = c(0.8,0,0,0))

for(i in 5:7) {
  Y <- rndsps[ind.ord][[i]]

  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0,0.5), scale = 0.5, clust = 100, ylab = b)
  if(i==7) {
    segments(1900, 0.25, 2000, 0.25, lwd = 2)
    text(1950, 0.25, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)
dev.off()

png(file = paste0(output.dir, stamp, "_RNDsp_supp2.png"), width = 10, height = 8, res = 600, units = "in")

## far
par(mfrow=c(4,1), mar = c(1,5,1,0), omi = c(0.8,0,0,0))

for(i in 9:12) {
  Y <- rndsps[ind.ord][[i]]

  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0,0.8), scale = 0.5, clust = 100, ylab = b)
  if(i==12) {
    segments(1900, 0.25, 2000, 0.25, lwd = 2)
    text(1950, 0.25, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)

dev.off()
