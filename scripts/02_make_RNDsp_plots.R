### RNDsp plots
load("./data/13V-22sc-Di-Nomiss-GT.rda")
source("./scripts/Resequencing-data-analysis-functions.R")
require(GenomicRanges)

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%M%d_%H%M%S")

date()
?timestamp
span.loess <- 0.01

vole.col <- c(rep("#228B22", 4), RColorBrewer::brewer.pal(9, "Set1"))

gene.positions <- read.table("./data/marv_genes_augustus.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

homologs <- read.table("./data/homology-marv-mmus.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

homologs <- homologs[order(homologs$eval), ]
homologs <- homologs[!duplicated(homologs$id_arv), ]

rownames(gene.positions) <- gene.positions$id
gene.positions$id_mouse <- NA
gene.positions[homologs$id_arv,]$id_mouse <- homologs$id_mouse

annot <- gene.positions[!is.na(gene.positions$id_mouse),]

win.size <- 0
annot.gr <- GRanges(
  seqnames=annot$sc,
  ranges = IRanges(start=annot$st - win.size, end=annot$en + win.size),
  id=annot$id_mouse
)

map.table <- read.table("./data/uniprot-to-name.txt", header = TRUE, sep = "\t")
map.table <- map.table[!duplicated(map.table$From), ]
rownames(map.table) <- map.table$From

annot.gr$name <- as.character(map.table[annot.gr$id, ]$To)


# get input data
dxys <- res[, grep("dxy", colnames(res))]
av.dxys <- colMeans(dxys, na.rm = TRUE)

myod <- dxys[, grep("pop9", colnames(dxys))]

samp.names <- c("M. arvalis (W)","M. arvalis (I)","M. arvalis (E)","M. arvalis (C)",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")

dspec <- list()
for(i in as.character(1:13)[-9]) {
  popu <- i
  combs <- do.call("rbind", strsplit(gsub("dxy.pop", "", names(av.dxys)), ".pop")) 
  ids <- apply(combs, 1, function(x) {
    any(x == popu)
  })
  dao <- names(which.max(av.dxys[ids]))
  dab <- names(which.min(av.dxys[ids]))
  
  n.dbo <- unlist(strsplit(gsub("dxy.pop", "", c(dao, dab)), ".pop"))
  n.dbo <- n.dbo[n.dbo != popu]
  dbo <- c(paste0("dxy.pop", n.dbo[1], ".pop", n.dbo[2]),
           paste0("dxy.pop", n.dbo[2], ".pop", n.dbo[1]))
  
  da <- (dxys[, dab] + dxys[, dao] - dxys[, colnames(dxys) %in% dbo]) / 2
  dspec[[i]] <- da
}
names(dspec) <- paste0("dxysp.", names(dspec))

d.av.myo <- rowMeans(myod)
rndsps <- lapply(dspec, function(x) {x / d.av.myo})
cond <- d.av.myo <= quantile(d.av.myo, 0.025)
# cond <- d.av.myo <= quantile(d.av.myo, 0.0)
ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)

# 3 x 4 individuals
## arvalis
pdf(file = paste0(output.dir, stamp, "_RNDsp.pdf"), paper = "a4", width = 10, height = 8)
par(mfrow=c(4,1), mar = c(1,5,1,0), omi = c(0.8,0,0,0))
for(i in 1:4) {
  Y <- rndsps[ind.ord][[i]]
  Y[cond] <- NA
  
  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0,0.4), scale = ifelse(i == 1, 0.3, 0.4), clust = 100, ylab = b)
  if(i==4) {
    segments(1900, 0.25, 2000, 0.25, lwd = 2)
    text(1950, 0.25, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)

## close
for(i in 5:8) {
  Y <- rndsps[ind.ord][[i]]
  Y[cond] <- NA
  
  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0,0.5), scale = 0.5, clust = 100, ylab = b)
  if(i==4) {
    segments(1900, 0.25, 2000, 0.25, lwd = 2)
    text(1950, 0.25, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)

## far
for(i in 9:12) {
  Y <- rndsps[ind.ord][[i]]
  Y[cond] <- NA
  
  a <- samp.names[-9][ind.ord][i]
  b <- bquote(paste('RND'['sp'], " (", .(a), ")"))
  
  plotOut(Y, ylim = c(0,0.8), scale = 0.5, clust = 100, ylab = b)
  if(i==4) {
    segments(1900, 0.25, 2000, 0.25, lwd = 2)
    text(1950, 0.25, "100 Mb",pos = 1,offset = 0.5, adj = 0.5)
  }
}
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.7)
title(xlab = "Scaffold", outer = TRUE)

dev.off()
