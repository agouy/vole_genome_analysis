## new figure

source("./scripts/00_utils.R")
require(GenomicRanges)

load("./data/annotation.rda")
load("./data/pi_rndsp.rda")


span.loess <- 0.01

vole.col <- c(rep("#228B22", 4), RColorBrewer::brewer.pal(9, "Set1"))

samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")

rndsps <- stats[, grep("rndsp", colnames(stats))]

getGeneList(rndsps[[1]], stats, annot.gr, 0.9995)
getGeneList(rndsps[[2]], stats, annot.gr, 0.9995)
getGeneList(rndsps[[3]], stats, annot.gr, 0.9995)
getGeneList(rndsps[[4]], stats, annot.gr, 0.9999)
# getGeneList(rndsps[[10]], stats, annot.gr, 0.9999)

#### ORs

plotDivDiv <- function(stats = stats, scaffold = "ScOZjSD_3553", center = 30533341, main = "", leg=FALSE, winsize=4e5) {
  toplo <- stats[stats$sc == scaffold & stats$st > center - winsize & stats$en < center + winsize, ]
  toplo <- toplo[!(rowSums(is.na(toplo)) > 0), ]
  
  X.val <- (toplo$st+toplo$en)/2/1e6
  par(mar = c(3,5,2,2))
  matplot(
    X.val,
    toplo[, grep("pi.", colnames(toplo))],
    pch = 16, 
    type = "o", 
    cex = 0.6, 
    lty = 1,
    bty = "l",
    lwd = 1,
    main = "",
    ylab = "Diversity >",
    xlab = "",
    col = couls
  )
  abline(v=center/1e6, lty = 2, col = "grey")
  
  title(main, adj = 0)
  par(mar = c(7,5,0,2))
  matplot(
    X.val,
    -toplo[, grep("rndsp.", colnames(toplo))],
    pch = 16, 
    cex = 0.6,
    bty = "n",
    type = "o", 
    lty = 1, 
    lwd = 1,
    main = "",
    ylab = "< Divergence",
    xlab = "",
    col = couls,
    xaxt="n",
    yaxt = "n"
  )
  as.val <- pretty(-unlist(toplo[, grep("rndsp.", colnames(toplo))]))
  as.lab <- -as.val
  up.val <- pretty(X.val)
  axis(2, at=as.val,labels = as.lab)
  axis(3, at=up.val,labels = rep("", length(up.val)))
  segments(center/1e6, 0, center/1e6, -0.5, lty = 2, col = "grey")
  
  ov <- subsetByOverlaps(
    annot.gr,
    GRanges(seqnames = scaffold, ranges = IRanges(start = center - 2e5, end = center + 2e5))
  )
  ov <- ov[lengths(ov) > 1e3]
  par(xpd=TRUE)
  segments(start(ov)/1e6, min(as.val), end(ov)/1e6, min(as.val), lwd = 2, col = "slategrey")
  text((start(ov)+end(ov))/2/1e6, 1.02*min(as.val), ov$name, pos=2, cex = 0.9, srt=90, offset = c(0,-0.))
  par(xpd=FALSE)
  
  if(leg) {
    par(xpd=TRUE)
    legend("bottomright", cex = 1.,
           inset=c(-0.05,-1.), legend=samp.names[species], lty = 1, pch = 16, col = couls) 
    par(xpd=FALSE)
  }
}


output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%m%d_%H%M%S")

# pdf(file = paste0(output.dir, stamp, "_Zooms_v2.pdf"), width = 10, height = 5)
png(file = paste0(output.dir, stamp, "_Zooms_v3.png"), width = 9, height = 7, res = 600, units = "in")
layout(matrix(c(1,2,5,6,3,4,7,8), nrow = 4))

species <- c(1:4, 11)
# couls <- c("#228B22","chartreuse1","chartreuse2","lightgreen","#984EA3")
couls <- RColorBrewer::brewer.pal(5, "Paired")
columns <- c("sc", "st", "en", paste0("pi.", species), paste0("rndsp.", species))

annot.gr[annot.gr$name%in%c("Olfr140"), ]
annot.gr[annot.gr$name%in%c("Olfr1019"), ]

annot.gr[annot.gr$name%in%c("H2-Eb1","H2-Ea"), ]
annot.gr[annot.gr$name %in% c("Cd200r1l"), ]

# Order: same as main text
plotDivDiv(stats[, columns], leg = TRUE,
           scaffold = "ScOZjSD_3553", center = 166250000, main = "a. Olfr140 Region")

plotDivDiv(stats[, columns],
           scaffold = "ScOZjSD_3553", center = 30500000, main = "b. Olfr1019 Region")

plotDivDiv(stats[,columns],
           scaffold = "ScOZjSD_2657", center = 17650000, main = "c. Histocompatibility-2 Cluster Region")

plotDivDiv(stats[,columns],
           scaffold = "ScOZjSD_3553", center = 156650000, main = "d. CD200 Region")


dev.off()


