## new figure

source("./scripts/Resequencing-data-analysis-functions.R")

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%M%d_%H%M%S")

getGeneList(rndsps[[1]], res, annot.gr[,], quanti = 0.9999)


str(res)

getGeneDiv <- function(test.var, res, annot.gr, quanti = 0.99, id.samp) {
  pis <- res[,grep("pi", colnames(res))]
  out.gr <- GRanges(
    seqnames = res$sc,
    ranges = IRanges(start = res$st, end = res$en),
    stat = test.var,
    pi = pis[,id.samp]
  )
  ovL <- findOverlaps(annot.gr, out.gr)
  annot.all <- annot.gr[queryHits(ovL), ]
  
  annot.all$stat <- out.gr[subjectHits(ovL),]$stat
  annot.all$pi <- out.gr[subjectHits(ovL),]$pi
  # qt <- quantile(annot.all$stat, quanti, na.rm = TRUE)
  qt <- quantile(test.var, quanti, na.rm = TRUE)
  gL <- annot.all$pi[annot.all$stat >= qt]
  g.n <- annot.all$name[annot.all$stat >= qt]
  gL <- gL[!duplicated(g.n)]
  names(gL) <- g.n[!duplicated(g.n)]
  return(gL)
}


par(mfrow=c(2,2))
annot.gr[annot.gr$name%in%"Skint1", ]
center <- 77238854
winsize <- 1e6
toplo <- res[res$sc == "ScOZjSD_3553" & res$st > center - winsize & res$en < center + winsize, ]

matplot(toplo[, grep("pi.", colnames(toplo))], pch = 16, type = "l")
matplot(toplo[, grep("rndsp.", colnames(toplo))], pch = 16, type = "l")


annot.gr[annot.gr$name%in%"H2-Eb1", ]
center <- 17600000
winsize <- 1e6
toplo <- res[res$sc == "ScOZjSD_2657"  & res$st > center - winsize & res$en < center + winsize, ]

matplot(toplo[, grep("pi", colnames(toplo))], pch = 16, type = "l")
matplot(toplo[, grep("rndsp", colnames(toplo))], pch = 16, type = "l")



#### ORs

par(mfrow=c(2,2))
annot.gr[annot.gr$name%in%"Olfr1019", ]
sca <- names(sort(table(res$sc), decreasing = TRUE)[3])

center <- 30533341
winsize <- 5e5
toplo <- res[res$sc == "ScOZjSD_3553" & res$st > center - winsize & res$en < center + winsize, ]

layout(1)
par(bty = "l", cex = 0.8, pch = 16)
matplot(
  toplo$st/1e6,
  toplo[, grep("pi.", colnames(toplo))],
  pch = 16, 
  type = "o", 
  lty = 1,
  lwd = 1.5,
  main = "",
  ylab = "Value",
  xlab = "Position (Mb)"
)
title("a. Diversity OR cluster", adj = 0)
matplot(
  toplo$st/1e6,
  toplo[, grep("rndsp.", colnames(toplo))],
  pch = 16, 
  type = "o", 
  lty = 1, 
  lwd = 1.5,
  main = "",
  ylab = "Value",
  xlab = "Position (Mb)"
)
title("b. Divergence OR cluster", adj = 0)


