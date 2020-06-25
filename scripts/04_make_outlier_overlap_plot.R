# plot overlap between gene lists

source("./scripts/Resequencing-data-analysis-functions.R")
require(GenomicRanges)

load("./data/annotation.rda")
load("./data/pi_rndsp.rda")

samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")


pis <- stats[, grep("pi", colnames(stats))]
rndsps <- stats[, grep("rndsp", colnames(stats))]

tab <- list()
for(i in 1:12) {
  x <- rndsps[[i]]
  pi.out <- getGeneDiv(rndsps[[i]], stats, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99, 1)
  qt.pi <- colMeans(sapply(pi.out, ">", pis[,i]))
  tab[[i]] <- qt.pi[!is.na(qt.pi)]
}
uniq.name <- unique(unlist(lapply(tab, function(x) names(x))))
df <- matrix(NA, ncol=length(tab), nrow=length(uniq.name))
rownames(df) <- uniq.name
for(i in 1:12) df[names(tab[[i]]), i] <- tab[[i]]

p.ord <- order(rowSums(!is.na(df)), decreasing = TRUE)
ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)
df <- df[p.ord, ][, ind.ord]

getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "PiYG"))

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%m%d_%H%M%S")

pdf(file = paste0(output.dir, stamp, "_Outlier_diversity.pdf"), width = 9.2, height = 7)

layout(matrix(c(1,2)))

par(mar=c(4,6,1,2), oma=c(1,1,1,1))
plot(NA, xlim = c(1, nrow(df)), ylim = c(1, ncol(df) + 2), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1:nrow(df), 1, 1:nrow(df), 12, col = rgb(.9,.9,.9))
segments(1, 1:12, nrow(df), 1:12, col = rgb(.95,.95,.95))
for(i in 1:12) {
  y <- df[, i]
  points(
    which(!is.na(y)), rep(i, sum(!is.na(y))),
    pch = 15, cex = 1.5,
    col = getPalette(100)[as.numeric(cut(y[!is.na(y)],breaks = seq(0,1,length.out = 100)))]
  )
}
axis(2, at = 1:ncol(df),
     labels = samp.names[-9][ind.ord],
     font = 3,
     lwd = 0, las = 2, cex.axis = 0.7, mgp = c(3, 0, 0))
axis(1, at = 1:nrow(df),
     labels = rownames(df),
     lwd = 0, las = 2, cex.axis = 0.6, mgp = c(3, 0, 0))

title("a. Outlier olfactory and vomeronasal receptors", adj = 0.4, cex.main = 0.9, outer = F)

tab <- list()
for(i in 1:12) {
  x <- rndsps[[i]]
  pi.out <- getGeneDiv(rndsps[[i]], stats, annot.gr[grep("H2|Skint|Clec|Cd200|Ig|Cfh|Cdk|Gbp", annot.gr$name),], quanti = 0.99, 1)
  qt.pi <- colMeans(sapply(pi.out, ">", pis[,i]))
  tab[[i]] <- qt.pi[!is.na(qt.pi)]
}
uniq.name <- unique(unlist(lapply(tab, function(x) names(x))))
df <- matrix(NA, ncol=length(tab), nrow=length(uniq.name))
rownames(df) <- uniq.name
for(i in 1:12) df[names(tab[[i]]), i] <- tab[[i]]

p.ord <- order(rowSums(!is.na(df)), decreasing = TRUE)
ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)
df <- df[p.ord, ][, ind.ord]



# par(mar=c(4,6,1,2), oma=c(1,1,1,1))
plot(NA, xlim = c(1, nrow(df)), ylim = c(1, ncol(df) + 2), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1:nrow(df), 1, 1:nrow(df), 12, col = rgb(.9,.9,.9))
segments(1, 1:12, nrow(df), 1:12, col = rgb(.95,.95,.95))
for(i in 1:12) {
  y <- df[, i]
  points(
    which(!is.na(y)), rep(i, sum(!is.na(y))),
    pch = 15, cex = 1.5,
    col = getPalette(100)[as.numeric(cut(y[!is.na(y)],breaks = seq(0,1,length.out = 100)))]
  )
}
axis(2, at = 1:ncol(df),
     labels = samp.names[-9][ind.ord],
     font = 3,
     lwd = 0, las = 2, cex.axis = 0.7, mgp = c(3, 0, 0))
axis(1, at = 1:nrow(df),
     labels = rownames(df),
     lwd = 0, las = 2, cex.axis = 0.6, mgp = c(3, 0, 0))

title("b. Outlier immunity genes", adj = 0.4, cex.main = 0.9, outer = F)

par(xpd=T)
text(nrow(df)-2, 15, expression(paste(pi, " quantile")), font = 2, cex = 0.8)
points(seq(nrow(df)-4, nrow(df), length.out = 100), rep(13.8, 100), col =  getPalette(100), pch = 15, cex = 2.5)
text(c(nrow(df)-5, nrow(df)+1), rep(13.8, 2), c(0,1), cex = 0.7)

dev.off()



### Supp figure

# par(mfrow=c(4, 3), mar=c(4,4,2,2))
# # for (i in 1:12) {
# #   hist(dfall[,i], xlim = c(0, 1), xlab = "Q(pi)", main = samp.names[-9][ind.ord][i], breaks = 100)
# # }
# # for (i in 1:12) {
# #   hist(dfolf[,i], xlim = c(0, 1), xlab = "Q(pi)", main = samp.names[-9][ind.ord][i], breaks = 100)
# # }
# for (i in 1:12) {
#   qqplot(dfolf[,i], dfall[,i], main = samp.names[-9][ind.ord][i], xlim = c(0,1), ylim = c(0,1),
#          xlab = "OR quantiles", ylab = "Non OR quantiles", col = 1, pch = 16, cex = 0.5, asp=1)
#   par(new=TRUE)
#   qqplot(dfolf[,i], ppoints(500), main = "", xlim = c(0,1), ylim = c(0,1),
#           xlab = "", ylab = "", add = TRUE, col = 2, pch = 16, cex = 0.5, asp=1)
#   abline(0,1)
# }
# 
