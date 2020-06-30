# MAke figure 2
load("./data/13V-22sc-Di-Nomiss-GT.rda")
source("./scripts/00_utils.R")

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%m%d_%H%M%S")

LT <- c(1:3, rep(1, 10))
# samp.names <- c("M. arv. F", "M. arv. CH", "M. arv. CZ", "Eugene",
#                 "M. agr.", "M. duo.", "M. oec.", "M. bra.", "M. gla.",
#                 "M. pen.", "M. cab.", "M. lus.", "M. ros.")
samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")

span.loess <- 0.01

vole.col <- c(rep("#228B22", 4), RColorBrewer::brewer.pal(9, "Set1"))
vole.col[c(2,3,1,4,11)] <- RColorBrewer::brewer.pal(5, "Paired")


divs <- res[, grep("pi", colnames(res))]

pdf(file = paste0(output.dir, stamp, "_History_v2.pdf"), width = 10, height = 7.5)

layout(matrix(c(1,2,3,3), ncol=2))
par(xpd=FALSE)
par(mar = c(4,10,4,2))
o <- order(colMeans(divs), decreasing = TRUE)
barplot(colMeans(divs)[o], names.arg = samp.names[o],
        col = vole.col[o], horiz = TRUE, las = 1, cex.axis = 0.8,
        xlab = expression(pi), border = F)
title("a. Average nucleotide diversity", adj = 0)

### B. PSMC

l.f <- list.files("./data/psmc/", full.names = TRUE)
files <- l.f[grep("Mar|Mcab", l.f)]

pp <- lapply(files, psmc.result, s=100, g=0.5, mu=0.8e-8)
pp <- lapply(pp, function(x) {x[,2] <- x[,2]/1e4; return(x)})
couls <- RColorBrewer::brewer.pal(5, "Paired")
LTI <- c(1,1,1,1,1)
plot(pp[[1]], type = "s", main = "", log = c("x"), lwd = 2,lty=LTI[1],col=couls[1],
     ylim = c(0, 65), bty = "l", xlim = c(1e3, 1.5e6),
     xlab = "Years ago", ylab = "Effective population size (x 10e4)")
title("b. PSMC-inferred demographic changes", adj = 0)

for(i in 2:5) {
  lines(pp[[i]], col = couls[i],lty=LTI[i], type = "s", lwd = 2)
}
legend("topright", lwd=2, lty=LTI,
       col=couls,
       c("M. arvalis I", "M. arvalis E", "M. arvalis W", "M. arvalis C", "M. cabrerae"))

### C. PHYLO TREE

library(ape)
library(phytools)
raxml_file <- "C:/gouy/data/voles/make_phylogeny/RAxML_bipartitionsBranchLabels.20200619-vole-ml-gtrgamma-bootstrap-AVX2"
raxml <- treeio::read.raxml(raxml_file)

tree <- raxml@phylo
tree$tip.label
tree$tip.label <- c("M. arvalis I", "M. arvalis E", "M. arvalis C", "M. levis", "M. duodecimcostatus",
                    "M. lusitanicus", "M. cabrerae", "M. agrestis", "M. pennsylvanicus",
                    "M. glareolus", "M. oeconomus", "M. brandtii", "M. arvalis W")
tree <- reroot(tree, node.number = 10, position = max(tree$edge.length) / 3.2)

tree$node.labels <- raxml@data$bootstrap

par(mar = c(4,1,4,1))
plot(tree, label.offset = .005, type="phylogram",
     show.node.label = F,
     align.tip.label = F, font = 2, edge.width = 2, edge.col = "slategrey")

tiplabels(pch = 16, col = "slategrey")
nodelabels(text = raxml@data$bootstrap, cex = .6, bg = "white",
           adj = c(1.2, 1.2), frame = "n", font = 2)
axisPhylo(backward = FALSE, lwd = 2, col = "slategrey")
title("c. Maximum likelihood phylogenetic tree", adj = 0)
title(xlab = "Nucleotide substitutions per site", font.lab=2)

dev.off()

# 
# 
par(mfrow=c(2,2))
plot(tree, label.offset = .005,type="fan",
     show.node.label = F,
     align.tip.label = F, font=2, edge.width=2, edge.col="slategrey")
plot(tree, label.offset = .005,type="radial",
     show.node.label = F,
     align.tip.label = F, font=2, edge.width=2, edge.col="slategrey")
plot(tree, label.offset = .005,type="phylogram",
     show.node.label = F,
     align.tip.label = F, font=2, edge.width=2, edge.col="slategrey")
plot(tree, label.offset = .005,type="unrooted",
     show.node.label = F,
     align.tip.label = F, font=2, edge.width=2, edge.col="slategrey")

# 
# 
# ### PSMC plot
# 
# ### SUPPLEMENTARY: DIVERSITY DISTRIBUTION
# 
# # par(mfrow=c(1, 2), mar = c(4,6,4,4))
# plot(density(log10(res$pi.1+1e-5)), col = vole.col[1], lwd = 1, ylim = c(0,2.2),
#      main = "", lty = LT[1], xlim = c(-5.5,-1.5),
#      xlab = expression(log[10](pi + 10^-5)), bty = "L")
# title("b. Nucleotide diversity distributions", adj = 0)
# 
# for(i in 2:13) {
#   to.plot <- eval(parse(text=paste0("res$pi.", i)))
#   lines(density(log10(to.plot + 1e-5)), col = vole.col[i], lwd = 1, lty = LT[i])
# }
# par(xpd=TRUE)
# legend(-2, 2, samp.names, col = vole.col, lwd = 2, cex = 0.8, lty = LT)
# par(xpd=FALSE)