### Script to perform GO enrichment tests
# and visualize results

run.GO <- FALSE
getStringDB <- FALSE

library(GenomicRanges)
library(STRINGdb)
source("./scripts/00_utils.R")

output.dir <- "./plots"
if(!dir.exists(output.dir)) dir.create(output.dir)
stamp <- format(Sys.time(), "/%Y%m%d_%H%M%S")

samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")

load("./data/annotation.rda")
load("./data/pi_rndsp.rda")
rndsps <- stats[, grep("rndsp", colnames(stats))]

if(run.GO) {
  ### Load annotation
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
  
  ### Load divergence measures
  load("./data/13V-22sc-Di-Nomiss-GT.rda")
  dspec <- list()
  dxys <- stats[, grep("dxy", colnames(stats))]
  av.dxys <- colMeans(dxys, na.rm = TRUE)
  
  myod <- dxys[, grep("pop9", colnames(dxys))]
  
  samp.names <- c("M. arvalis W","M. arvalis I","M. arvalis E","M. arvalis C",
                  "M. agrestis","M. duodecimcostatus","M. oeconomus","M. brandti","M. glareolus",
                  "M. pennsylvanicus","M. cabrerae","M. lusitanicus","M. levis")
  
  ### Compute RNDsp
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
  
  # prepare annotation
  if(getStringDB) {
    string_db <- STRINGdb$new(version="11", species = 10090, score_threshold = 400, input_directory = "")
    bg <- string_db$map(data.frame(gene=annot.gr$name), "gene")
    string_db$set_background(bg$STRING_id)
    save(string_db, file = "./data/string_db.rda")
  } else {
    load(string_db, file = "./data/string_db.rda")
  }
  
  # perform GO test
  tab <- list()
  for(i in 1:12) {
    print(i)
    eval(parse(text=paste0("YY <- rndsps[[", i, "]]")))
    YY[cond] <- NA
    GL <- unique(getGeneList(unlist(unname(c(YY))), stats, annot.gr, quant = 0.99))
    # GL <- unique(getGeneList2(unlist(unname(c(YY))), stats, annot.gr))
    ddf <- data.frame(gene = GL)
    mapped <- string_db$map(ddf, "gene")
    hits <- mapped$STRING_id
    enrichmentGO <- string_db$get_enrichment(hits, category = "Process", methodMT = "fdr", iea = TRUE)
    tab[[i]] <- enrichmentGO[enrichmentGO$fdr < 0.2, ]
  }
  names(tab) <- samp.names[-9]
  
  save(tab, file = "./data/GO-table.rda")
  
} else {
  
  load(file = "./data/GO-table.rda")
  
}

fdr.threshold <- 0.05
tab <- lapply(tab, function(x) x[x$fdr < fdr.threshold, ])

uniq.name <- unique(unlist(lapply(tab, function(x) x[, "description"])))
go.id <- (unlist(lapply(tab, function(x) x[, "term"])))
go.size <- (unlist(lapply(tab, function(x) x[, "number_of_genes_in_background"])))

go.size <- go.size[!duplicated(go.id)]
go.id <- go.id[!duplicated(go.id)]

df <- matrix(NA,ncol=length(tab),nrow=length(uniq.name))
rownames(df) <- uniq.name
pvs <- matrix(NA,ncol=length(tab),nrow=length(uniq.name))
rownames(pvs) <- uniq.name
hts <- matrix(NA,ncol=length(tab),nrow=length(uniq.name))
rownames(hts) <- uniq.name
g.names <- matrix(NA,ncol=length(tab),nrow=length(uniq.name))
rownames(g.names) <- uniq.name

sizes <- matrix(NA,ncol=length(tab),nrow=length(uniq.name))
rownames(sizes) <- uniq.name

for(i in 1:12) {
  df[tab[[i]]$description, i] <- paste(tab[[i]]$number_of_genes, signif(tab[[i]]$p_value, 1), sep = "; ")
  pvs[tab[[i]]$description, i] <- tab[[i]]$fdr
  hts[tab[[i]]$description, i] <- tab[[i]]$number_of_genes
  g.names[tab[[i]]$description, i] <- tab[[i]]$preferredNames
}


shared.genes <- apply(g.names, 1, function(x) {
  gn <- unname(do.call(c, sapply(x, strsplit, ",")))
  shared <- table(gn)
  return(c(shared))
})

tot.g <- tab[[i]]$number_of_genes_in_background
shar.g <- lengths(shared.genes[tab[[i]]$description])
uniq.g <- lapply(shared.genes[tab[[i]]$description], function(x) sum(x==1))
GO.labels <- paste0("(", uniq.g, "|", shar.g, "|", tot.g, ")")


shared.genes <- unlist(unname(shared.genes))

shared.genes <- shared.genes[!duplicated(names(shared.genes))]

uniq.name <- names(shared.genes)
pis <- stats[, grep("pi", colnames(stats))]
tab <- list()
for(i in 1:12) {
  x <- rndsps[[i]]
  # getGeneList(x, stats, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99)
  # getGeneList(x, stats, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99)
  pi.out <- getGeneDiv(rndsps[[i]], stats, annot.gr[annot.gr$name %in% uniq.name,], quanti = 0.99, 1)
  
  # pi.out <- getGeneDiv(rndsps[[i]], stats, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99, 1)
  qt.pi <- colMeans(sapply(pi.out, ">", pis[,i]))
  tab[[i]] <- qt.pi[!is.na(qt.pi)]
}
div.quantiles <- matrix(NA, ncol=length(tab), nrow=length(uniq.name))
rownames(div.quantiles) <- uniq.name
for(i in 1:12) div.quantiles[names(tab[[i]]), i] <- tab[[i]]

# 


library(Hmisc)
rownames(df) <- capitalize(paste(rownames(df), GO.labels))

df2 <- df[unname(go.size >= 5 & go.size < 1500), ]
df2[!is.na(df2)] <- 1
df2[is.na(df2)] <- 0
df3 <- apply(df2, 2, as.numeric)

condi <- colSums(apply(df3, 1, "==", 0)) < 11
df3 <- df3[condi, ]
p.ord <- hclust(dist(df3))$order
df3 <- df3[p.ord, ]
ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)

pvs <- pvs[condi, ][p.ord, ][, ind.ord]
hts <- hts[condi, ][p.ord, ][, ind.ord]


df3 <- apply(df3, 2, "*", 1:nrow(df3))


df3[df3 == 0] <- NA
df3 <- df3[, ind.ord]

col.sc <- (-log10(pvs) + 20) / max(-log10(pvs) + 20, na.rm = TRUE)
cx.sc <- 2 * log10(hts + 1) / max(log10(hts + 1), na.rm = TRUE)

### Make the plot
# pdf(file = paste0(output.dir, stamp, "_GO_Enrichment_v2.pdf"), width = 9, height = 7)

# layout(matrix(c(1,1,1,1,2,3), nrow=2))

par(mar=c(5,30,1,4), oma=c(1,1,3,0), mgp=c(3, 1, 0), mfrow=c(1,1))
plot(NA, xlim = c(1, ncol(df3) + 2), ylim = c(1, nrow(df3)), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1, 1:nrow(df3), 12, 1:nrow(df3), col = rgb(.9,.9,.9))
segments(1:12, 1, 1:12, nrow(df3), col = rgb(.95,.95,.95))


getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "PiYG"))
mapCol <- function(y, n = 100) {
  getPalette(n)[as.numeric(cut(y[!is.na(y)], breaks = seq(0, 1, length.out = n)))]
}
# couleurs <- rgb(0, 0.192, 0.325, col.sc[,i][!is.na(col.sc[,i])]) # prussian blue



for(i in 1:12) {
  y <- t(df3[,i])[!is.na(t(df3[,i]))]
  points(
    rep(i, sum(!is.na(t(df3[,i])))), y, 
    pch = 16, 
    cex = cx.sc[,i][!is.na(cx.sc[,i])],
    col = mapCol(y)
  )
}
axis(1, at = 1:ncol(df3), labels = FALSE, lwd = 0, las = 2, cex.axis = 0.9)
text(par("usr")[1]+1:ncol(df3), 0, labels = samp.names[-9][ind.ord], srt = 45, pos = 2, xpd = TRUE, cex = 0.9)

axis(2, at = 1:nrow(df3), labels = rownames(df2)[condi][p.ord], lwd = 0, las = 2, cex.axis = 0.9)

par(xpd=TRUE)
Y.le <- 22
text(-27, 5, "Legend", font = 2, cex = 0.8)
aa <- range(cx.sc[!is.na(cx.sc)])
points(seq(-28,-26, length.out = 20), rep(3.2, 20), cex = seq(aa[1],aa[2], length.out = 20), pch = 16)
text(c(-29,-25), rep(3.2, 2), range(hts, na.rm = TRUE), cex = 0.7)
text(-27, 3.5, "Number of genes", pos = 3, cex = 0.8)
sc.fac <- 0.95
bb <- range(col.sc[!is.na(col.sc)])

points(seq(-28,-26, length.out = 20), rep(2.5, 20),
       col = rgb(0, 0.192, 0.325, seq(bb[1],bb[2], length.out = 20)),
       pch = 15)


(range((pvs), na.rm = TRUE))
text(c(13,14), rep(sc.fac*Y.le, 2), c("10^-31", 0.05), pos = 1, cex = 0.6)
text(13.5, sc.fac*Y.le, expression(-log[10](italic(q))), pos = 3, cex = 0.8)

title("Gene Ontology terms enriched in divergent loci", adj = 0.5, cex.main = 1, outer = TRUE)


# dev.off()
# pdf(file = paste0(output.dir, stamp, "_GO_Enrichment_Supp1.pdf"), width = 5, height = 7)

par(mfrow=c(2,1), mar=c(7,4,2,4))
hist(shared.genes, main = "", cex.axis = 0.9,
     xlab = "Number of species",  ylab = "Count",
     border = FALSE)

title("a. Genes occurence in different lineages", adj = 0., cex.main = 0.9)

hist(div.quantiles, main = "", cex.axis = 0.9, 
     xlab = "Nucleotide diversity quantile", ylab = "Count",
     border = FALSE)
title("b. Diversity distribution of divergent genes", adj = 0., cex.main = 0.9)


# dev.off()

### Write results as a CSV table

# to.write <- as.data.frame(df)
# colnames(to.write) <- samp.names[-9]
# to.write$GOid <- go.id
# to.write$GOsz <- go.size
# write.table(to.write, file = "GO-results.tsv", quote = FALSE, sep = "\t")
