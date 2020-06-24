### Script to perform GO enrichment tests
# and visualize results

run.GO <- FALSE
getStringDB <- FALSE

library(GenomicRanges)
library(STRINGdb)
source("./scripts/Resequencing-data-analysis-functions.R")

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
  dxys <- res[, grep("dxy", colnames(res))]
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
    GL <- unique(getGeneList(unlist(unname(c(YY))), res, annot.gr, quant = 0.99))
    # GL <- unique(getGeneList2(unlist(unname(c(YY))), res, annot.gr))
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

fdr.threshold <- 0.1
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

for(i in 1:12) {
  df[tab[[i]]$description, i] <- paste(tab[[i]]$number_of_genes, signif(tab[[i]]$p_value, 1), sep = "; ")
  pvs[tab[[i]]$description, i] <- tab[[i]]$fdr
  hts[tab[[i]]$description, i] <- tab[[i]]$number_of_genes
  g.names[tab[[i]]$description, i] <- tab[[i]]$preferredNames
}

g.names[1,]

# test <- apply(g.names[1:2,], 1, function(x) {
#   gn <- do.call(c, sapply(x, strsplit, ","))
#   shared2 <- tableduplicated(gn)]
#   shared.all <- 
#   tot <- length(unique(gn))
# })
# 
# i <- 1
# tab[[1]]
# df
### TO GET: 
# total genes
# total unique significant genes
# total shared significant genes



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

# visualise results
par(mfrow=c(1,1), mar=c(7,20,0,0), oma=c(0,0,1,0), mgp=c(3, 0, 0))
plot(NA, xlim = c(1, ncol(df3) + 2), ylim = c(1, nrow(df3)), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1, 1:nrow(df3), 12, 1:nrow(df3), col = rgb(.9,.9,.9))
segments(1:12, 1, 1:12, nrow(df3), col = rgb(.95,.95,.95))

couleurs <- rgb(0, 0.192, 0.325, col.sc[,i][!is.na(col.sc[,i])]) # prussian blue
for(i in 1:12) {
  points(
    rep(i, sum(!is.na(t(df3[,i])))), t(df3[,i])[!is.na(t(df3[,i]))], 
    pch = 16, 
    cex = cx.sc[,i][!is.na(cx.sc[,i])],
    col = couleurs
  )
}
axis(1, at = 1:ncol(df3), labels = samp.names[-9][ind.ord], lwd = 0, las = 2, cex.axis = 0.8)
axis(2, at = 1:nrow(df3), labels = rownames(df2)[condi][p.ord], lwd = 0, las = 2, cex.axis = 0.7)

Y.le <- round(0.9*max(df3, na.rm = TRUE))
text(c(13.5), 1.1*Y.le, "Legend", font = 2, cex = 0.8)
points(c(13,14), rep(Y.le, 2), cex = range(cx.sc[!is.na(cx.sc)]))
text(c(13,14), rep(Y.le, 2), range(hts, na.rm = TRUE), pos = 1, cex = 0.6)
text(c(13.5), Y.le, "# genes", pos = 3, cex = 0.8)
sc.fac <- 0.9
points(c(13,14), rep(sc.fac*Y.le, 2), col = rgb(0, 0.192, 0.325, range(col.sc[!is.na(col.sc)])), pch = 16)
text(c(13,14), rep(sc.fac*Y.le, 2), round(range(-log10(pvs), na.rm = TRUE), 1), pos = 1, cex = 0.6)
text(c(13.5), sc.fac*Y.le, expression(-log[10](italic(q))), pos = 3, cex = 0.8)

title("Gene Ontology terms enriched in divergent loci", adj = 0.5, cex.main = 0.9, outer = TRUE)


### Write results as a CSV table

# to.write <- as.data.frame(df)
# colnames(to.write) <- samp.names[-9]
# to.write$GOid <- go.id
# to.write$GOsz <- go.size
# write.table(to.write, file = "GO-results.tsv", quote = FALSE, sep = "\t")
