### Get RNDsp and annotation

# Get PopGenome output
load("./data/13V-22sc-Di-Nomiss-GT.rda")
source("./scripts/Resequencing-data-analysis-functions.R")

# Prepare annotation
require(GenomicRanges)

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
save(annot.gr, file="./data/annotation.rda")

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
names(dspec) <- paste0("rndsp.", names(dspec))

d.av.myo <- rowMeans(myod)
rndsps <- lapply(dspec, function(x) {x / d.av.myo})
cond <- d.av.myo <= quantile(d.av.myo, 0.025)

rndsps <- lapply(rndsps, function(x){
  x[cond] <- NA
  return(x)
})

pis <- res[, grep("pi", colnames(res))]
stats <- cbind(res[,1:3], pis, rndsps)

save(stats, file = "./data/pi_rndsp.rda")
