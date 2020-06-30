### Computing popgen stats from VCF file

library(PopGenome)
vcf.fn <- "C:/gouy/data/voles/13V-22sc-Di-Nomiss-GT.vcf.gz"

getRNDsp <- function(data) {
  
  dxys <- data[, grep("dxy", colnames(data))]
  av.dxys <- colMeans(dxys, na.rm = TRUE)
  myod <- dxys[, grep("pop9", colnames(dxys))]
  
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
  
  return(rndsps)
}



vcf_handle <- .Call("VCF_open", vcf.fn)
sample.names <- .Call("VCF_getSampleNames", vcf_handle)

# marv <- c(grep("Mar", sample.names), grep("Eugene", sample.names))
# pops <- list(
#   Marv = sample.names[marv],
#   Mothers = sample.names[-marv]
# )

pops <- list(
  c("MarFTh497_wh"),
  c("MarCHBo17"), 
  c("MarCZD02"), 
  c("Eugene"),
  c("MagCHRi22_wh"),
  c("MduP01"),
  c("MoePBi05"),
  c("MbrRClm01_wh"),
  c("MglCHBk23"),
  c("MpeMi11"),
  c("Mcab10_wh"),
  c("MluP01"),
  c("MroFl38")
)

sclength <- read.table("./data/vol_scaff_len.txt",
                       sep = ",", header = FALSE)
sclength <- sclength[order(sclength[,2], decreasing=TRUE),][1:22,]
rownames(sclength) <- sclength[, 1]
sc.list <- as.character(sclength[, 1])
sclength[sc.list,]

res <- list()

for(scaff in as.character(sc.list)) {
  message(scaff)
    
  GENOME <- readVCF(
    vcf.fn,
    numcols = 1e5, # SNP-chunk size
    frompos = 1,
    topos = sclength[scaff, 2],
    tid = scaff, # chromosome id
    approx = FALSE, # if true, 0/1 and 1/0 go to 1
    samplenames = sample.names
  )
  
  GENOME <- set.populations(GENOME, pops, diploid = TRUE)
  win.size <- 5e4
  step <- 5e4
  slide <- sliding.window.transform(GENOME, win.size, step, type = 2)
  
  slide <- F_ST.stats(slide, mode = "nucleotide", FAST = TRUE) # includes div stats within and between
  
  if(is.null(names(pops))) names(pops) <- seq_along((pops))
  
  div.within <- slide@nuc.diversity.within / win.size
  colnames(div.within) <- paste0("pi.", names(pops))
  
  div.between <- data.frame(t(slide@nuc.diversity.between / win.size))
  colnames(div.between) <- paste0("dxy.", colnames(div.between))
  
  # fst <- data.frame(t(slide@nuc.F_ST.pairwise))
  # colnames(fst) <- paste0("fst.", colnames(fst))
  
  pos <- do.call("rbind", strsplit(gsub(pattern = " :", "", slide@region.names), split = " - "))
  colnames(pos) <- c("st", "en")
  
  res[[scaff]] <- data.frame(sc = scaff, pos, div.within, div.between)
}

# data frame preprocessing
res <- do.call("rbind", res)
rownames(res) <- NULL
res$st <- as.numeric(as.character(res$st)) - 1 # saved as factor in the original file
res$en <- as.numeric(as.character(res$en)) - 1
res$pos.cumul <- cumsum(rep(win.size, nrow(res)))/1e6
res$col.scaff <- as.numeric(as.factor(res$sc)) %% 2

# alternance of scaffold colors for plots
res$col.scaff[res$col.scaff==0] <- "grey"
res$col.scaff[res$col.scaff==1] <- "black"

RNDsp <- getRNDsp(res)

# branch specific value
res <- cbind(res, RNDsp)

save(res, file="/home/gouy/vole-genome/13V-22sc-Di-Nomiss-GT-v2.rda")

