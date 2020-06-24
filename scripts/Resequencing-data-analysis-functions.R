### Functions for the resequencing data analysis

# function for genome wide plot with smoothing
plotRes <- function(data.fr, var.y, titre, y.lab, scaff = 0, span.loess = 0.01, ylim = NULL, xlim = NULL) {
  quant <- quantile(var.y, c(0.99), na.rm = TRUE)
  if(scaff != 0) {
    cond <- data.fr$sc == scaff
    data.fr <- data.fr[cond, ]
    var.y <- var.y[cond]
  }
  
  if(is.null(ylim) & is.null(xlim)) plot(data.fr$pos.cumul, var.y, bty = "n", pch = 16, cex = 0.5,
       col = data.fr$col.scaff, main = "", xlab = "Position (Mb)", ylab = y.lab)
  else plot(data.fr$pos.cumul, var.y, bty = "n", pch = 16, cex = 0.5,
                                         ylim = ylim, col = data.fr$col.scaff, main = "", xlab = "Position (Mb)", ylab = y.lab)
  
  mod <- loess(var.y ~ data.fr$pos.cumul, span = span.loess, na.action = na.exclude)
  lines(data.fr$pos.cumul, predict(mod), col = "firebrick", lwd = 2)
  abline(h= quant, lty = 2, col = "firebrick", lwd = 1.5)
  title(titre, adj = 0)
}


getGeneList <- function(test.var, res, annot.gr, quanti = 0.99) {
  out.gr <- GRanges(
    seqnames = res$sc,
    ranges = IRanges(start = res$st, end = res$en),
    stat = test.var
  )
  ovL <- findOverlaps(annot.gr, out.gr)
  annot.all <- annot.gr[queryHits(ovL), ]
  
  annot.all$stat <- out.gr[subjectHits(ovL),]$stat
  # qt <- quantile(annot.all$stat, quanti, na.rm = TRUE)
  qt <- quantile(test.var, quanti, na.rm = TRUE)
  gL <- unique(annot.all$name[annot.all$stat >= qt])
  return(gL)
}
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


plot.dxy.tree <- function(dxy.vec, samp.names, title = "Nucleotide divergence based phylogeny") {
  dist.mat <- matrix(NA, nrow=length(samp.names), ncol=length(samp.names))
  
  ids <- strsplit(gsub("dxy.pop", "", names(dxy.vec)), ".pop")
  for(i in seq_along(dxy.vec)) {
    dist.mat[as.numeric(ids[[i]][1]), as.numeric(ids[[i]][2])] <- dxy.vec[i]
    dist.mat[as.numeric(ids[[i]][2]), as.numeric(ids[[i]][1])] <- dxy.vec[i]
  }
  
  diag(dist.mat) <- 0
  hc <- hclust(as.dist(dist.mat), method = "average") # "average" <-> UPGMA
  # tree
  plot(hc, labels = samp.names,
       xlab = "", sub = "",
       ylab = expression(d[xy]),
       main = title) 
  
}

get.dist.mat <- function(dxy.vec) {
  dist.mat <- matrix(NA, nrow=length(samp.names), ncol=length(samp.names))
  
  ids <- strsplit(gsub("dxy.pop", "", names(dxy.vec)), ".pop")
  for(i in seq_along(dxy.vec)) {
    dist.mat[as.numeric(ids[[i]][1]), as.numeric(ids[[i]][2])] <- dxy.vec[i]
    dist.mat[as.numeric(ids[[i]][2]), as.numeric(ids[[i]][1])] <- dxy.vec[i]
  }
  
  diag(dist.mat) <- 0
  return(dist.mat) 
}

getDa <- function(data) {
  colnames(data)
  
  col.ids.dxy <- grep(pattern = "dxy.pop", colnames(data))
  
  all.das <- sapply(colnames(data[, col.ids.dxy]), function(x) {
    
    nums <- as.numeric(unlist(strsplit(gsub("dxy.pop", "", x), ".pop")))
    
    col.ids.pi <- c(
      paste0("pi.", nums[1]),
      paste0("pi.", nums[2])
    )
    
    da <- data[, colnames(data) == x] - ((data[, colnames(data) %in% col.ids.pi[1]] - data[, colnames(data) %in% col.ids.pi[2]]) / 2)
    return(da)
  })
  
  return(all.das)
}


getDf <- function(data) {
  colnames(data)
  
  col.ids.dxy <- grep(pattern = "dxy.pop", colnames(data))
  
  all.das <- sapply(colnames(data[, col.ids.dxy]), function(x) {
    
    nums <- as.numeric(unlist(strsplit(gsub("dxy.pop", "", x), ".pop")))
    col.ids.pi <- c(
      paste0("pi.", nums[1]),
      paste0("pi.", nums[2])
    )
    
    da <- data[, colnames(data) == x] - rowSums(data[, colnames(data) %in% col.ids.pi])
    return(da)
  })
  
  return(all.das)
}


getRND <- function(data, outgroup = 9) {

  col.ids.dxy <- grep(pattern = "dxy.pop", colnames(data))
  col.ids.dxy.out <- c(
    grep(pattern = "dxy.pop9.pop", colnames(data)),
    grep(pattern = "dxy.pop[0-9]{1,2}.pop9$", colnames(data))
  )
  
  rnd <- sapply(colnames(data)[col.ids.dxy[!col.ids.dxy %in% col.ids.dxy.out]], function(i) {
    nums <- as.numeric(unlist(strsplit(gsub("dxy.pop", "", i), ".pop")))
    
    col.xy <- c(
      paste0("dxy.", paste0("pop", nums, collapse=".")),
      paste0("dxy.", paste0("pop", rev(nums), collapse="."))
    )
    col.out <- c(
      paste0("dxy.", paste0("pop", c(nums[1], 9), collapse=".")),
      paste0("dxy.", paste0("pop", rev(c(nums[1], 9)), collapse=".")),
      paste0("dxy.", paste0("pop", c(nums[2], 9), collapse=".")),
      paste0("dxy.", paste0("pop", rev(c(nums[2], 9)), collapse="."))    
    )
    
    dout <- rowMeans(data[, colnames(data) %in% col.out])
    dxy <- data[, colnames(data) %in% col.xy]
    rnd <- dxy / dout
    return(rnd)
  })
  colnames(rnd) <- gsub("dxy", "rnd", colnames(rnd))
  return(rnd)
}



cmyk <- function(C,M,Y,K) {
  
  C <- C / 100.0
  M <- M / 100.0
  Y <- Y / 100.0
  K <- K / 100.0
  
  n.c <- (C * (1-K) + K)
  n.m <- (M * (1-K) + K)  
  n.y <- (Y * (1-K) + K)
  
  r.col <- ceiling(255 * (1-n.c))
  g.col <- ceiling(255 * (1-n.m))
  b.col <- ceiling(255 * (1-n.y))
  
  x <- col2rgb(sprintf("#%02s%02s%02s",
                       as.hexmode(r.col), 
                       as.hexmode(g.col), 
                       as.hexmode(b.col)))
  
  return(rgb(t(x)/255))
}



plotOut <- function(y, ylim = c(0,1), scale = 0.7, clust = 200, quant = 0.99,
                    g.nam = TRUE, cx = 0.8, xlab = "", ylab = "") {
  y <- unlist(unname(c(y)))
  out.gr <- GRanges(
    seqnames = res$sc,
    ranges = IRanges(start = res$st, end = res$en),
    stat = y,
    pos.cumul = res$pos.cumul
  )
  ovL <- findOverlaps(annot.gr, out.gr)
  annot.all <- annot.gr[queryHits(ovL), ]
  
  annot.all$stat <- out.gr[subjectHits(ovL),]$stat
  annot.all$pos.cumul <- out.gr[subjectHits(ovL),]$pos.cumul
  Q <- quantile(annot.all$stat, quant, na.rm = TRUE)
  

  bool <- y>Q
  bool[is.na(bool)] <- FALSE

  bins <- cut(1:nrow(res), nrow(res)/3)
  condi <- tapply(bool, bins, sum) >= 2
  bool.all <- as.vector(condi[bins])

  ovL <- findOverlaps(annot.all, out.gr[bool.all,])
  annot.all <- annot.all[queryHits(ovL), ]
  
  
  coul <- res$col.scaff
  coul[bool & bool.all] <- "red"
  coul[bool & !bool.all] <- rgb(0.9,0.9,0.9,0.8)
  
  CXX <- rep(0.5, length(res$col.scaff))
  CXX[bool & bool.all] <- 0.75
  
  annot.all <- annot.all[!is.nan(annot.all$stat) & !is.na(annot.all$stat),]
  g.pos <- annot.all[annot.all$stat>Q,]$pos.cumul
  g.id <- annot.all[annot.all$stat>Q,]$name
  
  # getGeneList(y, res, annot.gr)
  
  cond <- !duplicated(g.id) & !is.na(g.id)
  
  g.id <- g.id[cond]
  g.pos <- g.pos[cond]
  g.id <- g.id[order(g.pos)]
  g.pos <- g.pos[order(g.pos)]
  hc <- hclust(dist(g.pos))
  y.pos <- unlist(tapply(g.id, cutree(hc, h = clust), function(x) 1:length(x)))
  y.pos <- scales::rescale(y.pos, rev(c(scale*ylim[2],ylim[2])))
  
  plot(res$pos.cumul, y, bty = "n", xlab = xlab, ylab = ylab,
       type = "p", pch = 16, cex = CXX, col = coul, ylim = ylim, xaxt = "n")
  segments(g.pos,ylim[1],g.pos,ylim[2], lty = 2, col = "pink")
  if(g.nam) text(g.pos, y.pos, g.id, cex = cx)
  
  mod <- loess(y ~ res$pos.cumul, span = 0.01)
  lines(res$pos.cumul[!is.na(y)], predict(mod), lwd = 2, col = "firebrick")
}
plotOut2 <- function(y, ylim = c(0,1), scale = 0.7, clust = 200, quant = 0.99,
                    g.nam = TRUE, cx = 0.8, xlab = "", ylab = "") {
  y <- unlist(unname(c(y)))

  coul <- res$col.scaff

  CXX <- rep(0.5, length(res$col.scaff))
 
  plot(res$pos.cumul, y, bty = "n", xlab = xlab, ylab = ylab,
       type = "p", pch = 16, cex = CXX, col = coul, ylim = ylim, xaxt = "n")
  mod <- loess(y ~ res$pos.cumul, span = 0.01)
  lines(res$pos.cumul[!is.na(y)], predict(mod), lwd = 2, col = "firebrick")
}


getGeneList2 <- function(y, res, annot.gr, quant = 0.99) {
  y <- unlist(unname(c(y)))
  out.gr <- GRanges(
    seqnames = res$sc,
    ranges = IRanges(start = res$st, end = res$en),
    stat = y,
    pos.cumul = res$pos.cumul
  )
  ovL <- findOverlaps(annot.gr, out.gr)
  annot.all <- annot.gr[queryHits(ovL), ]
  
  annot.all$stat <- out.gr[subjectHits(ovL),]$stat
  annot.all$pos.cumul <- out.gr[subjectHits(ovL),]$pos.cumul
  Q <- quantile(annot.all$stat, quant, na.rm = TRUE)
  
  
  bool <- y>Q
  bool[is.na(bool)] <- FALSE
  # bool.seq <- ave(bool, cumsum(!bool), FUN = cumsum)
  
  bins <- cut(1:nrow(res), nrow(res)/3)
  condi <- tapply(bool, bins, sum) >= 2
  bool.all <- as.vector(condi[bins])
  
  ovL <- findOverlaps(annot.all, out.gr[bool.all,])
  annot.all <- annot.all[queryHits(ovL), ]

  annot.all <- annot.all[!is.nan(annot.all$stat) & !is.na(annot.all$stat),]
  g.pos <- annot.all[annot.all$stat>Q,]$pos.cumul
  g.id <- annot.all[annot.all$stat>Q,]$name
  
  cond <- !duplicated(g.id) & !is.na(g.id)
  g.id <- g.id[cond]
  return(g.id)

}

getGeneScores <- function(y) {
  y <- unlist(unname(c(y)))
  out.gr <- GRanges(
    seqnames = res$sc,
    ranges = IRanges(start = res$st, end = res$en),
    stat = y,
    pos.cumul = res$pos.cumul
  )
  ovL <- findOverlaps(annot.gr, out.gr)
  annot.all <- annot.gr[queryHits(ovL), ]
  
  annot.all$stat <- out.gr[subjectHits(ovL),]$stat
  annot.all$pos.cumul <- out.gr[subjectHits(ovL),]$pos.cumul
  Q <- quantile(annot.all$stat, 0.99, na.rm = TRUE)
  
  bool <- y>Q
  bool[is.na(bool)] <- FALSE
  bool.seq <- ave(bool, cumsum(!bool), FUN = cumsum)<
    
    THRES <- 1
  bool.all <- bool.seq >= THRES
  id.bool <- which(bool.seq == THRES)
  bool.all[c(id.bool - 1)] <- TRUE
  bool.all <- as.vector(bool.all)
  ovL <- findOverlaps(annot.all, out.gr[bool.all,])
  annot.all <- annot.all[queryHits(ovL), ]
  out <- annot.all$stat
  names(out) <- annot.all$name
  return(out)
}


psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}

