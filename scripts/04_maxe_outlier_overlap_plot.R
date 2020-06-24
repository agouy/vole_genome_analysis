# plot overlap between gene lists
pis <- res[, grep("pi", colnames(res))]
tab <- list()
for(i in 1:12) {
  x <- rndsps[[i]]
  x[cond] <- NA
  # getGeneList(x, res, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99)
  # getGeneList(x, res, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99)
  
  pi.out <- getGeneDiv(rndsps[[i]], res, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99, 1)
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

# df <- apply(df, 2, "*", 1:nrow(df))
load("C:/PhD/def-par.rda")
par(pars)
par(mfrow=c(1,1), mar=c(4,4,4,0))
plot(NA, xlim = c(1, ncol(df) + 2), ylim = c(1, nrow(df)), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1, 1:nrow(df), 12, 1:nrow(df), col = rgb(.9,.9,.9))
segments(1:12, 1, 1:12, nrow(df), col = rgb(.95,.95,.95))
for(i in 1:12) {
  y <- df[,i]
  points(rep(i, sum(!is.na(y))), which(!is.na(y)), pch = 15, cex = 1.5, col = rgb(y[!is.na(y)], 0.5, 0.8*(1-y[!is.na(y)]), 0.9))
}
axis(1, at = 1:ncol(df), labels = samp.names[-9][ind.ord], lwd = 0, las = 2, cex.axis = 0.8)
axis(2, at = 1:nrow(df), labels = rownames(df), lwd = 0, las = 2, cex.axis = 0.5)

Y.le <- round(0.9*nrow(df))
text(c(13.5), 1.1*Y.le, expression(paste(pi, " quantile")), font = 2, cex = 0.8)
points(c(13,13.5,14), rep(1.05*Y.le, 3), col =  rgb(c(0,0.5,1), 0.5, 0.8*(1-c(0,0.5,1))), pch = 15, cex = 2)
text(c(13,13.5,14), rep(1.05*Y.le, 3), c(0,0.5,1), pos = 1, cex = 0.6)

title("Outlier olfactory and vomeronasal receptors", adj = 0.4, cex.main = 0.9, outer = F)

dfolf <- df


# plot overlap between gene lists
tab <- list()
for(i in 1:12) {
  x <- rndsps[[i]]
  x[cond] <- NA
  # getGeneList(x, res, annot.gr, quanti = 0.)
  # getGeneList(x, res, annot.gr[c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99)
  
  # pi.out <- getGeneDiv(rndsps[[i]], res, annot.gr[-c(grep("Olf", annot.gr$name), grep("Vmn", annot.gr$name)),], quanti = 0.99, 1)
  pi.out <- getGeneDiv(rndsps[[i]], res, annot.gr, quanti = 0.99, 1)
  qt.pi <- colMeans(sapply(pi.out, ">", pis[,i]))
  tab[[i]] <- qt.pi[!is.na(qt.pi)]
  print(i)
}
uniq.name <- unique(unlist(lapply(tab, function(x) names(x))))
df <- matrix(NA, ncol=length(tab), nrow=length(uniq.name))
rownames(df) <- uniq.name
for(i in 1:12) df[names(tab[[i]]), i] <- tab[[i]]

dfall <- df


p.ord <- order(rowSums(!is.na(df)), decreasing = TRUE)
ind.ord <- c(4,1,2,3,12,6,11,10,5,7,8,9)
df <- df[p.ord, ][, ind.ord]

# df <- apply(df, 2, "*", 1:nrow(df))
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(NA, xlim = c(1, ncol(df) + 2), 
     ylim = c(1, nrow(df)), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
segments(1, 1:nrow(df), 12, 1:nrow(df), col = rgb(.9,.9,.9))
segments(1:12, 1, 1:12, nrow(df), col = rgb(.95,.95,.95))
for(i in 1:12) {
  y <- df[,i]
  points(rep(i, sum(!is.na(y))), which(!is.na(y)), pch = 15, cex = .5, col = rgb(y[!is.na(y)], 0.5, 0.8*(1-y[!is.na(y)]), 0.9))
}
axis(1, at = 1:ncol(df), labels = samp.names[-9][ind.ord], lwd = 0, las = 2, cex.axis = 0.8)
axis(2, at = 1:nrow(df), labels = rownames(df), lwd = 0, las = 2, cex.axis = 0.5)

Y.le <- round(max(df3, na.rm = TRUE) / 2)
text(c(13.5), 1.1*Y.le, "Q(pi)", font = 2, cex = 0.8)
points(c(13,13.5,14), rep(Y.le, 3), col =  rgb(c(0,0.5,1), 0.5, 0.8*(1-c(0,0.5,1))), pch = 15, cex = 2)
text(c(13,13.5,14), rep(Y.le, 3), c(0,0.5,1), pos = 1, cex = 0.6)

title("Outlier olfactory and vomeronasal receptors", adj = 0.4, cex.main = 0.9, outer = F)



par(mfrow=c(4, 3), mar=c(4,4,2,2))
# for (i in 1:12) {
#   hist(dfall[,i], xlim = c(0, 1), xlab = "Q(pi)", main = samp.names[-9][ind.ord][i], breaks = 100)
# }
# for (i in 1:12) {
#   hist(dfolf[,i], xlim = c(0, 1), xlab = "Q(pi)", main = samp.names[-9][ind.ord][i], breaks = 100)
# }
for (i in 1:12) {
  qqplot(dfolf[,i], dfall[,i], main = samp.names[-9][ind.ord][i], xlim = c(0,1), ylim = c(0,1),
         xlab = "OR quantiles", ylab = "Non OR quantiles", col = 1, pch = 16, cex = 0.5, asp=1)
  par(new=TRUE)
  qqplot(dfolf[,i], ppoints(500), main = "", xlim = c(0,1), ylim = c(0,1),
          xlab = "", ylab = "", add = TRUE, col = 2, pch = 16, cex = 0.5, asp=1)
  abline(0,1)
}




dfneut <- df