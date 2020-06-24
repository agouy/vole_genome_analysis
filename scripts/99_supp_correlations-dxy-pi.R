## LOAD res object

dxys <- res[, grep("dxy", colnames(res))]
head(dxys)

todo <- combn(1:13, 2)

dxnames <- apply(todo, 2, function(x) {
  paste0("dxy.pop", x[1], ".pop", x[2])
})

cor.mat <- cor(dxys, method = "spearman")
dx.mean <- colMeans(dxys)

cor.2 <- cor.mat
cor.2[cor.2 > -1] <- NA
for(i in rownames(cor.mat)) {
  for(j in colnames(cor.mat)) {
    nms <- combn(as.numeric(c(unlist(strsplit(gsub("dxy.pop", "", i), ".pop")),
                              unlist(strsplit(gsub("dxy.pop", "", j), ".pop")))), 2)
    nms <- nms[,1-colSums(apply(nms, 2, duplicated))]
    
    
    maxi <- dx.mean[paste0("dxy.pop", nms[1,], ".pop", nms[2,])]
    cor.2[i,j] <- maxi[which.max(maxi)]
  }
}

# per species
# plot dxy sp focale vs other as a function od dxy

for(i in 1:13) {
  print(i)
  abc <- 1:13
  cond <- names(dx.mean) %in% c(paste0("dxy.pop", i, ".pop", abc[-i]),  paste0("dxy.pop", abc[-i], ".pop", i))
  cor.sp <- cor.mat[rownames(cor.mat) %in% names(dx.mean)[cond], colnames(cor.mat) %in% names(dx.mean)[cond]]

  matplot(dx.mean[cond], t(cor.sp), type = "p", pch = 16, col = rainbow(20))
}

str(rnds)
rnds <- getRND(dxys)

### correlations between nucleotide diversities

pis <- res[, grep("pi", colnames(res))]
cor.pi <- cor(pis, method= "spearman")
heatmap(cor.pi)

dx.mean <- colMeans(dxys)
mat <- matrix(NA, ncol=dim(cor.pi)[1], nrow=dim(cor.pi)[1])
mat[lower.tri(mat)] <- dx.mean
mat[upper.tri(mat)] = t(mat)[upper.tri(mat)]
mat[is.na(mat)] <- 0

coul <- rainbow(13)
for(i in 1:13) {
  if(i == 1) plot(mat[i,], cor.pi[i,], col = coul[i], pch = 16)
  else points(mat[i,], cor.pi[i,], col = coul[i], pch = 16)
}
matplot(cor.pi, mat)

diag(mat) <- NA
diag(cor.pi) <- NA

par(mfrow=c(2,1))
plot((mat), (cor.pi), pch = 16, cex = 1.2, bty = "n", xlab = "Average divergence", 
     ylab = "Correlation between diversities")

newx <- seq(0.0001,0.004,0.00001)
pi2 <- as.dist(mat)
lines(newx, predict(loess(as.dist(cor.pi) ~ pi2, family = "symmetric", span = 1), newdata = data.frame(pi2 = newx)),
      type = "l", lwd = 3, lty = 1, col = "firebrick")
title("Correlations between divergence and diversity")


layout(matrix(c(1,2,3,4,3,4), nrow = 2))
hist(cor.mat, breaks = 50, main = "", xlab = expression(d[xy]))
title("a.", adj = 0)

plot(cor.2, cor.mat, type = "p", pch = 16, col = rgb(0,0,0,0.3), bty = "n", xlab = "Average divergence",
     ylab = "Correlation between dxys")
newx <- seq(0, 0.015, 0.00001)
lines(newx, predict(loess(c(cor.mat)~c(cor.2), span = 1, family = "symmetric"), 
                    newdata=data.frame(cor.2 = newx)), type = "l", lwd = 3, lty = 1, col = "firebrick")
title("b.", adj = 0)

plotOut2(dxys[,"dxy.pop5.pop6"], ylim=c(0,0.015), g.nam = F, ylab = expression(d[xy]))
title("c. dxy - M. agrestis vs. M. lusitanicus", adj = 0)
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.9)

plotOut2(dxys[,"dxy.pop7.pop8"], ylim=c(0,0.02), g.nam = F, ylab = expression(d[xy]))
title("d. dxy - M. oeconomus vs. M. brandtii", adj = 0)
ax <- tapply(res$pos.cumul, res$sc, mean)
axis(1, at = ax, labels = 1:22, lwd.ticks = 1, lwd = 0, cex.axis = 0.9)


###

tt <- apply(dxys, 2, function(x) {
  ret <- loess(x ~ res$pos.cumul, span = 0.01)
  ret <- predict(ret)
  return(ret)
})

matplot(res$pos.cumul, tt, type = "l")
pairs(tt[,1:20])
rnds <- getRND(res, 9)
cor.rnd <- cor(rnds, use = "complete.obs", method = "spearman")

par(mfrow=c(1,1))
boxplot(c(cor.pi), c(cor.mat), c(cor.rnd), main = "Correlation coefficients distributions",
        names = c("Pi", "Dxy", "RND"))
points(1+rnorm(length(cor.pi),0,0.1), c(cor.pi), pch = 16, col = rgb(0,0,0,0.2))
points(2+rnorm(length(cor.mat),0,0.1), c(cor.mat), pch = 16, col = rgb(0,0,0,0.1))
points(3+rnorm(length(cor.rnd),0,0.1), c(cor.rnd), pch = 16, col = rgb(0,0,0,0.1))

