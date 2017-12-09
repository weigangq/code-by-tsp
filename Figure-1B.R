# Plot TSP path for SGC (Figure 1B)
aa.feat<-read.table("../Google Drive/Monte-Carlo-Club-2017/aa.txt2", sep="\t", header = T, row.names = 1)
paths <- read.table("../Google Drive/Monte-Carlo-Club-2017/codon-paths.txt", sep = "\t", header=F)
path.list <-split(paths, paths[,1])
path.pol.df <- path.list$'CTC' # highest cor

# SGC
library(Biostrings);
cd <- GENETIC_CODE;
cd[which(cd=="*")] <- 'X'
codes <- names(cd);
mat <- matrix(NA, nrow = 21, ncol = 64);
colnames(mat) <- path.pol.df[,2];
col.idx <- grep(paste("^p", sep=""), colnames(aa.feat));
aa.ordered <- aa.feat[order(aa.feat[,col.idx]),col.idx];
names(aa.ordered) <- rownames(aa.feat)[order(aa.feat[,col.idx])];
rownames(mat) <- names(aa.ordered);
  for(r in 1:nrow(mat)) {
    for (c in 1:ncol(mat)) {
      mat[r,c] <- ifelse(cd[colnames(mat)[c]]==rownames(mat)[r], 1, 0)
    }
  }

map.aa <- integer(64);
map.codon <- integer(64);
k <- 1;
for(c in 1:ncol(mat)) {
  for (r in 1:nrow(mat)) {
    if (mat[r,c]==1) {
      map.aa[k] <- c;
      map.codon[k] <- r;
      k <- k + 1;
    }
  }
}

# shuffle code by 21 synonymous blocks
cd.ran.df <- read.table("../Google Drive/Monte-Carlo-Club-2017/sgc.s4", sep = "\t", row.names = 1, header=F)
cd.ran <- as.character(cd.ran.df[,1])
names(cd.ran) <- rownames(cd.ran.df)
mat.ran <- matrix(NA, nrow = 21, ncol = 64);
colnames(mat.ran) <- path.pol.df[,2];
rownames(mat.ran) <- names(aa.ordered);
for(r in 1:nrow(mat.ran)) {
  for (c in 1:ncol(mat.ran)) {
    mat.ran[r,c] <- ifelse(cd.ran[colnames(mat.ran)[c]]==rownames(mat.ran)[r], 1, 0)
  }
}

map.aa.ran <- integer(64);
map.codon.ran <- integer(64);
k <- 1;
for(c in 1:ncol(mat.ran)) {
  for (r in 1:nrow(mat.ran)) {
    if (mat.ran[r,c]==1) {
      map.aa.ran[k] <- c;
      map.codon.ran[k] <- r;
      k <- k + 1;
      #        points(c, r, col=2, pch=16)
    }
  }
}
  
# plot paths
plot(1:10, 1:10, xlim=c(0,65), ylim=c(0,23), type="n", las=1, xlab="Codon", ylab="AA", main="SGC", cex.main=0.8)
abline(h=1:21, col="gray")
abline(v=1:64, col="gray")
text(1:64, rep(22,64), colnames(mat), cex=0.75, srt=90)
text(rep(0, 21), 1:21, rownames(mat), cex=0.75)
lines(map.aa, map.codon, type = "b", pch=16)
lines(map.aa.ran, map.codon.ran, type = "b", pch=1, col=2, lty=2)
#legend(20,10, c("SGC", "Permutated"), lty=c(1,2), col=c(1,2), pch=c(16,1), cex=0.75)
  
    
