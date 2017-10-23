# Plot TSP path for SGC (Figure 1B)
aa.feat<-read.table("../Google Drive/Monte-Carlo-Club-2017/aa.txt2", sep="\t", header = T, row.names = 1)
#path <- scan("../Google Drive/Monte-Carlo-Club-2017/codon-path.TGT", what = "ch")
path.pol.df <- path.pol.list$'CTC' # highest cor
path.hydro.df <- path.hydro.list$'CTA';

library(Biostrings);
cd <- GENETIC_CODE;
cd[which(cd=="*")] <- 'X'
codes <- names(cd);
mat <- matrix(NA, nrow = 21, ncol = 64);
#rownames(mat) <- sort(names(cd));
rownames(mat) <- path.hydro.df[,2];
plot.code <- function (feat = "p", title="SGC") {
  col.idx <- grep(paste("^", feat, sep=""), colnames(aa.feat));
  #  cat(col.idx);
  aa.ordered <- aa.feat[order(aa.feat[,col.idx]),col.idx];
  names(aa.ordered) <- rownames(aa.feat)[order(aa.feat[,col.idx])];
  colnames(mat) <- names(aa.ordered);
  for(r in 1:nrow(mat)) {
    for (c in 1:ncol(mat)) {
      mat[r,c] <- ifelse(cd[rownames(mat)[r]]==colnames(mat)[c], 1, 0)
    }
  }
  
  plot(1:10, 1:10, xlim=c(0,23), ylim=c(0,75), type="n", las=1, xlab="AA", ylab="Codon", main=title, cex.main=0.8)
  abline(h=1:64, col="gray")
  abline(v=1:21, col="gray")
  text(1:21, rep(66,21), colnames(mat))
  text(rep(22, 64), 1:64, rownames(mat), cex=0.5)
  for(r in 1:nrow(mat)) {
    for (c in 1:ncol(mat)) {
      if (mat[r,c]==1) {
        points(c, r, col=2, pch=16)
      }
    }
  }
  
  abline(h=68)
  #  segments(1:21, rep(68,21), 1:21, polar.ordered/2 + 68, col=4, lwd=4)
  # points(1:21, polar.ordered/2 + 68, col=4, pch=16)
  text(3,75, paste("odered by", colnames(aa.feat)[col.idx]), col=4)
}