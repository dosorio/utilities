scQC <- function(X){
  lSize <- apply(X,2,sum)
  lLimits <- quantile(lSize, c(0.15,0.95))
  mtCounts <- apply(X[grepl("^MT",toupper(rownames(X))),,drop=FALSE],2,sum)
  mtRate <- mtCounts/lSize
  #X <- t(t(X)/lSize)*1e6
  X <- X[,lSize > lLimits[1] & lSize < lLimits[2] & mtRate < 0.1]
  X <- X[apply(X != 0, 1 ,mean) > 0.1,]
  return(X)
}
