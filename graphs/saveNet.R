saveNet <- function(X, outFile){
  require(igraph)
  X <- as.matrix(X[,])
  X <- X[rowSums(X) > 0,colSums(X) > 0]
  write.csv(X, quote = FALSE, file = outFile)
}
