cleanNet <- function(X){
  X <- as.matrix(X[,])
  X <- X[apply(X,1,sum) > 0, apply(X,2,sum) > 0]
  X <- reshape2::melt(X)[,1:2]
  colnames(X) <- c("from", "to")
  X <- igraph::graph_from_data_frame(X)
  return(X)
}
