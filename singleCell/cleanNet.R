cleanNet <- function(X){
    X <- as.matrix(X[,])
    X <- reshape2::melt(X)
    colnames(X) <- c("from", "to", "weight")
    X <- X[abs(X[,3]) > 0,]
    X <- igraph::graph_from_data_frame(X)
    return(X)
  }
