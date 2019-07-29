  makeGraph <- function(X, qThreshold = 0.9){
    X <- as.matrix(X)
    diag(X) <- NA
    qT <- quantile(abs(X), qThreshold, na.rm = TRUE)
    X[abs(X) < qT] <- NA
    X <- reshape2::melt(X)
    X <- X[complete.cases(X), ]
    X <- igraph::graph_from_data_frame(X)
    return(X)
  }
