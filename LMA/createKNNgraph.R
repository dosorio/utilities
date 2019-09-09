createKNNgraph <- function(X){
  n1 <- ncol(X)
  W1 <- Matrix::Matrix(0, nrow = n1, ncol = n1)
  for(i in seq_len(n1)){
    W1[i,X[,i]] <- 1
  }
  W <- (W1 + t(W1))/2
  return(W)
}
