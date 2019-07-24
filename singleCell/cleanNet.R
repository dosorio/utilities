cleanNet <- function(X){
  X <- as.matrix(X[,])
  X <- X[apply(X,1,sum) > 0, apply(X,2,sum) > 0]
}
