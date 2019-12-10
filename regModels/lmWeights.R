lmWeights <- function(X){
   X <- cbind(1,X)
   X <- X %*% solve((t(X)%*%X)) %*% t(X)
   X <-  diag(X)
   return(X)
}
