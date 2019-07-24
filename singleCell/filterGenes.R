filterGenes <- function(X){
  X <- X[apply(X!=0,1,sum) > 0.1,]
  X <- X[apply(X,1,var) > 0,]
  return(X)
}
