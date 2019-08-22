pcNet <- function(X, nCom = 3, scaleScores = TRUE, symmetric = FALSE, q = 0){
  gNames <- rownames(X)
  X <- (scale(t(X)))  
  n <- ncol(X)
  A <- 1-diag(n)
  getCoefficients <- function(K){
    y <- X[,K]
    Xi <- X
    Xi <- Xi[,-K]
    coeff <- RSpectra::svds(Xi, nCom)$v
    score <- Xi %*% coeff
    score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    return(Beta)
  }
  B <- pbapply::pbsapply(seq_len(n), getCoefficients)
  B <- t(B)
  for(K in seq_len(n)){
    A[K,A[K,] == 1] = B[K,]
  }
  if(isTRUE(symmetric)){
    A <- (A + t(A))/2  
  }
  absA <- abs(A)
  if(isTRUE(scaleScores)){
    A <- (A/max(absA))  
  }
  A[absA < quantile(absA,q)] <- 0
  diag(A) <- 0
  colnames(A) <- rownames(A) <- gNames
  A <- Matrix::Matrix(A, sparse = TRUE)
  return(A)
}
