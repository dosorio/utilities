pcNet <- function(X, nCom = 3){
  gNames <- rownames(X)
  X <- (scale(t(X)))
  n <- ncol(X)
  A <- 1-diag(n)
  getCoefficients <- function(K){
    y <- X[,K]
    Xi <- X
    Xi <- Xi[,-K]
    coeff <- try(RSpectra::svds(Xi, nCom)$v, silent = TRUE)
    if(class(coeff) != 'try-error'){
      score <- Xi %*% coeff
      score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
      Beta <- colSums(y * score)
      Beta <- coeff %*% (Beta)
      Beta <- round(Beta,5)
      return(Beta)
    } else{
      return(rep(NA, ncol(Xi)))
    }
  }
  B <- pbapply::pbsapply(seq_len(n), getCoefficients)
  B <- t(B)
  for(K in seq_len(n)){
    A[K,A[K,] == 1] = B[K,]
  }
  diag(A) <- 1
  colnames(A) <- rownames(A) <- gNames
  return(A)
}
