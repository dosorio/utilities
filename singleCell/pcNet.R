pcNet <- function(X, nCom = 3, nCores = 1){
  require(RSpectra)
  gNames <- rownames(X)
  X <- (scale(t(X)))
  n <- ncol(X)
  A <- 1-diag(n)
  if(nCores > 1){
    require(parallel)
    cl <- makeCluster(getOption("cl.cores", nCores))
    clusterExport(cl,"X", envir = environment())
    clusterExport(cl,"nCom", envir = environment())
    B <- parSapply(cl, seq_len(n), function(K){
      y <- X[,K]
      Xi <- X
      Xi <- Xi[,-K]
      coeff <- RSpectra::svds(Xi, nCom)$v
      score <- Xi %*% coeff
      score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
      Beta <- colSums(y * score)
      return(coeff %*% (Beta))
    })
    stopCluster(cl)
  } else {
    B <- sapply(seq_len(n), function(K){
      y <- X[,K]
      Xi <- X
      Xi <- Xi[,-K]
      coeff <- RSpectra::svds(Xi, nCom)$v
      score <- Xi %*% coeff
      score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
      Beta <- colSums(y * score)
      return(coeff %*% (Beta))
    })
  }
  B <- t(B)
  for(K in seq_len(n)){
    A[K,A[K,] == 1] = B[K,]
  }
  diag(A) <- 1
  A <- round(A,5)
  colnames(A) <- rownames(A) <- gNames
  return(A)
}
