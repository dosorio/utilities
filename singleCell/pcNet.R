pcNet <- function(X, nCom = 3, nCores = 1){
  gNames <- rownames(X)
  X <- (scale(t(X)))
  n <- ncol(X)
  A <- 1-diag(n)
  if(nCores > 1){
    if(grepl("Windows", Sys.info()[[1]])){
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
    } else {
      cl <- parallel::makeCluster(getOption("cl.cores", nCores),  type = "FORK")
    }
    parallel::clusterExport(cl,"X", envir = environment())
    parallel::clusterExport(cl,"nCom", envir = environment())
    B <- pbapply::pbsapply(seq_len(n), function(K){
      y <- X[,K]
      Xi <- X
      Xi <- Xi[,-K]
      coeff <- RSpectra::svds(Xi, nCom)$v
      if(class(coeff) != 'try-error'){
        score <- Xi %*% coeff
        score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
        Beta <- colSums(y * score)
        return(coeff %*% (Beta))
      } else{
        return(rep(NA, ncol(Xi)))
      }
    }, cl = cl)
    parallel::stopCluster(cl)
  } else {
    B <- pbapply::pbsapply(seq_len(n), function(K){
      y <- X[,K]
      Xi <- X
      Xi <- Xi[,-K]
      coeff <- try(RSpectra::svds(Xi, nCom)$v, silent = TRUE)
      if(class(coeff) != 'try-error'){
        score <- Xi %*% coeff
        score <- t(t(score)/(apply(score,2,function(X){sqrt(sum(X^2))})^2))
        Beta <- colSums(y * score)
        return(coeff %*% (Beta))
      } else{
        return(rep(NA, ncol(Xi)))
      }
    })
  }
  B <- t(B)
  for(K in seq_len(n)){
    A[K,A[K,] == 1] = B[K,]
  }
  diag(A) <- 1
  colnames(A) <- rownames(A) <- gNames
  return(A)
}
