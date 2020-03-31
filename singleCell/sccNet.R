sccNet <- function(X, q = 0.95, nCell = 500, nNet = 25, K = 2, denoiseNet = TRUE){
  X <- t(t(X)/colSums(X))
  nGenes <- nrow(X)
  gList <- rownames(X)
  set.seed(1)
  oNet <- pbapply::pbsapply(seq_len(nNet), function(Z){
    tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
    tNet <- cor(as.matrix(tNet), method = 'sp')
    while(any(is.na(tNet))){
      tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
      tNet <- cor(as.matrix(tNet), method = 'sp')
    }
    tNet <- round(tNet,1)
    diag(tNet) <- 0
    tNet[abs(tNet) < quantile(abs(tNet), q, na.rm = TRUE)] <- 0
    tNet <- Matrix::Matrix(tNet)
    return(tNet)  
  })
  
  aNet <- matrix(0, nGenes, nGenes)
  rownames(aNet) <- colnames(aNet) <- gList
  for(i in seq_along(oNet)){
    tNet <- oNet[[i]]
    if(denoiseNet){
      if(K == ncol(X)){
        tNet <- svd(tNet)
        tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)  
        rownames(tNet) <- colnames(tNet) <- gList  
      } else if(K == 1){
        tNet <- RSpectra::svds(tNet,K)
        tNet <- tNet$u %*% (tNet$d) %*% t(tNet$v)  
        rownames(tNet) <- colnames(tNet) <- gList  
      } else {
        tNet <- RSpectra::svds(tNet,K, maxitr = 1e6)
        tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)  
        rownames(tNet) <- colnames(tNet) <- gList  
      }
      
    }
    tNet <- Matrix::Matrix(tNet)
    aNet <- aNet + tNet[gList, gList]
  }
  aNet <- aNet/length(oNet)
  aNet <- Matrix::Matrix(aNet)
  return(aNet)
}
