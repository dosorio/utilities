sccNet <- function(X, q = 0.95, nCell = 500, nNet = 50, denoiseNet = TRUE){
  nGenes <- nrow(X)
  gList <- rownames(X)
  oNet <- pbapply::pbsapply(seq_len(nNet), function(Z){
    tNet <- Matrix::t(X[,sample(colnames(X), nCell)])
    tNet <- cor(as.matrix(tNet), method = 'sp')
    tNet <- tNet/max(abs(tNet))
    tNet <- round(tNet,1)
    diag(tNet) <- 0
    tNet[abs(tNet) < quantile(abs(tNet), q)] <- 0
    tNet <- Matrix::Matrix(tNet)
    return(tNet)  
  })
  
  aNet <- matrix(0, nGenes, nGenes)
  rownames(aNet) <- colnames(aNet) <- gList
  for(i in seq_along(oNet)){
    tNet <- oNet[[i]]
    if(denoiseNet){
      tNet <- RSpectra::svds(tNet,2)
      tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)  
      rownames(tNet) <- colnames(tNet) <- gList
    }
    tNet <- Matrix::Matrix(tNet)
    aNet <- aNet + tNet[gList, gList]
  }
  aNet <- aNet/length(oNet)
  aNet <- aNet/abs(max(aNet))
  aNet <- round(aNet,1)
  aNet <- Matrix::Matrix(aNet)
  return(aNet)
}
