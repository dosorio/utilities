getPseudoTime <- function(cMatrix, simplified = TRUE, nDim = 100){
  cMatrix <- cMatrix[rowSums(cMatrix) > 0,]
  if(isTRUE(simplified)){
    require(RSpectra)
    require(Matrix)
    nMatrix <- log1p(t(t(cMatrix)/colSums(cMatrix)) * 1e4)
    nMatrix <- scale(t(nMatrix))
    nMatrix <- t(svds(nMatrix, nDim)$v)
    colnames(nMatrix) <- colnames(cMatrix)
    rownames(nMatrix) <- paste0('g', seq_len(nDim))
    cMatrix <- nMatrix
  }
  require(monocle)
  fd <- data.frame('gene_short_name' = rownames(cMatrix))
  rownames(fd) <- rownames(cMatrix)
  fd <- new("AnnotatedDataFrame", data = fd)
  cds <- newCellDataSet(as.matrix(cMatrix), featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 2)
  cds <- orderCells(cds)
  pseudoTiveV <- pData(cds)
  return(pseudoTiveV[,2])
}
