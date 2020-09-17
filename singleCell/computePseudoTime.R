computePseudoTime <- function(cMatrix, simplified = TRUE, nDim = 100){
  #nDim <- ifelse(nDim < ncol(cMatrix), ncol(cMatrix)-1, nDim)
  cMatrix <- cMatrix[rowSums(cMatrix) != 0,]
  if(isTRUE(simplified)){
    require(RSpectra)
    require(Matrix)
    nMatrix <- Matrix::t(cMatrix)
    nMatrix <- nMatrix/rowSums(cMatrix) * 1e4
    nMatrix <- log1p(nMatrix)
    nMatrix <- Matrix::t(scale(nMatrix))
    nMatrix <- Matrix::t(svds(nMatrix, nDim)$v)
    colnames(nMatrix) <- colnames(cMatrix)
    rownames(nMatrix) <- paste0('g', seq_len(nDim))
    cMatrix <- nMatrix
    remove(nMatrix)
  }
  require(monocle)
  fd <- data.frame('gene_short_name' = rownames(cMatrix))
  rownames(fd) <- rownames(cMatrix)
  fd <- new("AnnotatedDataFrame", data = fd)
  cds <- newCellDataSet(as.matrix(cMatrix), featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 2)
  cds <- orderCells(cds)
  o <- pData(cds)[,2]
  attr(o, which = 'DDRTree') <- t(cds@reducedDimS)
  colnames(attr(o, which = 'DDRTree')) <- paste0('DDRTree', 1:2)
  return(o)
}
