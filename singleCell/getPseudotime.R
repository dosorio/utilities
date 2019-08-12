getPseudoTime <- function(cMatrix){
  cMatrix <- cMatrix[rowSums(cMatrix) > 0,]
  require(monocle)
  fd <- data.frame('gene_short_name' = rownames(cMatrix))
  rownames(fd) <- rownames(cMatrix)
  fd <- new("AnnotatedDataFrame", data = fd)
  cds <- newCellDataSet(as.matrix(cMatrix), featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 2)
  cds <- orderCells(cds)
  pseudoTiveV <- pData(cds)
  return(pseudoTiveV)
}
