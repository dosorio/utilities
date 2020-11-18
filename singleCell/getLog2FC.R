# Takes a Seurat object and computes the log2 Fold Change for each cluster compared with the others.
getLog2FC <- function(X){
  require(Matrix)
  allIdent <- levels(Idents(X))
  countData <- X@assays$RNA@counts
  countData <- (t(t(countData)/colSums(countData)))*1e4
  lapply(allIdent, function(i){
    iMean <- rowMeans(countData[,Idents(X) %in% i])
    oMean <- rowMeans(countData[,!(Idents(X) %in% i)])
    sGenes <- (iMean > 0) & (oMean > 0)
    iMean <- iMean[sGenes]
    oMean <- oMean[sGenes]
    FC <- iMean/oMean
    log2FC <- log2(FC)
    return(log2FC)
  })
}
