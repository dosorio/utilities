getCellClusters <- function(X){
  X <- CreateSeuratObject(as.matrix(X))
  X <- NormalizeData(X)
  X <- ScaleData(X)
  X <- FindVariableFeatures(X)
  X <- RunPCA(X, verbose = FALSE)
  X <- FindNeighbors(X)
  X <- FindClusters(X, resolution = 0.8)
  return(Idents(X))
}
