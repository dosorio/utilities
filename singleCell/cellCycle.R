 cellCycle <- function(X){
  require(Seurat)
  X <- CreateSeuratObject(X)
  X <- NormalizeData(X)
  X <- ScaleData(X)
  set.seed(1)
  X <- CellCycleScoring(X, 
                        s.features = cc.genes$s.genes, 
                        g2m.features = cc.genes$g2m.genes)
  X <- as.vector(X$Phase)
  return(X)
}
