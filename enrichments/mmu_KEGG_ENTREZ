hsa_KEGG_ENTREZ <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  eOut <- clusterProfiler::enrichKEGG(gene = X,organism = "mmu")
  as.data.frame(eOut)
}
