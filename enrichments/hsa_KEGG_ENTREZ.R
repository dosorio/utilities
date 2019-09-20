hsa_KEGG_ENTREZ <- function(X, ...){
  require(clusterProfiler)
  require(AnnotationHub)
  eOut <- clusterProfiler::enrichKEGG(gene = X,organism = "hsa", ...)
  as.data.frame(eOut)
}
