hsa_GO_ENTREZ <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  out <- clusterProfiler::enrichGO(gene = X, keyType = "ENTREZID", OrgDb = hsa)
  as.data.frame(out)
}
