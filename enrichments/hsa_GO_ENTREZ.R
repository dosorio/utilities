hsa_GO_ENTREZ <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  db <- AnnotationHub()
  hsa <- db[["AH70572"]
  out <- clusterProfiler::enrichGO(gene = X, keyType = "ENTREZID", OrgDb = hsa)
  as.data.frame(out)
}
