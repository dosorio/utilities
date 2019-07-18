mmu_GO_ENTREZ <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  out <- clusterProfiler::enrichGO(gene = X, keyType = "ENTREZID", OrgDb = mmu)
  as.data.frame(out)
}
