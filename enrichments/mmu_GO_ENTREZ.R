mmu_GO_ENTREZ <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  db <- AnnotationHub()
  mmu <- db[["AH70573"]]
  out <- clusterProfiler::enrichGO(gene = X, keyType = "ENTREZID", OrgDb = mmu, ont = "ALL")
  as.data.frame(out)
}
