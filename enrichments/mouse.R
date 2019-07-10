makeMouseEnrichmentGO <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  db <- AnnotationHub()
  mmu <- db[["AH70573"]]
  eIds <- select(mmu,X, "ENTREZID", "SYMBOL")[,2]
  eOut <- clusterProfiler::enrichGO(gene = eIds, keyType = "ENTREZID", OrgDb = mmu)
  as.data.frame(eOut)
}

makeMouseEnrichmentKEGG <- function(X){
  require(clusterProfiler)
  require(AnnotationHub)
  db <- AnnotationHub()
  mmu <- db[["AH70573"]]
  eIds <- select(mmu,X, "ENTREZID", "SYMBOL")[,2]
  eOut <- clusterProfiler::enrichKEGG(gene = eIds,organism = "mmu")
  as.data.frame(eOut)
}
