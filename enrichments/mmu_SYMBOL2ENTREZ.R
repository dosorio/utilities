mmu_SYMBOL2ENTREZ <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  db <- AnnotationHub()
  mmu <- db[["AH70573"]]
  eIds <- select(mmu,X, "ENTREZID", "SYMBOL")
  eIds <- eIds[complete.cases(eIds),]
  return(eIds)
}
