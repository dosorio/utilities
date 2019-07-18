mmu_ENTREZ2SYMBOL <- function(X){
  require(AnnotationHub)
  db <- AnnotationHub()
  mmu <- db[["AH70573"]]
  select(x = mmu, keys = X, columns = "SYMBOL", keytype = "ENTREZID")
}
