hsa_SYMBOL2ENTREZ <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  db <- AnnotationHub()
  hsa <- db[["AH70572"]]
  eIds <- select(hsa,X, "ENTREZID", "SYMBOL")
  eIds <- eIds[complete.cases(eIds),]
  return(eIds)
}
