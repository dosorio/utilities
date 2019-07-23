hsa_ENSEMBL2SYMBOL <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  db <- AnnotationHub()
  hsa <- db[["AH70572"]]
  eIds <- select(hsa,X, "SYMBOL", "ENSEMBL")
  eIds <- eIds[complete.cases(eIds),]
  return(eIds)
}
