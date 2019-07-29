hsa_ENSEMBL2SYMBOL <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  db <- AnnotationHub()
  hsa <- db[["AH73881"]]
  eIds <- select(hsa,X, "SYMBOL", "GENEID")
  eIds <- eIds[eIds[,2] != "",]
  return(eIds)
}
