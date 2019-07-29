ptr_ENSEMBL2SYMBOL <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  require(ensembldb)
  db <- AnnotationHub()
  ptr <- db[["AH73952"]]
  eIds <- select(ptr,X, "SYMBOL", "GENEID")
  eIds <- eIds[eIds[,2] != "",]
  return(eIds)
}
