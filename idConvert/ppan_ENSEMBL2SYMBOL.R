ppan_ENSEMBL2SYMBOL <- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  require(ensembldb)
  db <- AnnotationHub()
  ppan <- db[["AH73946"]]
  eIds <- select(ppan,X, "SYMBOL", "GENEID")
  eIds <- eIds[eIds[,2] != "",]
  return(eIds)
}
