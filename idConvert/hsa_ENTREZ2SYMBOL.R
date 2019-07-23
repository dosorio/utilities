 hsa_ENTREZ2SYMBOL<- function(X){
  require(AnnotationHub)
  require(AnnotationDbi)
  db <- AnnotationHub()
  hsa <- db[["AH70572"]]
  eIds <- select(hsa,X,"SYMBOL","ENTREZID")
  eIds <- eIds[complete.cases(eIds),]
  return(eIds)
}
