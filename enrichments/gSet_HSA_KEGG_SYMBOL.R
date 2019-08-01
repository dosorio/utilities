gSet_HSA_KEGG_SYMBOL <- function(){
  hsaKEGG <- gSets_KEGG("hsa")
  hsaKEGG <- lapply(hsaKEGG, function(X){
    hsa_ENTREZ2SYMBOL(X)[,2]
  })
  return(hsaKEGG)
}
