gSet_HSA_KEGG_SYMBOL <- function(){
  source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/gSets_KEGG.R')
  source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa_ENTREZ2SYMBOL.R')
  hsaKEGG <- gSets_KEGG("hsa")
  hsaKEGG <- lapply(hsaKEGG, function(X){
    hsa_ENTREZ2SYMBOL(X)[,2]
  })
  return(hsaKEGG)
}
