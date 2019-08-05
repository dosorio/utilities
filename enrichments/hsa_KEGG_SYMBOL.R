hsa_KEGG_SYMBOL <- function(X){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa_SYMBOL2ENTREZ.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_KEGG_ENTREZ.R")
  hsa_KEGG_ENTREZ(hsa_SYMBOL2ENTREZ(X)[,2])
}
