hsa_GO_SYMBOL <- function(X){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa_SYMBOL2ENTREZ.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_ENTREZ.R")
  hsa_GO_ENTREZ(hsa_SYMBOL2ENTREZ(X)[,2])
}
