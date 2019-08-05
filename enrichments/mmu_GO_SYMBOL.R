mmu_GO_SYMBOL <- function(X){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/mmu_SYMBOL2ENTREZ.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/mmu_GO_ENTREZ.R")
  mmu_GO_ENTREZ(mmu_SYMBOL2ENTREZ(X)[,2])
}
