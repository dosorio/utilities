mmu_KEGG_SYMBOL <- function(X){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/mmu_SYMBOL2ENTREZ.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/mmu_KEGG_ENTREZ.R")
  mmu_KEGG_ENTREZ(mmu_SYMBOL2ENTREZ(X)[,2])
}
