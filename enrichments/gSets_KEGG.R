gSets_KEGG <- function(sp){
  require(KEGGREST)
  keggSP <- keggList("pathway",sp)
  keggSP <- as.data.frame(keggSP)
  gSetsSP <- lapply(rownames(keggSP), function(X){
    X <- try(keggLink(X)[,2], silent = TRUE)
    if(class(X)!="try-error"){
      X <- X[grepl(paste0(sp,":"),X)]
      return(gsub(paste0(sp,":"),"",X) )
    } else{
      return(NA)
    }
  })
  names(gSetsSP) <- as.character(keggSP[,1])
  return(gSetsSP)
}
