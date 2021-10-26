collapseEnrichrPathways <- function(geneList, dbList, jaccardThreshold = NA, nPathways = NA){
  require(enrichR)
  E <- enrichr(genes = geneList, databases = dbList)
  E <- do.call(rbind.data.frame, E)
  E <- E[E$Adjusted.P.value < 0.05,]
  E <- E[order(E$Adjusted.P.value, decreasing = FALSE),]
  G <- strsplit(E$Genes, split = ';')
  names(G) <- E$Term
  G <- UpSetR::fromList(G)
  if(!is.na(jaccardThreshold)){
    E$G <- cutree(hclust(proxy::dist(G, by_rows = FALSE, method = "Jaccard")), h = jaccardThreshold)  
  }
  if(!is.na(nPathways)){
    E$G <- cutree(hclust(proxy::dist(G, by_rows = FALSE, method = "Jaccard")), k = nPathways)
  }
  mainPathways <- lapply(split(E, E$G), function(X){X[1,]})
  mainPathways <- do.call(rbind.data.frame, mainPathways)
  mainPathways$P <- lapply(split(E, E$G), function(X){paste0(unique(X$Term[-1]), collapse = ';')})
  mainPathways <- mainPathways[,c('Term', 'Overlap', 'Adjusted.P.value', 'Genes', 'P')]
  return(mainPathways)
}
