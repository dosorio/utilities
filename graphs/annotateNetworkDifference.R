annotateNetworkDifference <- function(X, Y, type="KEGG"){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/mouse.R")
  if(type == "KEGG"){
    D <- difference(X,Y)
    louvainClusters <- cluster_louvain(D)
    annotationKEGG <- lapply(unique(louvainClusters$membership), function(C){
      out <- try(makeMouseEnrichmentKEGG(louvainClusters$names[louvainClusters$membership %in% C]), silent = TRUE)
      if(class(out) != "try-error"){return(out)}
    })
    annotationKEGG <- do.call(rbind.data.frame, annotationKEGG)
    annotationKEGG <- annotationKEGG[order(annotationKEGG$p.adjust),]
    rownames(annotationKEGG) <- NULL
    return(annotationKEGG)
  }
  if(type == "GO"){
    D <- difference(X,Y)
    louvainClusters <- cluster_louvain(D)
    annotationGO <- lapply(unique(louvainClusters$membership), function(C){
      out <- try(makeMouseEnrichmentGO(louvainClusters$names[louvainClusters$membership %in% C]), silent = TRUE)
      if(class(out) != "try-error"){return(out)}
    })
    annotationGO <- do.call(rbind.data.frame, annotationGO)
    annotationGO <- annotationGO[order(annotationGO$p.adjust),]
    rownames(annotationGO) <- NULL
    return(annotationGO)
  }
}
