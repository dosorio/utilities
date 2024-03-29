networkScoring <- function(rNetwork, geneSets){
  
  # A specificity value is computed for each gene depending on how many gene sets it is included. 
  # Retrieves a value between 0 and 1
  geneSpecificty <- 1-(table(unlist(geneSets))/length(geneSets))
  geneSpecificty <- geneSpecificty/max(geneSpecificty)
  
  # The names of the genes included in the network are collected
  geneNames <- names(V(rNetwork))
  
  # The score is computed using the averages shorted paths between vertices multiplied by the 
  # specificity and the proportion observed of each gene set. 
  geneSetScores <- lapply(geneSets, function(gSet){
    setSize <- length(gSet)
    gSet <- gSet[gSet %in% geneNames]
    obsProportion <- length(gSet)/setSize
    setSpecificity <- mean(geneSpecificty[gSet])
    setDistances <- induced_subgraph(rNetwork, gSet)
    setDistances <- distances(setDistances)
    setDistances[!is.finite(setDistances)] <- 0
    gmSetDistances <- exp(mean(log(setDistances[setDistances > 0])))
    setScore <- (gmSetDistances * obsProportion) / setSpecificity
    return(setScore)
  })
  
  # Scores are normalized using the BOX-COX normalization method
  geneSetScores <- unlist(geneSetScores)
  geneSetScores[!is.finite(geneSetScores)] <- 0
  
  BC <- MASS::boxcox(geneSetScores[geneSetScores >  0] ~ 1, plotit = FALSE)
  BC <- 1 + abs(BC$x[which.max(BC$y)])
  geneSetScores <- geneSetScores ^ BC
  
   # Scores are standardized
  geneSetScores <- scale(geneSetScores)[,1]
  
  # Return
  return(geneSetScores)
}
