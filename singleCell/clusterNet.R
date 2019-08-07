clusterNet <- function(X){
  source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/getCellClusters.R")
  clustersID <- getCellClusters(X)
  clustersCells <- lapply(unique(clustersID), function(C){
    cCells <- X[,clustersID == C]
    cCells <- cCells[apply(cCells,1,sum) > 0,]
    return(cCells)
  })
  clustersNet <- pblapply(clustersCells, pcNet)
  clustersNet <- pblapply(clustersNet, makeGraph)
  consensusNet <- clustersNet[[1]]
  genesChange <- unlist(lapply(clustersNet[-1], function(Z){
    cNet_sDF <<- cleanNet(intersection(consensusNet, Z, keep.all.vertices = FALSE))
    length(V(cNet_sDF))
  }))
  return(consensusNet)
}
