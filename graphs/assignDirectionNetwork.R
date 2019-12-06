assignDirectionNetwork <- function(igraphNetwork, countMatrix, bootR= 10){
   bnModel <- boot.strength(countMatrix, algorithm = 'hc', R = bootR)
   bnModel$direction <- round(bnModel$direction,1)
   bnModel <- bnModel[bnModel$direction >= 0.5,]
   bY <- graph_from_data_frame(bnModel)
   bY <- as.matrix(bY[])
   igraphNetwork <- as.matrix(igraphNetwork[])
   bY <- bY[rownames(igraphNetwork),colnames(igraphNetwork)]
   igraphNetwork[!(igraphNetwork != 0 & bY != 0)] <- 0
   igraphNetwork <- (graph_from_adjacency_matrix(igraphNetwork, weighted = TRUE))
   return(igraphNetwork)
}
