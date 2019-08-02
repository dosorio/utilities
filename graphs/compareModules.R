#  cX <- cluster_walktrap(netX, steps = 100)
#  cY <- cluster_walktrap(netY, steps = 100)

compareModules <- function(cX, cY){
  S <- sapply(unique(cX$membership), function(mX){
    sapply(unique(cY$membership), function(mY){
      gX <- cX$names[cX$membership == mX]
      gY <- cY$names[cY$membership == mY]
      allG <- unique(c(gX,gY))
      sMatrix <- matrix(FALSE, nrow = length(allG), ncol = 2)
      rownames(sMatrix) <- allG
      sMatrix[gX,1] <- TRUE
      sMatrix[gY,2] <- TRUE
      round(mean(apply(sMatrix,1,all)),2)
    })
  })
  return(S)
}
