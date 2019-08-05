#  cX <- cluster_walktrap(netX, steps = 100)
#  cY <- cluster_walktrap(netY, steps = 100)

compareModules <- function(cX, cY){
  cX <- groups(cX)
  cY <- groups(cY)
  S <- sapply(cX, function(gX){
    sapply(cY, function(gY){
      allG <- unique(c(gX,gY))
      sMatrix <- matrix(FALSE, nrow = length(allG), ncol = 2)
      rownames(sMatrix) <- allG
      sMatrix[gX,1] <- TRUE
      sMatrix[gY,2] <- TRUE
      round(mean(apply(sMatrix,1,all)),2)
    })
  })
  S <- S[apply(S,1,function(X){any(X > 0)}),apply(S,2,function(X){any(X > 0)})] 
  return(t(S))
}
