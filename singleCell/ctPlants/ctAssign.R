ctAssign <- function(X){
  require(pbapply)
  
  # Data
  M <- read.csv('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/ctPlants/M.csv')
  S <- read.csv('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/ctPlants/S.csv')
  
  # Matching
  tM <- M[M$Locus %in% rownames(X),]
  tS <- S[S$Locus %in% rownames(X),]
  
  # Filtering
  xM <- X[M$Locus,]
  xS <- X[S$Locus,]
  
  # Function
  pbsapply(seq_len(ncol(xM)), function(cell){
    outValues <- c()
    for(cT in 2:ncol(tM)){
      outValues[(cT-1)] <- mean(xS * tS[,cT]) * mean(xM[tM[,cT] >0 ,] >0)
    }
    outValues <- outValues/sum(outValues)
    return(outValues)
  })
}