  makeBnet <- function(X, qThreshold = 0.9, nCell = 1000, nCom = 3){
    source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/pcNet.R")
    cellNames <- colnames(X)
    bX <- X[, sample(cellNames, nCell, replace = TRUE)]
    bX <- bX[rowMeans(bX != 0) > 0, ]
    oN <- pcNet(bX, nCom = nCom)
    oN <- makeGraph(oN, qThreshold = qThreshold)
    return(oN)
  }
