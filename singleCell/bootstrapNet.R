bootstrapNet <- function(X, B = 100, nCell = 1000, nCom = 3, nCores = 1) {
  makeBnet <- function(X, nCell = 1000, nCom = 3, nCores = 1) {
    cellNames <- colnames(X)
    bX <- X[, sample(cellNames, nCell)]
    bX <- bX[rowMeans(bX != 0) > 0.05, ]
    oN <- pcNet(bX, nCom = nCom,  nCores = nCores)
    diag(oN) <- NA
    qT <- quantile(abs(oN), 0.9, na.rm = TRUE)
    oN[abs(oN) < qT] <- NA
    oN <- reshape2::melt(oN)
    oN <- oN[complete.cases(oN), ]
    oN <- igraph::graph_from_data_frame(oN)
    return(oN)
  }
  oN <- makeBnet(X, nCell, nCom, nCores)
  b <- 1
  while (b <= B) {
    #message(b)
    nN <- try(makeBnet(X, nCell, nCom, nCores), silent = TRUE)
    if (class(nN) == 'igraph') {
      oN <- igraph::intersection(oN, nN, keep.all.vertices = FALSE)
      b <- b + 1
    }
  }
  oN <- oN[,]
  oN <- oN[apply(oN,1,sum) > 0, apply(oN,2,sum)]
  oN <- reshape2::melt(as.matrix(oN))
  oN <- oN[oN[,3] != 0,]
  oN <- igraph::graph_from_data_frame(oN)
  return(oN)
}
