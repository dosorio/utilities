bootstrapNet <- function(X, B = 100, nCell = 1000, nCom = 3) {
  makeGraph <- function(X){
    diag(X) <- NA
    qT <- quantile(abs(X), 0.9, na.rm = TRUE)
    X[abs(X) < qT] <- NA
    X <- reshape2::melt(X)
    X <- X[complete.cases(X), ]
    X <- igraph::graph_from_data_frame(X)
    return(X)
  }
  makeBnet <- function(X, nCell = 1000, nCom = 3){
    cellNames <- colnames(X)
    bX <- X[, sample(cellNames, nCell, replace = TRUE)]
    bX <- bX[rowMeans(bX != 0) > 0.05, ]
    oN <- pcNet(bX, nCom = nCom)
    oN <- makeGraph(oN)
    return(oN)
  }
  oN <- pcNet(X, nCom)
  oN <- makeGraph(oN)
  sapply(seq_len((B)), function(b){
    nN <- try(makeBnet(X, nCell, nCom), silent = TRUE)
    if (class(nN) == 'igraph') {
      oN <<- igraph::intersection(oN, nN, keep.all.vertices = FALSE)
    }
  })
  oN <- oN[,]
  oN[oN != 0] <- 1
  oN <- oN[apply(oN,1,sum) > 0, apply(oN,2,sum)]
  oN <- reshape2::melt(as.matrix(oN))
  oN <- oN[oN[,3] != 0,]
  oN <- igraph::graph_from_data_frame(oN)
  return(oN)
}
