bootstrapNet <- function(X, B = 100, nCell = 1000, nCom = 3) {
  makeBnet <- function(X, nCell = 1000, nCom = 3){
    cellNames <- colnames(X)
    bX <- X[, sample(cellNames, nCell, replace = TRUE)]
    bX <- bX[rowMeans(bX != 0) > 0.05, ]
    oN <- pcNet(bX, nCom = nCom)
    diag(oN) <- NA
    qT <- quantile(abs(oN), 0.9, na.rm = TRUE)
    oN[abs(oN) < qT] <- NA
    oN <- reshape2::melt(oN)
    oN <- oN[complete.cases(oN), ]
    oN <- igraph::graph_from_data_frame(oN)
    return(oN)
  }
  oN <- makeBnet(X, nCell, nCom)
  sapply(seq_len((B-1)), function(b){
    nN <- try(makeBnet(X, nCell, nCom), silent = TRUE)
    if (class(nN) == 'igraph') {
      oN <<- igraph::intersection(oN, nN, keep.all.vertices = FALSE)
    }
  })
  oN <- oN[,]
  oN <- oN[apply(oN,1,sum) > 0, apply(oN,2,sum)]
  oN <- reshape2::melt(as.matrix(oN))
  oN <- oN[oN[,3] != 0,]
  oN <- igraph::graph_from_data_frame(oN)
  return(oN)
}
