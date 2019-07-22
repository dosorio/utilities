bootstrapNet <- function(X, B = 100, qThreshold = 0.9, nCell = 1000, nCom = 3) {
  source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/makeBnet.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/makeGraph.R")
  oN <- pcNet(X, nCom)
  oN <- makeGraph(oN, qThreshold = qThreshold)
  b <- 0
  while (b < B) {
    nN <- try(makeBnet(X, qThreshold = qThreshold, nCell, nCom), silent = TRUE)
    if (class(nN) == 'igraph') {
      oN <- igraph::intersection(oN, nN, keep.all.vertices = FALSE)
      b <- b + 1
    }
  }
  oN <- oN[,]
  oN <- oN[apply(oN,1,sum) > 0, apply(oN,2,sum)> 0]
  oN <- reshape2::melt(as.matrix(oN))
  oN <- oN[oN[,3] != 0,]
  oN <- igraph::graph_from_data_frame(oN)
  return(oN)
}
