bootstrapNet <- function(X, B = 100, qThreshold = 0.9, nCell = 1000, nCom = 3) {
  source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/makeBnet.R")
  source("https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/makeGraph.R")
  oN <- pcNet(X, nCom)
  oN <- makeGraph(oN, qThreshold = qThreshold)
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
