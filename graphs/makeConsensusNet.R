makeConsensusNet <- function(fileList, qThreshold=0.9){
  require(igraph)
  makeNet <- function(X, qThreshold = qThreshold){
    rNet <- read.csv(X, header = TRUE, row.names = 1)
    diag(rNet) <- NA
    thresholdValue <- quantile(abs(rNet), qThreshold, na.rm = TRUE)
    rNet[rNet < thresholdValue] <- NA
    rNet[upper.tri(rNet, diag = TRUE)] <- NA
    rNet <- reshape2::melt(as.matrix(rNet))
    rNet <- rNet[complete.cases(rNet),]
    return(rNet)
  }
  message(1)
  oNet <- makeNet(fileList[1], qThreshold = qThreshold)
  oNet <- graph_from_data_frame(oNet, directed = FALSE)
  for(i in seq_along(fileList)[-1]){
    message(i)
    nNet <- makeNet(fileList[i], qThreshold = qThreshold)
    nNet <- graph_from_data_frame(nNet, directed = FALSE)
    oNet <- intersection(oNet, nNet, keep.all.vertices = FALSE)
  }
  oNet <- oNet[,]
  oNet <- oNet[apply(oNet, 1, sum) > 0, apply(oNet, 2, sum) > 0]
  oNet <- graph_from_adjacency_matrix(oNet)
  return(oNet)
}
