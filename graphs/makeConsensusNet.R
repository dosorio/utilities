makeConsensusNet <- function(fileList, threshold=0.9){
  require(igraph)
  makeNet <- function(X, threshold){
    rNet <- read.csv(X, header = TRUE, row.names = 1)
    diag(rNet) <- NA
    thresholdValue <- quantile(abs(rNet), threshold, na.rm = TRUE)
    rNet[isTRUE(rNet < thresholdValue)] <- NA
    rNet[upper.tri(rNet, diag = TRUE)] <- NA
    rNet <- reshape2::melt(as.matrix(rNet))
    rNet <- rNet[complete.cases(rNet),]
    return(rNet)
  }
  message(1)
  oNet <- makeNet(fileList[1])
  oNet <- graph_from_data_frame(oNet, directed = FALSE)
  for(i in seq_along(fileList)[-1]){
    message(i)
    nNet <- makeNet(fileList[i])
    nNet <- graph_from_data_frame(nNet, directed = FALSE)
    oNet <- intersection(oNet, nNet, keep.all.vertices = FALSE)
  }
  return(oNet)
}
