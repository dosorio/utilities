makeConsensusNet <- function(fileList){
  makeNet <- function(X){
    rNet <- read.csv(X, header = TRUE, row.names = 1)
    diag(rNet) <- NA
    thresholdValue <- quantile(abs(rNet), 0.9, na.rm = TRUE)
    rNet[isTRUE(rNet < thresholdValue)] <- NA
    rNet[upper.tri(rNet, diag = TRUE)] <- NA
    rNet <- reshape2::melt(as.matrix(rNet))
    rNet <- rNet[complete.cases(rNet),]
    return(rNet)
  }
  oNet <- makeNet(fileList[1])
  oNet <- graph_from_data_frame(oNet, directed = FALSE)
  for(i in seq_along(fileList)[-1]){
    nNet <- makeNet(fileList[i])
    nNet <- graph_from_data_frame(nNet, directed = FALSE)
    oNet <- intersection(oNet, nNet, keep.all.vertices = FALSE)
  }
  return(oNet)
}
