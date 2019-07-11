makeConsensusNet <- function(fileList){
  oNet <- read.table(fileList[1], header = TRUE)
  oNet <- graph_from_data_frame(oNet, directed = FALSE)
  for(i in seq_along(fileList)[-1]){
    nNet <- read.table(fileList[i], header = TRUE)
    nNet <- graph_from_data_frame(nNet, directed = FALSE)
    oNet <- intersection(oNet, nNet, keep.all.vertices = FALSE)
  }
  return(oNet)
}
