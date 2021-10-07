list2matrix <- function(X){
  fList <- unique(unlist(X))
  X <- sapply(CMAP, function(i){fList %in% i})
  rownames(X) <- fList
  return(X)
}
