list2matrix <- function(X){
  fList <- unique(unlist(X))
  X <- sapply(X, function(i){fList %in% i})
  rownames(X) <- fList
  return(X)
}
