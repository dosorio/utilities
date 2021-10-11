list2matrix <- function(X){
  fList <- unique(unlist(X))
  O <- sapply(X, function(i){fList %in% i})
  rownames(O) <- fList
  return(O)
}
