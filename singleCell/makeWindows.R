# W is the number of desired Windows
makeWindows <- function(X, W=10){
  nCol <- ncol(X)
  W <- W+1
  O <- nCol/W
  zoo::rollapply(X, W, by=O, c)
}
