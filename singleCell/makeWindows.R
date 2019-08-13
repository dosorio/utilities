makeWindows <- function(cMatrix, nWindows, ovlpProportion = 0.5){
  nCol <- ncol(cMatrix)
  makeM <- function(X){
    O <- (X * ovlpProportion)
    M <- zoo::rollapply(seq_len(nCol), X, by=O, c)
    abs(dim(M)[1] - nWindows)
  }
  Z <- optimize(f = makeM, interval = seq_len(nCol))
  zoo::rollapply(seq_len(nCol), Z$minimum, by=Z$minimum*ovlpProportion, c)
}
