makeWindows <- function(X, W=100, O=50){
  nCol <- ncol(X)
  zoo::rollapply(seq_len(nCol), W, by = O, c)
}
