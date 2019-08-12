makeWindows <- function(X, W=10, O=50){
  nCol <- ncol(X)
  zoo::rollapply(seq_len(nCol), ceiling(nCol/W), by = O, c)
}
