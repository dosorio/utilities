 <- function(a,b, df=0){
  if(nrow(a) == 1){
    a <- rbind(0,a)
    b <- rbind(0,b)
  }
  aa <- colSums(a * a)
  bb <- colSums(b * b)
  ab <- t(a) %*% b
  D <- sqrt(matrix(data = rep(aa, each=ncol(b)), ncol= ncol(b), byrow = TRUE) + matrix(data = rep(bb, each=ncol(a)), ncol= ncol(a), byrow = FALSE) - 2*ab)
  if(df == 1){
    D <- D * (1-diag(length(diag(D))))
  }
  return(D)
}
