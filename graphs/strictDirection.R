strictDirection <- function(X, lambda = 1){
  S <- X
  S[abs(S) < abs(t(S))] <- 0
  O <- (((1-lambda) * X) + (lambda * S))
  return(O)
}
