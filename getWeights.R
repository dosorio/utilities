getLeverage <- function(X){
  diag(X %*% solve(t(X) %*% X) %*% t(X))
}
