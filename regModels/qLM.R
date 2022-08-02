qLM <- function(X, y){
  X <- cbind(1,X)
  n <- nrow(X)
  p <- ncol(X)
  m <- solve(t(X) %*% X) %*% t(X) %*% t(y)
  yHat <- X %*% m
  y <- t(y)
  MSM <- rowSums((t(yHat) - colMeans(y))^2)/(p-1)
  MSE <- colSums((y - yHat)^2)/(n-p)
  P <- pf(MSM/MSE, df1 = p-1, df2 = n-p, lower.tail = FALSE)
  O <- data.frame(beta = t(m)[,2],P)
  return(O)
}
