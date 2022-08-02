qLM <- function(X, y){
  X <- cbind(1,X)
  m <- solve(t(X) %*% X) %*% t(X) %*% y
  yHat <- X %*% m
  MSM <- rowSums((t(yHat) - colMeans(y))^2)/(p-1)
  MSE <- colSums((y - yHat)^2)/(n-p)
  P <- pf(MSM/MSE, df1 = p-1, df2 = n-p, lower.tail = FALSE)
  O <- data.frame(beta = t(m)[,1],P)
  return(O)
}
