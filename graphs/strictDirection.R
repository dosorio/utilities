strictDirection <- function(X){
  lT <- X[lower.tri(X)]
  uT <- t(X)[lower.tri(X)]
  all(lT == uT)
  lT[abs(lT) < abs(uT)] <- 0
  uT[abs(uT) < abs(lT)] <- 0
  X[lower.tri(X)] <- lT
  X[upper.tri(X)] <- uT
  return(X)
}
