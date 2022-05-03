 matrixFT <- function(x,y, p.adj = TRUE){
  rs_x = rowSums(x)
  rs_y = rowSums(y)
  RESULT = tcrossprod(x, y)
  RESULT = as(RESULT, 'dgTMatrix')
  RESULT <- reshape2::melt(as.matrix(RESULT))
  b = rowSums(y)[RESULT[,2]] - RESULT[,3]
  c = rowSums(x)[RESULT[,1]] - RESULT[,3]
  d = RESULT[,3]
  a = tcrossprod(x!=1, y!=1)
  a = c(as.matrix(a))
  o <- HighSpeedStats::ultrafastfet(a, b, c, d)
  if(p.adj){
    o <- matrix(p.adjust(o), nrow = nrow(x), ncol = nrow(y), byrow = FALSE)  
  } else {
    o <- matrix(o, nrow = nrow(x), ncol = nrow(y), byrow = FALSE)
  }
  return(o)
}
