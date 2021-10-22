df2matrix <- function(df, columns, rows, values){
  require(Matrix)
  df <- df[,c(columns, rows, values)]
  sList <- unique(df[,columns])
  fList <- unique(df[,rows])
  oMatrix <- Matrix(matrix(0, nrow = length(fList), ncol = length(sList)))
  colnames(oMatrix) <- sList
  rownames(oMatrix) <- fList
  O <- pbapply::pbapply(df,1,function(X){
    oMatrix[X[2],X[1]] <<- as.numeric(X[3])
  })
  return(oMatrix)
}
