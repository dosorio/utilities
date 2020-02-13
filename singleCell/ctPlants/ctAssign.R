library(Seurat)
PTI <- readRDS('PTI.combined.rds')
PTI <- PTI[,1:1000]
UMAPPlot(PTI)

ctAssign <- function(X){
  require(pbapply)
  
  # Data
  M <- read.csv('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/ctPlants/M.csv')
  S <- read.csv('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/ctPlants/S.csv')
  
  # Matching
  tM <- M[M$Locus %in% rownames(X),]
  tS <- S[S$Locus %in% rownames(X),]
  
  # Filtering
  xM <- X[tM$Locus,]
  xS <- X[tS$Locus,]
  
  # Function
  outValues <- pbsapply(seq_len(ncol(xM)), function(cell){
    outValues <- c()
    for(cT in 2:ncol(tM)){
      outValues[(cT-1)] <- mean(xS[,cell] * tS[,cT], na.rm = TRUE) * mean(xM[tM[,cT] >0 ,cell] >0, na.rm = TRUE)
    }
    outValues <- outValues/sum(outValues,na.rm = TRUE)
    return(outValues)
  })
  outValues <- t(outValues)
  tissueNames <- colnames(tM)[2:ncol(tM)]
  colnames(outValues) <- tissueNames
  
  out <- list()
  out$Values <- outValues
  cType <- apply(outValues,1,function(X){tissueNames[which.max(X)]})
  cType[lengths(cType) < 1] <- 'NA'
  out$cType <- cType
  return(out)
}

cT <- ctAssign(PTI@assays$RNA@data)
Idents(PTI) <- cT[[2]]
UMAPPlot(PTI)
