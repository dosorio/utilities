ctAssign <- function(X){
  require(pbapply)
  require(Matrix)
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
  outValues <- pbsapply(seq_len(ncol(tM))[-1], function(tissue){
    outValues <- (colMeans(t(t(xS)*tS[,tissue])) * colMeans(xM[tM[,tissue] > 0,] > 0))
    return(outValues)
  })
  tissueNames <- colnames(tM)[2:ncol(tM)]
  colnames(outValues) <- tissueNames
  
  outValues <- outValues/rowSums(outValues)
  plot(outValues[2,])
  out <- list()
  out$Values <- outValues
  cType <- apply(outValues,1,function(X){tissueNames[which.max(X)]})
  cType[lengths(cType) < 1] <- 'NA'
  out$cType <- unlist(cType)
  return(out)
}

# library(Seurat)
# 
# PTI <- readRDS('PTI.combined.rds')
# #PTI <- NormalizeData(PTI, normalization.method = 'RC', scale.factor = 1e6)
# PTI <- RunICA(PTI)
# PTI <- RunUMAP(PTI, dims = 1:50, reduction = 'ica')
# PTI <- RunTSNE(PTI, reduction = 'ica', perplexity = 1000)
# TSNEPlot(PTI)
# UMAPPlot(PTI)
# 
# cT <- ctAssign(PTI@assays$RNA@data)
# cT$cType[cT$cType == 'NA'] <- NA
# Idents(PTI) <- cT[[2]]
# UMAPPlot(PTI)
# TSNEPlot(PTI)
# 
# PTI <- RunTSNE(PTI, dims = 1:20, perplexity =1500)
# TSNEPlot(PTI)
# 
# library(phateR)
# O <- phate(t(as.matrix(PTI@assays$RNA@data)))
# plot(O$embedding, col = as.factor(cT$cType))
