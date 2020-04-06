sccTenifoldKNK <- function(X, Y, qc_mtThreshold = 0.1, qc_minLSize = 1000, qc_minNvalues = 100, nc_q = 0.95, nc_nCell = 500, nc_nNet = 25, nc_K = 2, nc_denoiseNet = TRUE, ma_nDim = 30){
  scQC <- function(X, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize){
    if(class(X) == 'Seurat'){
      countMatrix <- X@assays$RNA@counts
    } else {
      countMatrix <- X
    }
    librarySize <- Matrix::colSums(countMatrix)
    countMatrix <- countMatrix[,librarySize >= minLSize]
    librarySize <- Matrix::colSums(countMatrix)
    mtGenes <- grep('^MT-',toupper(rownames(countMatrix)), ignore.case = TRUE)
    nGenes <- Matrix::colSums(countMatrix != 0)
    
    genesLM <- lm(nGenes~librarySize)
    genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
    
    if(isTRUE(length(mtGenes) > 0)){
      mtCounts <- Matrix::colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
      mtProportion <- mtCounts/librarySize
      mtLM <- lm(mtCounts~librarySize)
      mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
      selectedCells <- mtCounts > mtLM$lwr & mtCounts < mtLM$upr & nGenes > genesLM$lwr & nGenes < genesLM$upr & mtProportion <= mtThreshold & librarySize < 2 * mean(librarySize)
    } else {
      selectedCells <- nGenes > genesLM$lwr & nGenes < genesLM$upr & librarySize < 2 * mean(librarySize)
    }
    if(class(X) == 'Seurat'){
      selectedCells <- colnames(countMatrix)[selectedCells]
      X <- subset(X, cells = selectedCells)
    } else {
      selectedCells <- which(selectedCells)
      X <- countMatrix[,selectedCells]
    }
    return(X)
  }
  geneFilter <- function(X, minNvalues = qc_minNvalues){
    X[Matrix::rowSums(X!=0) >= minNvalues,]
  }
  cpmNormalization <- function(X){
    Matrix::t(Matrix::t(X)/Matrix::colSums(X)) * 1e6
  }
  sccNet <- function(X, q = nc_q, nCell = nc_nCell, nNet = nc_nNet, K = nc_K, denoiseNet = nc_denoiseNet){
    nGenes <- nrow(X)
    gList <- rownames(X)
    set.seed(1)
    oNet <- pbapply::pbsapply(seq_len(nNet), function(Z){
      tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
      tNet <- cor(as.matrix(tNet), method = 'sp')
      while(any(is.na(tNet))){
        tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
        tNet <- cor(as.matrix(tNet), method = 'sp')
      }
      tNet <- round(tNet,1)
      diag(tNet) <- 0
      tNet[abs(tNet) < quantile(abs(tNet), q, na.rm = TRUE)] <- 0
      tNet <- Matrix::Matrix(tNet)
      return(tNet)  
    })
    
    aNet <- matrix(0, nGenes, nGenes)
    rownames(aNet) <- colnames(aNet) <- gList
    for(i in seq_along(oNet)){
      tNet <- oNet[[i]]
      if(denoiseNet){
        if(K == ncol(X)){
          tNet <- svd(tNet)
          tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)  
          rownames(tNet) <- colnames(tNet) <- gList  
        } else if(K == 1){
          tNet <- RSpectra::svds(tNet,K)
          tNet <- tNet$u %*% (tNet$d) %*% t(tNet$v)  
          rownames(tNet) <- colnames(tNet) <- gList  
        } else {
          tNet <- RSpectra::svds(tNet,K, maxitr = 1e6)
          tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)  
          rownames(tNet) <- colnames(tNet) <- gList  
        }
        
      }
      tNet <- Matrix::Matrix(tNet)
      aNet <- aNet + tNet[gList, gList]
    }
    aNet <- aNet/length(oNet)
    aNet <- Matrix::Matrix(aNet)
    return(aNet)
  }
  manifoldAlignment <- function(X, Y, d = ma_nDim){
    sharedGenes <- intersect(rownames(X), rownames(Y))
    X <- X[sharedGenes, sharedGenes]
    Y <- Y[sharedGenes, sharedGenes]
    L <- diag(length(sharedGenes))
    wX <- X+1
    wY <- Y+1
    mu <- 0.9
    eps <- 1e-8
    wXY <- L
    wXY <- mu * (sum(wX) + sum(wY)) / (2 * sum(wXY)) * wXY
    W <- rbind(cbind(wX, wXY), cbind(Matrix::t(wXY), wY))
    nNodes <- nrow(W)
    diag(W) <- 0
    diag(W) <- Matrix::rowSums(W)
    #E <- suppressWarnings(RSpectra::eigs(W, ncol(W)))
    E <- suppressWarnings(RSpectra::eigs(W, d, 'SM'))
    newOrder <- order(E$values)
    eVal <- suppressWarnings(as.numeric(E$values[newOrder]))
    eVec <- E$vectors[,newOrder]
    eVec <- suppressWarnings(apply(eVec,2,as.numeric))
    selectedV <- (eVal > eps)
    eVal <- eVal[selectedV]
    eVec <- eVec[,selectedV]
    cNorms <- apply(eVec,2,function(X){norm(as.matrix(X))})
    eVec <- t(t(eVec)/cNorms)
    eVec[((nNodes/2)+1):nNodes,] <- -1 * eVec[((nNodes/2)+1):nNodes,]
    alignedNet <- eVec[,seq_len(d)]
    colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
    rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
    return(alignedNet)
  }
  dRegulation <- function(manifoldOutput){
    
    geneList <- rownames(manifoldOutput)
    geneList <- geneList[grepl('^X_', geneList)]
    geneList <- gsub('^X_','', geneList)
    nGenes <- length(geneList)
    
    eGenes <- nrow(manifoldOutput)/2
    
    eGeneList <- rownames(manifoldOutput)
    eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
    eGeneList <- gsub('^Y_','', eGeneList)
    
    if(nGenes != eGenes){
      stop('Number of identified and expected genes are not the same')
    }
    if(!all(eGeneList == geneList)){
      stop('Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order')
    }
    
    dMetric <- sapply(seq_len(nGenes), function(G){
      X <- manifoldOutput[G,]
      Y <- manifoldOutput[(G+nGenes),]
      I <- rbind(X,Y)
      O <- dist(I)
      O <- as.numeric(O)
      return(O)
    })
    #dMetric[dMetric < 1e-10] <- 1e-10
    
    ### BOX-COX
    lambdaValues <- seq(-2,2,length.out = 1000)
    lambdaValues <- lambdaValues[lambdaValues != 0]
    BC <- MASS::boxcox(dMetric~1, plot=FALSE, lambda = lambdaValues)
    BC <- BC$x[which.max(BC$y)]
    if(BC < 0){
      nD <- 1/(dMetric ^ BC)
    } else {
      nD <- dMetric ^ BC
    }
    Z <- scale(nD)
    E <- mean(dMetric^2)
    FC <- dMetric^2/E
    pValues <- pchisq(q = FC,df = 1,lower.tail = FALSE)
    pAdjusted <- p.adjust(pValues, method = 'fdr')
    dOut <- data.frame(
      gene = geneList, 
      distance = dMetric,
      Z = Z,
      FC = FC,
      p.value = pValues,
      p.adj = pAdjusted
    )
    dOut <- dOut[order(dOut$p.value),]
    dOut <- as.data.frame.array(dOut)
    return(dOut)
  }
  
  X <- scQC(X, minLSize = qc_minLSize, mtThreshold = qc_mtThreshold)
  X <- geneFilter(X, minNvalues = qc_minNvalues)
  X <- cpmNormalization(X)
  X <- sccNet(X, q = nc_q, nCell = nc_nCell, nNet = nc_nNet, K = nc_K, denoiseNet = nc_denoiseNet)
  X <- as.matrix(X)
  diag(X) <- 0
  X[abs(X) < quantile(abs(X), nc_q)] <- 0
  Y <- X
  Y[,gKO] <- 0
  X <- Matrix::Matrix(X)
  Y <- Matrix::Matrix(Y)
  mA <- manifoldAlignment(X, Y, d = ma_nDim)
  dR <- dRegulation(mA)
  outputResult <- list()
  outputResult$WT <- X
  outputResult$KO <- Y
  outputResult$manifoldAlignment <- mA
  outputResult$diffRegulation <- dR
  return(outputResult)
}
