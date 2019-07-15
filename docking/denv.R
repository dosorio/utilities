args = commandArgs(trailingOnly=TRUE)
proteins <- list.files(args[1], full.names = "TRUE")
ligands <- list.files(args[2], full.names = "TRUE")

nRun <- 10

# proteins <- list.files("proteins/", full.names = TRUE)
# ligands <- list.files("ligands/", full.names = TRUE)

findCenter <- function (X) {
  fileContent <- readLines(X)
  fileContent <- fileContent[grepl("^ATOM", fileContent)]
  fileContent <- strsplit(fileContent, "[[:space:]]+")
  fileContent <- do.call(rbind.data.frame, fileContent)
  fileContent <- fileContent[, 7:9]
  fileContent <- t(apply(fileContent, 1, as.numeric))
  colMeans(fileContent)
}

findSize <- function (X) {
  fileContent <- readLines(X)
  fileContent <- fileContent[grepl("^ATOM", fileContent)]
  fileContent <- strsplit(fileContent, "[[:space:]]+")
  fileContent <- do.call(rbind.data.frame, fileContent)
  fileContent <- fileContent[, 7:9]
  fileContent <- t(apply(fileContent, 1, as.numeric))
  fileContent <- abs(t(t(fileContent) - colMeans(fileContent)))
  apply(fileContent, 2, max)
}

sapply(seq_len(nRun), function(run){
  dir.create(paste0("out/r", run))
  sapply(proteins, function(X){
    sapply(ligands, function(Y){
      CM <- round(findCenter(X),1)
      D <- ceiling(findSize(X))+5
      outFile <- paste0(getwd(),"/out/r",run,"/",gsub(".pdbqt$", "", basename(X)),"_",gsub(".pdbqt$", "", basename(Y)),".pdbqt")
      if(!file.exists(outFile)){
        configFile <- c(
          paste0("receptor = ", X),
          paste0("ligand = ", Y),
          paste0("center_x = ", CM[1]),
          paste0("center_y = ", CM[2]),
          paste0("center_z = ", CM[3]),
          paste0("size_x = ", D[1]),
          paste0("size_y = ", D[2]),
          paste0("size_z = ", D[3]),
          paste0("exhaustiveness =", 100),
          paste0("out =", outFile)
        )
        writeLines(configFile, "vinaDocking.conf")
        system("vina --config vinaDocking.conf")
        unlink("vinaDocking.conf") 
      }
    })
  })
})
outP <- sapply(proteins, function(X){
  outL <- sapply(ligands, function(Y){
    outFiles <- paste0(getwd(),"/out/r",seq_len(nRun),"/",gsub(".pdbqt$", "", basename(X)),"_",gsub(".pdbqt$", "", basename(Y)),".pdbqt")
    fileContent <- lapply(outFiles, readLines)
    fileContent <- lapply(fileContent, function(X){
      X[seq_len(grep("^ENDMDL", X)[1])]
    })
    mEnergy <- unlist(lapply(fileContent, function(X){
      X <- X[grepl("RESULT", X)]
      as.numeric(unlist(strsplit(X, "[[:space:]]+"))[4])
    }))
    cmLigand <- lapply(fileContent, function(X){
      atoms <- (grepl("^ATOM",X) | grepl("^HETATM", X))
      cmLigand <- X[atoms]
      cmLigand <- strsplit(cmLigand, "[[:space:]]+")
      cmLigand <- do.call(rbind.data.frame, cmLigand)
      cmLigand <- cmLigand[,6:8]
      cmLigand <- t(apply(cmLigand, 1, as.numeric))
      cmLigand <- colMeans(cmLigand)
    })
    cmLigand <- do.call(rbind.data.frame, cmLigand)
    colnames(cmLigand) <- c("X", "Y", "Z")
    dLigand <- apply(cmLigand,2,sd)
    mLigand <- apply(cmLigand,2,mean)
    dEnergy <- sd(mEnergy)
    mEnergy <- mean(mEnergy)
    X <- gsub(".pdbqt", "", basename(X))
    Y <- gsub(".pdbqt", "", basename(Y))
    outValues <- data.frame(PROTEIN = X, 
                            LIGAND= Y, 
                            X = round(mLigand[1],2), 
                            sdX = round(dLigand[1],2),
                            Y = round(mLigand[2],2), 
                            sdY = round(dLigand[2],2),
                            Z = round(mLigand[3],2), 
                            sdZ = round(dLigand[3],2),
                            ENERGY = round(mEnergy,2),
                            sdENERGY = round(dEnergy,2))
    return(outValues)
  }, simplify = FALSE)
  outL <- do.call(rbind.data.frame, outL)
  return(outL)
}, simplify = FALSE)
outP <- do.call(rbind.data.frame, outP)
outP <- outP[order(outP$sdENERGY/outP$ENERGY, decreasing = FALSE),]
write.csv(outP, quote = FALSE, row.names = FALSE, file = "results.csv")
