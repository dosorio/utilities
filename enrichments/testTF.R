# testTF <- function(X){
#   
# }

load('~/../Desktop/GC1.RData')
X <- O$diffRegulation$Z
names(X) <- toupper(O$diffRegulation$gene)
X <- X[X > 0]

CHEA2016 <- readLines('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2016')
CHEA2016 <- strsplit(CHEA2016, '\t')
namesCHEA2016 <- lapply(CHEA2016, function(X){X[1]})
CHEA2016 <- lapply(CHEA2016, function(X){X[-c(1,2)]})
names(CHEA2016) <- namesCHEA2016
CHEA2016[2]

JASPAR2016 <- readLines('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=TRANSFAC_and_JASPAR_PWMs')
JASPAR2016 <- strsplit(JASPAR2016, '\t')
namesJASPAR2016 <- lapply(JASPAR2016, function(X){X[1]})
JASPAR2016 <- lapply(JASPAR2016, function(X){X[-c(1,2)]})
names(JASPAR2016) <- namesJASPAR2016

ALL <- c(CHEA2016, JASPAR2016)
library(fgsea)

E <- fgsea(pathways = ALL, stats = X, nperm = 1e4)
E <- fgseaMultilevel(pathways = ALL, stats = X)
fgsea::plotEnrichment(CHEA2016[['FOXM1 23109430 ChIP-Seq U2OS Human']], stats = X)
