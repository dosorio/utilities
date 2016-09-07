#' @export extractGenes
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Return list of unique genes from a set of GPR
extractGenes <- function(GPR){
  GPR <- as.vector(GPR)
  GPR <- unique(unlist(strsplit(gsub("\\(|and|or|\\)","",GPR),"[[:blank:]]+")))
  return(GPR)
}