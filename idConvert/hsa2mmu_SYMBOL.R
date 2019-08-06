hsa2mmu_SYMBOL <- function(Z){
  Z <- sapply(Z, function(I){
    I <- unlist(strsplit(I,''))
    I <- c(I[1],tolower(I[-1]))
    I <- paste0(I, collapse = "")
    return(I)
  })
  Z <- as.character(Z)
  return(Z)
}
