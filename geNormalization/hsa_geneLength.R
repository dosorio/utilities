hsa_geneLenght <- function(){
  require(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <-  useDataset("hsapiens_gene_ensembl",mart=ensembl)
  gLength <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "gene_biotype"), mart = ensembl)
  gL <- abs(gLength[,3] - gLength[,4])
  gLength$gL <- gL
  gLength[,c(1,2,3,4,6,5)]
}
