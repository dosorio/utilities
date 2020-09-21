downloadHarmonizeGeneSet <- function(url, tableList = c(1)){
  require(xml2)
  require(XML)
  tFile <- tempfile()
  download_html(url, tFile)
  D <- readHTMLTable(tFile)
  D <- sapply(seq_along(D)[(tableList+1)], function(tNumber){
    D[[tNumber]]$Symbol
  })
  D <- unique(unlist(D))[-1]
  return(D)
}
