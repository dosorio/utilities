library(clustermole)
library(fgsea)

mmPanglao <- clustermole_markers(species = 'mm')
mmPanglao <- mmPanglao[mmPanglao$species %in% 'Mouse',]
mmPanglao <- mmPanglao[mmPanglao$db %in% 'PanglaoDB',]
ctList <- unique(mmPanglao$celltype)
ctList <- sapply(ctList, function(X){
  paste0(c(X, unique(mmPanglao$gene[mmPanglao$celltype %in% X])), collapse = '\t')
})
writeLines(ctList, 'mmPanglaoDB.gmt')


hsaPanglao <- clustermole_markers(species = 'hs')
hsaPanglao <- hsaPanglao[hsaPanglao$species %in% 'Human',]
hsaPanglao <- hsaPanglao[hsaPanglao$db %in% 'PanglaoDB',]
ctList <- unique(hsaPanglao$celltype)
ctList <- sapply(ctList, function(X){
  paste0(c(X, unique(hsaPanglao$gene[hsaPanglao$celltype %in% X])), collapse = '\t')
})
writeLines(ctList, 'hsaPanglaoDB.gmt')
