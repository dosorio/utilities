plotDR <- function(KO, labelGenes = 'FDR', boldGenes = NULL, title = NULL, subtitle = NULL){
  library(ggplot2)
  library(ggrepel)
  o <- -log10(KO$diffRegulation$p.value)
  e <- -log10(pchisq(sort(rchisq(nrow(KO$diffRegulation), df = 1), decreasing = TRUE), df = 1, lower.tail = FALSE))
  dF <- data.frame(X = e, Y = o, geneID = KO$diffRegulation$gene)
  dF <- as.data.frame.array(dF)
  if(labelGenes == 'FDR'){
    dF$geneID[KO$diffRegulation$p.adj > 0.05] = ''
  }
  if(labelGenes == 'P'){
    dF$geneID[KO$diffRegulation$p.value > 0.05] = ''
  }
  geneColor <- ifelse(KO$diffRegulation$p.value < 0.05, 'red', 'black')
  genePoint <- 16
  plotQQ <- ggplot(dF, aes(X,Y, label = geneID)) + 
    geom_point(color = geneColor, pch = genePoint) + 
    theme_bw() + 
    geom_text_repel(segment.color = 'gray60', segment.alpha = 0.5, max.iter = 1e4, aes(fontface = ifelse(dF$geneID %in% boldGenes,2,1)), box.padding = .08) + 
    labs(y=expression(-log[1*0]*" (Observed P-values)"),x=expression(-log[1*0]*" (Expected P-values)"), title = title, subtitle = subtitle)
  print(plotQQ)
}
