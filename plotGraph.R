plotNetwork <- function(matrix, symmetric, outFile, ...) {
  require(igraph)
  require(reshape2)
  diag(matrix) <- NA
  threshold <- quantile(abs(matrix), 0.9, na.rm = TRUE)
  matrix[abs(matrix) < threshold] <- NA
  if (isTRUE(symmetric)) {
    matrix[upper.tri(matrix, diag = TRUE)] <- NA
  }
  matrix <- melt(as.matrix(matrix))
  matrix <- matrix[complete.cases(matrix), ]
  colnames(matrix) <- c("from", "to", "weight")
  if (isTRUE(symmetric)) {
    graph <- graph_from_data_frame(matrix, directed = FALSE)
  } else {
    graph <- graph_from_data_frame(matrix, directed = TRUE)
  }
  oLayout <- layout.fruchterman.reingold(graph, niter = 1000)
  vColor <- rgb(0, 0.5, 1, 0.3)
  eColor <- ifelse(E(graph)$weight > 0, yes = "red", no = "blue")
  png(
    outFile,
    width = 3000,
    height = 3000,
    res = 300,
    pointsize = 20
  )
  par(mar = c(0, 0, 0, 0))
  plot(
    graph,
    layout = oLayout,
    vertex.size = 10,
    edge.arrow.size = 0.1,
    edge.color = eColor,
    vertex.label.cex = 0.5,
    vertex.color = vColor,
    vertex.frame.color = vColor,
    vertex.label.family = "Arial",
    vertex.label.color = "black",
    ...
  )
  dev.off()
}
