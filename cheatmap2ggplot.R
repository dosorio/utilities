cheatmap2ggplot <- function(X){
  ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(X)))
}
