#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @export plotFluxes
#' 
plotFluxes <- function(model,maxiter=100,...){
  S <- as.matrix(t(model@S))
  colnames(S) <- model@met_id
  rownames(S) <- model@react_id
  D <- rep(0,model@react_num)
  fDistribution <- getFluxDist(optimizeProb(model))
  fluxes <- matrix(fDistribution,ncol = 1,dimnames = list(model@react_id,c()))
  D[fluxes<1] <- -1
  D[fluxes>1] <- 1
  S <- S*D
  S[S < 0] <- -1
  S[S > 0] <- 1
  S <- S[fDistribution!=0,]
  S <- S[,colSums(S!=0)!=0]
  types <- structure(c(rep(1,length(rownames(S))),rep(0,length(colnames(S)))),names=c(rownames(S),colnames(S)))
  edges <- NULL
  for(i in names(types)[types==1]){
    n_mets <- names(S[i,])
    if (length(grep(2,S[i,]))>0){
      edges <- c(edges,unlist(sapply(n_mets[S[i,]!=0],function(met){c(met,i)})))
      edges <- c(edges,unlist(sapply(n_mets[S[i,]!=0],function(met){c(i,met)})))
    } else{
      edges <- c(edges,unlist(sapply(n_mets[S[i,]<0],function(met){c(met,i)})))
      edges <- c(edges,unlist(sapply(n_mets[S[i,]>0],function(met){c(i,met)})))
    }
  }
  edges <- as.vector(unlist(as.vector(sapply(edges, function(edge){which(names(types)==edge)}))))
  g <- make_bipartite_graph(types = types,edges = edges,directed = TRUE)
  V(g)$name <- names(types)
  reactions <- types==1
  edges <- get.edgelist(g)
  input <- edges[edges[,2]%in% V(g)[reactions]$name,]
  output <- edges[edges[,1]%in% V(g)[reactions]$name,]
  match1 <- match(output[,2],input[,1])
  match2 <- match(input[,1],output[,2])
  networkConnections <- rbind(
    cbind(output[,],input[match1,2]),
    cbind(output[match2,], input[,2])
  )
  networkConnections <- networkConnections[complete.cases(networkConnections),]
  networkConnections <- unique(networkConnections)
  networkConnections <- data.frame(networkConnections[,c(1,3,2)])
  names(networkConnections) <- c("from", "to", "metabolite")
  g <-  simplify(graph.data.frame(networkConnections),edge.attr.comb="c")
  flux <- abs(fluxes[V(g)$name,])
  V(g)$color <- ifelse(V(g)$name%in%findExchReact(model)@react_id,"green","gray90")
  topology <- layout_with_dh(g,weight.node.edge.dist=(flux/max(flux))*10)
  plot(g, vertex.size=(flux/max(flux))*10, edge.arrow.size=0.4, layout=topology)
}
