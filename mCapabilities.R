mCapabilities <- function(model){
  capabilities <- matrix(data = 0,
                         nrow = model@met_num,
                         ncol = model@met_num,
                         dimnames = list(sort(model@met_id),sort(model@met_id))
                         )
  for(i in model@met_id){
    for(j in model@met_id){
      model <- addReact(model = model,
                        id = "MC",
                        Scoef = c(-1,1),
                        met = c(i,j),
                        reversible = FALSE,
                        ub = 1000,
                        obj = 1)
      capabilities[i,j] <- optimizeProb(model)@lp_obj
    }
  }
  return (capabilities)
}