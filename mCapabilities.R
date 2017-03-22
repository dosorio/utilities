library(foreach)
library(doMC)
library(sybilSBML)
registerDoMC(4)

mCapabilities <- function(model){
  mc <- data.frame(foreach(i=model@met_id,.combine = rbind) %:% foreach(j=model@met_id) %dopar% {
    model@obj_coef <- rep(0,model@react_num)
    model <- addReact(model = model,
                      id = "MC",
                      Scoef = c(-1,1),
                      met = c(i,j),
                      reversible = FALSE,
                      ub = 1000,
                      obj = 1)
    suppressMessages(optimizeProb(model)@lp_obj)
  })
  colnames(mc) <- model@met_id
  rownames(mc) <- model@met_id
  return(mc)
}