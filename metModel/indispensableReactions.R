indispensableReactions <- function(model){
  indispensable <- NULL
  pb <- txtProgressBar(min = 1,max = model@react_num,style=3)
  for (ID in model@react_id){
    reactionNumber <- (react_id(model) == ID)
    setTxtProgressBar(pb, grep("TRUE",reactionNumber))
    lb <- lowbnd(model)[reactionNumber]
    ub <- uppbnd(model)[reactionNumber]
    lowbnd(model)[reactionNumber] <- 0
    uppbnd(model)[reactionNumber] <- 0
    if (optimizeProb(model)@lp_obj==0){
      indispensable <- c(indispensable,ID)
    }
    lowbnd(model)[reactionNumber] <- lb
    uppbnd(model)[reactionNumber] <- ub
  }
  close(pb)
  return(indispensable)
}
