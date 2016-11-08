library(sybilSBML)
library(foreach)
library(doMC)
registerDoMC(4)
model <- readSBMLmod("Documents/masterThesis/Results/matureAstrocyte.xml")

pOF <- function(model,R){
  t0 <- optimizeProb(model)@lp_obj
  set.seed(1234)
  V <- sapply(seq_len(model@react_num), function(x){runif(model@react_num,model@lowbnd[x],model@uppbnd[x])})
  diag(V) <- getFluxDist(optimizeProb(model))
  OF <- function(w){
    f <- mapply(function(x,y){V[x,y]},x=w,y=seq_len(model@react_num))
    model@lowbnd <- f
    model@uppbnd <- f
    optimizeProb(model)@lp_obj
  }
  b <- unlist(foreach(i=1:R) %dopar% {
    OF(sample(seq_len(model@react_num)))
  })
  return(mean(b<=t0))
}

pOF(model,1000)

sdM <- function(model,R){
  t0 <- optimizeProb(model)@lp_obj
  set.seed(1234)
  FVA <- fluxVar(model)
  V <- sapply(seq_len(model@react_num), function(x){runif(n = FVA@react@react_num,min = min(c(FVA@lp_obj[x],FVA@lp_obj[(model@react_num+x)])),max = max(c(FVA@lp_obj[x],FVA@lp_obj[(model@react_num+x)])) )})
  diag(V) <- getFluxDist(optimizeProb(model))
  OF <- function(w){
    f <- mapply(function(x,y){V[x,y]},x=w,y=seq_len(model@react_num))
    model@lowbnd <- f
    model@uppbnd <- f
    optimizeProb(model)@lp_obj
  }
  b <- unlist(foreach(i=1:R) %dopar% {
    OF(sample(seq_len(model@react_num)))
  })
  return(c(BIAS=mean(mean(b)-b)^2,SD=sd(b)))
}


tM <- function(model1,model2,R){
  V <- function(model,R){
    t0 <- optimizeProb(model)@lp_obj
    set.seed(1234)
    FVA <- fluxVar(model)
    V <- sapply(seq_len(model@react_num), function(x){runif(n = FVA@react@react_num,min = min(c(FVA@lp_obj[x],FVA@lp_obj[(model@react_num+x)])),max = max(c(FVA@lp_obj[x],FVA@lp_obj[(model@react_num+x)])) )})
    diag(V) <- getFluxDist(optimizeProb(model))
    OF <- function(w){
      f <- mapply(function(x,y){V[x,y]},x=w,y=seq_len(model@react_num))
      model@lowbnd <- f
      model@uppbnd <- f
      optimizeProb(model)@lp_obj
    }
    b <- unlist(foreach(i=1:R) %dopar% {
      OF(sample(seq_len(model@react_num)))
    })
  }
  M1 <- V(model1,R)
  M2 <- V(model2,R)
  t.test(M1,M2)
}

healthy <- readSBMLmod("Documents/masterThesis/Results/matureAstrocyte.xml")
inflammated <- healthy
lowbnd(inflammated)[inflammated@react_id=="EX_hdca(e)"] <- -0.208
uppbnd(inflammated)[inflammated@react_id=="EX_hdca(e)"] <- -0.208
tM(healthy,inflammated,1000)
