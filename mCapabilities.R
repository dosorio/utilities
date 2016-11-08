library(foreach)
library(doMC)
library(sybilSBML)
registerDoMC(32)

mCapabilities <- function(model){
  foreach(i=model@met_id,.combine = rbind) %:% foreach(j=model@met_id) %dopar% {
    model@obj_coef <- rep(0,model@react_num)
    model <- addReact(model = model,
                      id = "MC",
                      Scoef = c(-1,1),
                      met = c(i,j),
                      reversible = FALSE,
                      ub = 1000,
                      obj = 1)
    suppressMessages(optimizeProb(model)@lp_obj)
  }
}

healthy <- readSBMLmod("matureAstrocyte.xml")
healthyR <- mCapabilities(healthy)
write.table(healthyR,file = "healthy.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

inflammated <- readSBMLmod("matureAstrocyte.xml")
lowbnd(inflammated)[react_id(inflammated) == 'EX_hdca(e)'] <- -0.208
uppbnd(inflammated)[react_id(inflammated) == 'EX_hdca(e)'] <- -0.208
inflammatedR <- mCapabilities(inflammated)
write.table(inflammatedR,file = "inflammated.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

tibolone <- readSBMLmod("matureAstrocyte_Tibolone.xml")
lowbnd(tibolone)[react_id(tibolone) == 'EX_hdca(e)'] <- -0.208
uppbnd(tibolone)[react_id(tibolone) == 'EX_hdca(e)'] <- -0.208
mCapabilities(tibolone,"tibolone.txt")