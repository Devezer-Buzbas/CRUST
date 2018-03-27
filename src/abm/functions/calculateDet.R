################
##
## @description Calculate the deterministic part of the model
##
## @param model   Model in matrix format
## @param xset    X values randomly generated
## @param weights Beta weights
## @param betas   Random betas
##
## @return Deterministic value of the model
##
## @lastChange 2018-03-01
##
## @changes
##   Fixed the use of wrong beta values [2018-03-01]
##   Included parameter weights [2018-03-01]
##   Adjust model [2017-02-13]
##
################
calculateDet <- function(model, xset, weights, betas){
  
  if(!is.matrix(model)){
    model <- t(as.matrix(model))
  }
  
  k <- length(model[1,])
  
  deterministic <- 0
  f <- 10^((k - 1):0)
  for(r in 1:nrow(model)){
    x <- rowProds(as.matrix(xset[, model[r,] == 1]))
    
    index <- weights[, 1] == sum(as.numeric(model[r,] == 1) * f)
    
    deterministic <- deterministic + (betas[which(index == TRUE)] * x)
  }
  
  return(deterministic)
}
