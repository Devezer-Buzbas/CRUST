################
##
## @description Generate random betas
##
## @param model   Model in matrix format
## @param weights Beta weights
## @param sigma   Data error variance
##
## @return Random betas
##
## @lastChange 2016-12-28
##
## @changes
##
################
getBetas <- function(model, weights, sigma){
  
  if(!is.matrix(model)){
    model <- t(matrix(model))
  }
  
  k <- length(model[1,])
  
  nBetas <- length(weights[, 1])
  
  ## Get the weight of Betas
  f <- 10^((k - 1):0)
  betasWeight <- rep(0, length(weights[, 1]))
  for(r in 1:nrow(model)){
    index <- weights[, 1] == sum(as.numeric(model[r,] == 1) * f)
    betasWeight[index] <- weights[index, 2]
  }
  
  ## Calculate random betas
  calcBetas <- sign(runif(nBetas, -1, 1)) * sqrt(rdirichlet(1, betasWeight)) * (1 - sigma)
  
  return(calcBetas)
}
