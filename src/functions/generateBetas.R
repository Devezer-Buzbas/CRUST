################
##
## @description Set the weights of all betas
##
## @param models List of models in matrix format
##
## @return Beta weights
##
## @lastChange 2017-02-14
##
## @changes
##
################
generateBetas <- function(models){
  ## Identify the number of terms the full model has
  fModel <- t(as.matrix(models[[1]]))
  for(model in models){
    
    if(!is.matrix(model)){
      model <- t(as.matrix(model))
    }
    
    if(nrow(model) > nrow(fModel)){
      fModel <- model
    }
  }
  beta <- c(1, rep(1, nrow(fModel)))
  betas <- cbind(0, 1)
  
  k <- length(fModel[1,])
  
  f <- 10^((k-1):0)
  for(r in 1:nrow(fModel)){
    index <- sum(fModel[r,] * f)
    betas <- rbind(betas, cbind(as.numeric(index), as.numeric(beta[r + 1])))
  }
  
  return(betas)
}
