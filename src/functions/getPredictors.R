################
##
## @description Get predictors of all models
##
## @param models List of models in matrix format
## @param k      Number of factors
##
## @return List of predictors
##
## @lastChange 2016-12-28
##
## @changes
##
################
getPredictors <- function(models){
  predictors <- array(0, dim=length(models))
  fPredictors <- list()
  for(i in 1:length(models)){
    model <- models[[i]]
    if(!is.matrix(model)){
      model <- t(as.matrix(models[[i]]))
    }
    
    k <- length(model[1,])
    
    nTerms <- 0
    factors <- c()
    for(r in 1:nrow(model)){
      if(sum(model[r,]) == 1){
        nTerms <- nTerms + 1
        factors <- c(factors, sum(model[r,] * 1:k))
      }
    }
    predictors[i] <- nTerms
    factors <- paste0(list(factors))
    if(nchar(factors) == 1){
      fPredictors[[i]] <- factors
    } else {
      fPredictors[[i]] <- substr(factors, start=3, stop=(nchar(factors) - 1))
    }
  }
  
  return(list(predictors, fPredictors))
}
