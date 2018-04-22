################
##
## @description Compare two models
##
## @param model1 Model 1 in matrix format
## @param model2 Model 2 in matrix format
##
## @return True if models are equal, False otherwise
##
## @lastChange 2016-12-28
##
## @changes
##
################
compareModels <- function(model1, model2){
  
  if(!is.matrix(model1)){
    model1 <- t(as.matrix(model1))
  }
  
  if(!is.matrix(model2)){
    model2 <- t(as.matrix(model2))
  }
  
  k <- min(length(model1[1,]), length(model2[1,]))
  
  f <- 10^((k-1):0)
  
  match <- FALSE
  if(nrow(model1) == nrow(model2)){
    match <- TRUE
    for(r in 1:nrow(model1)){
      l <- sum(model1[r,] * f)
      found <- FALSE
      for(j in 1:nrow(model2)){
        if(l == sum(model2[j,] * f)){
          found <- TRUE
        }
      }
      
      if(!found){
        match <- FALSE
      }
    }
  }
  
  return(match)
}
