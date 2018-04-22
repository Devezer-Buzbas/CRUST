################
##
## @description Search index of the model
##
## @param model  Model in matrix format
## @param models List of all possible models
##
## @return Index of the model
##
## @lastChange 2016-12-28
##
## @changes
##
################
searchModel <- function(model, models){
  
  if(!is.matrix(model)){
    model <- t(as.matrix(model))
  }
  
  for(i in 1:length(models)){
    
    if(is.matrix(models[[i]])){
      ms <- models[[i]]
    } else {
      ms <- t(as.matrix(models[[i]]))
    }
    
    if(compareModels(model, ms)){
      return(i)
    }
  }
  
  return(0)
}
