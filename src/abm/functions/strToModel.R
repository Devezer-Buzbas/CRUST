################
##
## @description Convert a model from string to matrix format
##
## @param modelStr  Model in string format
## @param k         Number of factors
##
## @return Model in matrix format
##
## @lastChange 2016-12-28
##
## @changes
##
################
strToModel <- function(modelStr, k){
  tokens <- unlist(strsplit(tolower(modelStr), "[+]"))
  model <- matrix(0, nrow=length(tokens), ncol=k)
  i <- 1
  for(terms in tokens){
    for(factor in unlist(strsplit(terms, "x"))){
      if(factor != ""){
        model[i, as.numeric(as.character(factor))] <- 1
      }
    }
    i <- i + 1
  }
  
  return(model)
}
