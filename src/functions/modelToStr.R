################
##
## @description Convert a model from matrix to string format
##
## @param model  Model in matrix format
##
## @return Model in string format
##
## @lastChange 2016-12-28
##
## @changes
##
################
modelToStr <- function(model){
  f <- NULL
  
  if(!is.matrix(model)){
    model <- t(matrix(model))
  }
  
  for(term in 1:nrow(model)){
    s <- NULL
    index <- 1
    for(item in model[term,]){
      if(item == 1){
        if(is.null(s)){
          s <- paste0("X", index)
        } else {
          s <- paste0(s, ":X", index)
        }
      }
      index <- index + 1
    }
    if(is.null(f)){
      f <- c("Y ~", s)
    } else {
      f <- c(f, c("+", s))
    }
  }
  
  return(paste0(f, collapse=" "))
}
