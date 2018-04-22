################
##
## @description Calculate the distance between the betas of two models
##
## @param betas1 Betas of model 1
## @param betas2 Betas of model 2
##
## @return Distance between models
##
## @lastChange 2016-12-28
##
## @changes
##
################
calculateDistance <- function(betas1, betas2){
  
  dist <- sum((betas1 - betas2)^2)
  
  return(dist)
}
