################
##
## @description Generate stochastic data under the True Model
##
## @param deterministic Deterministic value
## @param sigma         Data error variance
##
## @return Stochastic data
##
## @lastChange 2016-12-28
##
## @changes
##
################
generateY <- function(deterministic, sigma){
  
  EY <- deterministic + rnorm(length(deterministic), 0, sigma)
  
  normalizedEY <- (EY - mean(EY) ) / sd(EY)
  
  return(normalizedEY)
}
