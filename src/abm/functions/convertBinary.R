################
##
## @description Convert number into binary format
##
## @param v Vector of positions to put value 1
## @param k Number of factors
##
## @return Binary vector
##
## @lastChange 2016-12-28
##
## @changes
##
################
convertBinary <- function(v, k){
  r <- rep(0, k)
  r[v] <- 1
  
  return(r)
}
