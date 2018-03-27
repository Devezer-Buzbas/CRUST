################
##
## @description Generate the predictors values
##
## @param n           Sample size
## @param k           Number of factors
## @param correlation Sample correlation
##
## @return Predictors value
##
## @lastChange 2016-12-28
##
## @changes
##
################
generateXSet <- function(n, k, correlation){
  
  xperm <- matrix(data=0, nrow=n, ncol=k)
  x <- seq(1, n, 1)
  xperm[, 1] <- x
  xperm[, 1] <- (xperm[, 1] - mean(xperm[, 1]) ) / sd(xperm[, 1])
  
  for(i in 2:k){
    newx <- cbind(x, runif(n)) # fixed and some random data values in a matrix
    xctr <- scale(newx, center=TRUE, scale=FALSE) # center both x and newdata
    Q <- qr.Q(qr(xctr[,1,drop=FALSE])) # QR decomposition, take just Q
    P <- tcrossprod(Q) # QQ', projection onto space defined by x
    xort <- (diag(n) - P) %*% xctr[,2] # make x centered orthogonal to new data centered
    xmat <- cbind(xctr[,1], xort)
    xscal <- xmat %*% diag(1 / sqrt(colSums(xmat^2))) # scale vectors to length 1
    xperm[, i] <- xscal[,2] + (1 / tan(acos(correlation))) * xscal[,1] # correlated vector to x
    
    xperm[, i] <- (xperm[, i] - mean(xperm[, i]) ) / sd(xperm[, i])
  }
  
  xpermsorted <- data.frame(xperm)
  
  return(xpermsorted)
}
