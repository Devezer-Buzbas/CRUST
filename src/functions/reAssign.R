reAssign <- function(x){
  # Description:
  # ReAssign models for k=3 according to complexity defined by 
  # by geometric shapes, from Gustavo's way of creating models
  # for k=3.
  # x: L x L matrix of statistics, where
  # Model 11 should replace model 10
  # Model 12 should replace model 11
  # Model 13 should replace model 12
  # Model 14 should replace model 13
  # Model 10 should replace model 14
  
  xx = x # keep the input matrix intact
  xx[10:13,] = x[11:14,]
  xx[14,] = x[10,]
  reassignmatrix = xx
  reassignmatrix[, 10:13] = xx[, 11:14]
  reassignmatrix[, 14] = xx[, 10]
  
  return(reassignmatrix)
}