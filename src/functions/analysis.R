################
##
## @description Calculate the statistics for the Selected and Global models
##              assuming data generated under the True Model (Yset),
##              X values randomly generated (Xset) and betas weight
##              (betas)
##
## @param sModel  Selected Model in matrix format
## @param gModel  Global Model in matrix format
## @param yset    Data generated under the True Model
## @param xset    X values randomly generated
## @param weights Betas weight
##
## @return Selected model (model) and Global model (gModel) information:
##           betasEst    Betas estimate
##           betasErr    Betas error
##           tstatistics T Statistic
##           rsq         R-Squared
##           arsq        Adjusted R-Squared
##           aic         Aikeka Information Criterion
##           bic         Bayesian Information Criterion
##
## @lastChange 2018-02-28
##
## @changes
##    Added BIC analysis [2018-02-28]
##    Fixed betas estimate and error [2017-01-18]
##    Added comments [2016-12-28]
##
################
analysis <- function(sModel, gModel, yset, xset, weights){
  
  data <- cbind(Y=yset, xset)
  
  ## Linear regression
  slmy <- lm(modelToStr(sModel), data=data)
  glmy <- lm(modelToStr(gModel), data=data)
  
  ## Summary of model fitting
  slm.res <- summary(slmy)
  glm.res <- summary(glmy)
  
  ## Selected Model Betas
  k <- ncol(xset)
  f <- 10^((k - 1):0)
  
  if(!is.matrix(sModel)){
    sModel <- t(as.matrix(sModel))
  }
  
  sCoef <- coef(slm.res)
  sName <- rownames(sCoef)
  
  sBetasEst <- rep(0, length(weights[, 1]))
  sBetasErr <- rep(0, length(weights[, 1]))
  for(n in sName){
    if(n == "(Intercept)"){
      sBetasEst[1] <- sCoef[1, 1]
      sBetasErr[1] <- sCoef[1, 2]
    } else {
      index <- weights[, 1] == sum(strToModel(gsub(":", "", n), k)[1,] * f)
      sBetasEst[index] <- sCoef[which(sName == n), 1]
      sBetasErr[index] <- sCoef[which(sName == n), 2]
    }
  }
  
  ## Global Model Betas
  if(!is.matrix(gModel)){
    gModel <- t(as.matrix(gModel))
  }
  
  gCoef <- coef(glm.res)
  gName <- rownames(gCoef)
  
  gBetasEst <- rep(0, length(weights[, 1]))
  gBetasErr <- rep(0, length(weights[, 1]))
  for(n in gName){
    if(n == "(Intercept)"){
      gBetasEst[1] <- gCoef[1, 1]
      gBetasErr[1] <- gCoef[1, 2]
    } else {
      index <- weights[,1] == sum(strToModel(gsub(":", "", n), k)[1,] * f)
      gBetasEst[index] <- gCoef[which(gName == n), 1]
      gBetasErr[index] <- gCoef[which(gName == n), 2]
    }
  }
  
  ## Selected Model statistics
  selectedModel <- list(
      betasEst = sBetasEst,
      betasErr = sBetasErr,
      tstatistics = coef(slm.res)[2, 3],
      rsq = slm.res$r.squared,
      arsq = slm.res$adj.r.squared,
      aic = AIC(slmy),
      bic = BIC(slmy))
  
  ## Global Model statistics
  globalModel <- list(
      betasEst = gBetasEst,
      betasErr = gBetasErr,
      tstatistics = coef(glm.res)[2, 3],
      rsq = glm.res$r.squared,
      arsq = glm.res$adj.r.squared,
      aic = AIC(glmy),
      bic = BIC(glmy))
  
  return(list("model"=selectedModel, "gModel"=globalModel))
}
