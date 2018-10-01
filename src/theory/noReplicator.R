################
##
## @description Population with No Replicator
##              Model: Markov chain (first order, satisfying Markov property)
##              Analysis: using probability calculus and
##              Monte Carlo estimates for model comparison
##
## @param None
##
## @return None
##
## @lastChange 2018-03-26
##
## @changes
##
################
# begin all
#-------------------------------------------------------------
library(caTools)
library(data.table)
library(permute)
library(matrixStats)
library(MCMCpack)
library(lattice)
library(ggplot2)
library(directlabels)
library(tibble)
library(reshape2)
library(gridExtra)


#############
## PATHS
#############
baseDir <- "."
scriptDir <- paste0(baseDir, "/src/functions")
inputDir <- paste0(baseDir, "/data/modelComparison")
outputDir <- paste0(baseDir, "/data/plot")


#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/calculateDet.R"))
source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/strToModel.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/generateY.R"))
source(paste0(scriptDir, "/generateXSet.R"))
source(paste0(scriptDir, "/getBetas.R"))
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/reAssign.R"))


###################
## INPUT PARAMETERS
###################
## Number of factors
k <- 3

## Sigma (Error variance)
sigma <- 0.2 #(1:4) error variance to model expectation
sigmaTxt <- gsub("\\.", "", toString(sigma))

## Generate all models
models = generateModels(k)
models[[1]] = t(as.matrix(models[[1]]))

## Number of models
L = length(models)

#-----------------------------------
# Scientist model proposal matrices: element ij is the probability of
# proposing model j if the current global is model i

## Tess matrix
tessMatrix <- matrix(data=0.0, L, L)
for(index in 1:L){
  tessMatrix[index,] <- modelSimilarByTerm(models[[index]], models,
                        mode="all", modelSelection="soft")
  nSimModels <- sum(tessMatrix[index,])
  tessMatrix[index,] <- tessMatrix[index,] * (1 / (nSimModels + 1))
  tessMatrix[index, which(tessMatrix[index,] == 0)] <- (1 -
        sum(tessMatrix[index,])) / (length(tessMatrix[index,]) - nSimModels)
}

## Bo matrix
boMatrix <- matrix(data=0.0, L, L)
for(index in 1:L){
  boMatrix[index,] <- modelSimilarByInteraction(models[[index]], models,
                      mode="all", modelSelection="soft")
  nSimModels <- sum(boMatrix[index,])
  boMatrix[index,] <- boMatrix[index,] * (1 / (nSimModels + 1))
  boMatrix[index, which(boMatrix[index,] == 0)] <- (1 -
        sum(boMatrix[index,])) / (length(boMatrix[index,]) - nSimModels)
}

## Mave Matrix
maveMatrix <- matrix(data=0.0, L, L)
for(index in 1:L){
  maveMatrix[index,] <- 1 / length(maveMatrix[index,])
}

# Scientist populations

# Allscientists at equal proportion
all <- (tessMatrix + boMatrix + maveMatrix) * (1 / 3)

# Tess dominant population
tess <- (tessMatrix * 0.99) + (boMatrix + maveMatrix) * 0.005

# Bo dominant population
bo <- (boMatrix * 0.99) + (tessMatrix + maveMatrix) * 0.005

# Mave dominant population
mave <- (maveMatrix * 0.99) + (boMatrix + tessMatrix) * 0.005

#------------------------------------------------------------
# Load and modify the probability that a model's score compares with respect to
# another model's score, given a true model, when model comparison statistic
# is AIC or SC
#
# Paic: List whose element l is for true model l
# for model selection statistic AIC.
#
# Psc: List whose element l is for true model l
# for model selection statistic AIC.
#
# Each element of the list is a matrix of probabilities,
# element ij is the model j having lower model selection score than model i:
# that is P(S(M_i)>S(M_j)|true model l).
# Diagonal, probability from a model to self is 1 by convention.
#
Paic <- list()
Psc <- list()
for (l in 1:L) {
  # for AIC unzip the ModelComparisonProbabilitiesSigma02AIC
  # look up tables (generates with N=1e5 samples) in the appropriate directory.
  # Alternatively use
  # ModelComparisonProbabilitiesByMonteCarlo.R coupled with
  # the function GetAICConstant.R function to regenerate these files for a desired
  # precision.
  name <- paste0(inputDir,
      '/ModelComparisonProbabilitiesSigma', sigmaTxt, 'AIC', l, '.RData',
      sep="")
  load(name)

  # SC
  Paic[[l]] <- ProbMC
  name <- paste0(inputDir,
      '/ModelComparisonProbabilitiesSigma', sigmaTxt, 'Schwarz', l, '.RData',
      sep="")
  load(name)
  Psc[[l]] <- ProbMC
}
#------------------------------------------------------------
# Transition probability matrices for
# a: 4 scientist populations: All, Tess, Bo, Mave
# m: 14 true models
# P: (14,14,a,m) array. P(,,a,m) is the transition
# probability matrix for scientist population a, true model m,
# model selection statistic is AIC
#
# Q: (14,14,a,m) array. Q(,,a,m) is the transition
# probability matrix for scientist population a, true model m,
# model selection statistic is SC
#
P <- array(0, dim = c(L, L, 4, L)) # AIC
for (i in 1: L) {
  p <- all * Paic[[i]]; P[,,1,i] <- p / rowSums(p); # normalize so that probabilities sum to 1

  p <- tess*Paic[[i]]; P[,,2,i] <- p / rowSums(p);

  p <- bo*Paic[[i]]; P[,,3,i] <- p / rowSums(p);

  p <- mave*Paic[[i]]; P[,,4,i] <- p / rowSums(p);
}
#
Q <- array(0, dim = c(L, L, 4, L)) # SC
for (i in 1:L) {
  q <- all * Psc[[i]]; Q[,,1,i] = q / rowSums(q);

  q <- tess*Psc[[i]]; Q[,,2,i] = q / rowSums(q);

  q <- bo*Psc[[i]]; Q[,,3,i] = q / rowSums(q);

  q <- mave*Psc[[i]]; Q[,,4,i] = q/rowSums(q);
}
#------------------------------------------------------------
#------------------------------------------------------------
# A. Calculate the Transient Properties of the System
#------------------------------------------------------------
# begin Stickiness
#
# 1. Stickiness.
# Mean (over proposed models) probability of staying in a global model j
# conditional on a true model i. The diagonal is for the true model
# One matrix of probabilities for each scientist population.

MeanProbStayInModelAllAIC <- rep(0, L) # All scientists in equal proportions, using AIC
MeanProbStayInModelTessAIC <- rep(0, L)
MeanProbStayInModelBoAIC <- rep(0, L)
MeanProbStayInModelMaveAIC <- rep(0, L)
#
MeanProbStayInModelAllSC <- rep(0, L) # All scientists in equal proportions, using SC
MeanProbStayInModelTessSC <- rep(0, L)
MeanProbStayInModelBoSC <- rep(0, L)
MeanProbStayInModelMaveSC <- rep(0, L)
#
for (i in 1:L) {
# for each true model transpose of Paic (and PSc) gives
# the prob of model staying as global, not being beaten
# by the proposed model: Paic[[l]][i,j]
# is equal to P(S(M_j)<S(M_i)|l is true model)
#
MeanProbStayInModelAllAIC[i] <- rowSums(all * t(Paic[[i]]))[i]
MeanProbStayInModelTessAIC[i] <- rowSums(tess * t(Paic[[i]]))[i]
MeanProbStayInModelBoAIC[i] <- rowSums(bo * t(Paic[[i]]))[i]
MeanProbStayInModelMaveAIC[i] <- rowSums(mave * t(Paic[[i]]))[i]
#
MeanProbStayInModelAllSC[i] <- rowSums(all * t(Psc[[i]]))[i]
MeanProbStayInModelTessSC[i] <- rowSums(tess * t(Psc[[i]]))[i]
MeanProbStayInModelBoSC[i] <- rowSums(bo * t(Psc[[i]]))[i]
MeanProbStayInModelMaveSC[i] <- rowSums(mave * t(Psc[[i]]))[i]
}

#-----------------------------------------
# Produce Plots
modelnum <- c(1:L)
#
MeanProbStayInModelmatrixAIC <- matrix(c(MeanProbStayInModelAllAIC,
                                         MeanProbStayInModelBoAIC,
                                         MeanProbStayInModelMaveAIC,
                                         MeanProbStayInModelTessAIC),
                                         nrow = 4, ncol = L, byrow = TRUE)
#
MeanProbStayInModelmatrixSC <- matrix(c(MeanProbStayInModelAllSC,
                                        MeanProbStayInModelBoSC,
                                        MeanProbStayInModelMaveSC,
                                        MeanProbStayInModelTessSC),
                                        nrow = 4, ncol = L, byrow = TRUE)
#------------------------------------------------------------------
# Next bit of code rearranges models for the heatmap
# so that the model complexity is in increasing order
#
zAIC <- MeanProbStayInModelmatrixAIC
zAIC[, 14] <- MeanProbStayInModelmatrixAIC[, 10]
zAIC[, 10:13] <- MeanProbStayInModelmatrixAIC[, 11:14]
#
zSC <- MeanProbStayInModelmatrixSC
zSC[, 14] <- MeanProbStayInModelmatrixSC[, 10]
zSC[, 10:13] <- MeanProbStayInModelmatrixSC[, 11:14]
#------------------------------------------------------------------
# Heatmap
x <- L
y <- c(1:4)
zAIC <- t(zAIC)
zSC <- t(zSC)

myPanel <- function(x, y, z,...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,2))
}
#
StickinessAIC <- levelplot(zAIC, grid, panel=myPanel,
    col.regions = heat.colors(100)[length(heat.colors(90)):50])

StickinessSC <- levelplot(zSC, grid, panel=myPanel,
    col.regions = heat.colors(100)[length(heat.colors(90)):50])

grid.arrange(StickinessAIC,StickinessSC, nrow=2)
#
# end Stickiness
#---------------------------------------------
# begin Mean First Passage Time
#
# 2. Mean (over all models) passage time to
# the true model j, starting form model i.
# Solution to a system of L linear equations, with
# L unknowns.
# Given in the format of A*x=b
# where A is the modified transition probability matrix P with
# diagonals (p_{ii}-1), b is a (L x 1) vector of -1,
# x are tau_{i,T} of mean passage times.
#
TauAllAIC <- matrix(0, L, L)
TauTessAIC <- matrix(0, L, L)
TauBoAIC <- matrix(0, L, L)
TauMaveAIC <- matrix(0, L, L)
#
TauAllSC <- matrix(0, L, L)
TauTessSC <- matrix(0, L, L)
TauBoSC <- matrix(0, L, L)
TauMaveSC <- matrix(0, L, L)

D <- diag(rep(1, (L-1)))
b <- rep(-1, (L-1))
for (i in 1:L) {
  # By definition of passage time tau_{T,T}=0.
  # For the system of linear equations, we have (L-1) x (L-1)
  # matrix of coefficients
  vallaic <- solve((P[,,1,i][-i,-i] - D), b)
  vtessaic <- solve((P[,,2,i][-i,-i] - D), b)
  vboaic <- solve((P[,,3,i][-i,-i] - D), b)
  vmaveaic <- solve((P[,,4,i][-i,-i] - D), b)
  #
  vallsc <- solve((Q[,,1,i][-i,-i] - D), b)
  vtesssc <- solve((Q[,,2,i][-i,-i] - D), b)
  vbosc <- solve((Q[,,3,i][-i,-i] - D), b)
  vmavesc <- solve((Q[,,4,i][-i,-i] - D), b)
  #
  if (i == 1){
    TauAllAIC[i,] <- c(0, vallaic)
    TauTessAIC[i,] <- c(0, vtessaic)
    TauBoAIC[i,] <- c(0, vboaic)
    TauMaveAIC[i,] <- c(0, vmaveaic)
  } else if(i == 14) {
    TauAllAIC[i,] <- c(vallaic, 0)
    TauTessAIC[i,] <- c(vtessaic, 0)
    TauBoAIC[i,] <- c(vboaic, 0)
    TauMaveAIC[i,] <-  c(vmaveaic, 0)
  } else {
    TauAllAIC[i,] <- c(vallaic[1:(i-1)], 0, vallaic[i:(L-1)])
    TauTessAIC[i,] <- c(vtessaic[1:(i-1)], 0, vtessaic[i:(L-1)])
    TauBoAIC[i,] <- c(vboaic[1:(i-1)], 0, vboaic[i:(L-1)])
    TauMaveAIC[i,] <- c(vmaveaic[1:(i-1)], 0, vmaveaic[i:(L-1)])
  }
  #
  if (i == 1) {
    TauAllSC[i,] <- c(0, vallsc)
    TauTessSC[i,] <- c(0, vtesssc)
    TauBoSC[i,] <- c(0, vbosc)
    TauMaveSC[i,] <-  c(0, vmavesc)
  } else if( i==14 ) {
    TauAllSC[i,] <- c(vallsc, 0)
    TauTessSC[i,] <- c(vtesssc, 0)
    TauBoSC[i,] <- c(vbosc, 0)
    TauMaveSC[i,] <-  c(vmavesc, 0)
  } else {
    TauAllSC[i,] <- c(vallsc[1:(i-1)], 0, vallsc[i:(L-1)])
    TauTessSC[i,] <- c(vtesssc[1:(i-1)], 0, vtesssc[i:(L-1)])
    TauBoSC[i,] <- c(vbosc[1:(i-1)], 0, vbosc[i:(L-1)])
    TauMaveSC[i,] <- c(vmavesc[1:(i-1)], 0, vmavesc[i:(L-1)])
  }
}
#-----------------------------------------
# Produce Plots
#
# Model 10 is the most complex, switch with model 14
# to put in correct order of complexity
#
#
# Plot Results
modelnum <- c(1:L)

TauAllAIC <- reAssign(TauAllAIC)
TauTessAIC <- reAssign(TauTessAIC)
TauMaveAIC <- reAssign(TauMaveAIC)
TauBoAIC <- reAssign(TauBoAIC)
#
TauAllSC <- reAssign(TauAllSC)
TauTessSC <- reAssign(TauTessSC)
TauMaveSC <- reAssign(TauMaveSC)
TauBoSC <- reAssign(TauBoSC)

MFPTAllAIC.dat <- as.data.frame(TauAllAIC)
MFPTTessAIC.dat <- as.data.frame(TauTessAIC)
MFPTBoAIC.dat <- as.data.frame(TauBoAIC)
MFPTMaveAIC.dat <- as.data.frame(TauMaveAIC)
#
MFPTAllSC.dat <- as.data.frame(TauAllSC)
MFPTTessSC.dat <- as.data.frame(TauTessSC)
MFPTBoSC.dat <- as.data.frame(TauBoSC)
MFPTMaveSC.dat <- as.data.frame(TauMaveSC)

# Heatmap
maxcol <- max(max(TauAllAIC), max(TauTessAIC), max(TauBoAIC), max(TauMaveAIC),
             max(TauAllSC), max(TauTessSC), max(TauBoSC), max(TauMaveSC)) + 0.2
mincol <- min(min(TauAllAIC), min(TauTessAIC), min(TauBoAIC), min(TauMaveAIC),
             min(TauAllSC), min(TauTessSC), min(TauBoSC), min(TauMaveSC)) + 0.2

PlotMFTPTessAIC <- levelplot(t(MFPTTessAIC.dat),
          xlab = "To True Model (AIC)",
          ylab = "From Model",
          at = seq(mincol,maxcol,0.5),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Tess")

PlotMFTPMaveAIC <- levelplot(t(MFPTMaveAIC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,0.25),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Mave")

PlotMFTPBoAIC <- levelplot(t(MFPTBoAIC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,0.2),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Bo")

PlotMFTPAllAIC <- levelplot(t(MFPTAllAIC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,0.2),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "All")

PlotMFTPTessSC <- levelplot(t(MFPTTessSC.dat),
          xlab = "To True Model (SC)",
          ylab = "From Model",
          at = seq(mincol,maxcol,1),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Tess")

PlotMFTPMaveSC <- levelplot(t(MFPTMaveSC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,1),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Mave")

PlotMFTPBoSC <- levelplot(t(MFPTBoSC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,1),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "Bo")

PlotMFTPAllSC <- levelplot(t(MFPTAllSC.dat),
          xlab = "To True Model",
          ylab = "From Model",
          at = seq(mincol,maxcol,1),
          col.regions = heat.colors(100)[length(heat.colors(100)):10],
          main = "All")

p <- grid.arrange(PlotMFTPTessAIC, PlotMFTPMaveAIC, PlotMFTPBoAIC,
    PlotMFTPAllAIC, PlotMFTPTessSC, PlotMFTPMaveSC, PlotMFTPBoSC,
    PlotMFTPAllSC, ncol = 4, nrow = 2)

ggsave(filename = paste0(outputDir, "/Fig6.pdf"), plot = p,
    width = 15, height = 8, units = "in")

mean(TauTessAIC)
sd(TauTessAIC)
mean(TauMaveAIC)
sd(TauMaveAIC)
mean(TauBoAIC)
sd(TauBoAIC)
mean(TauAllAIC)
sd(TauAllAIC)
mean(TauTessSC)
sd(TauTessSC)
mean(TauMaveSC)
sd(TauMaveSC)
mean(TauBoSC)
sd(TauBoSC)
mean(TauAllSC)
sd(TauAllSC)

# end Mean First Passage Time
#------------------------------------------------------------
#------------------------------------------------------------
# B. Stationary Properties of the System
#------------------------------------------------------------
# begin Stationary Distribution
#

StatDistModelsSCAll <- matrix(0, L, L) # All equal proportion of scientists in the population, using SC as model selection statistic
StatDistModelsSCTess <- matrix(0, L, L)
StatDistModelsSCBo <- matrix(0, L, L)
StatDistModelsSCMave <- matrix(0, L, L)
#
StatDistModelsAICAll <- matrix(0, L, L)
StatDistModelsAICTess <- matrix(0, L, L)
StatDistModelsAICBo <- matrix(0, L, L)
StatDistModelsAICMave <- matrix(0, L, L)
#----------------------------------------------------------
# Get the transition matrix for the Markov chain

for (i in 1:L){

# Normalize to avoid numerical error (probability distribution sums to 1)

PAll <- P[,,1,i]; PTess <- P[,,2,i]; PBo <- P[,,3,i]; PMave <- P[,,4,i];

QAll <- Q[,,1,i]; QTess <- Q[,,2,i]; QBo <- Q[,,3,i]; QMave <- Q[,,4,i];

PPAll <- PAll%*%PAll; PPTess <- PTess%*%PTess; PPBo <- PBo%*%PBo; PPMave <- PMave%*%PMave;

QQAll <- QAll%*%QAll; QQTess <- QTess%*%QTess; QQBo <- QBo%*%QBo; QQMave <- QMave%*%QMave;

for (j in 1:30) {
# run for 30 time steps
# stationarity is reached quickly so we do not implement
# a numerical convergence check, but easy to run for each
# time step and check for numerical convergence if need be
  PPAll <- PPAll%*%PAll; PPTess <- PPTess%*%PTess;  PPBo <- PPBo%*%PBo;  PPMave <- PPMave%*%PMave;

  QQAll <- QQAll%*%QAll; QQTess <- QQTess%*%QTess;  QQBo <- QQBo%*%QBo;  QQMave <- QQMave%*%QMave;
}

StatDistModelsAICAll[i,] <- PPAll[1,]
StatDistModelsAICTess[i,] <- PPTess[1,]
StatDistModelsAICBo[i,] <- PPBo[1,]
StatDistModelsAICMave[i,] <- PPMave[1,]


StatDistModelsSCAll[i,] <- QQAll[1,]
StatDistModelsSCTess[i,] <- QQTess[1,]
StatDistModelsSCBo[i,] <- QQBo[1,]
StatDistModelsSCMave[i,] <- QQMave[1,]
}

StatDistModelsAICAll <- reAssign(StatDistModelsAICAll)
StatDistModelsAICTess <- reAssign(StatDistModelsAICTess)
StatDistModelsAICBo <- reAssign(StatDistModelsAICBo)
StatDistModelsAICMave <- reAssign(StatDistModelsAICMave)

StatDistModelsSCAll <- reAssign(StatDistModelsSCAll)
StatDistModelsSCTess <- reAssign(StatDistModelsSCTess)
StatDistModelsSCBo <- reAssign(StatDistModelsSCBo)
StatDistModelsSCMave <- reAssign(StatDistModelsSCMave)

#------------------------------------------------------

# For Plots create summary information:

# For three principal probabilities for each true model,
# and scientist population, create a matrix of model IDs (cols 1:3)
# and stationary probabilities for these models (cols 4:6) corresponding
# to these models
#
nclosemodels <- 3
AICAll <- matrix(0, L, nclosemodels * 2)
AICTess <- AICAll
AICBo <- AICAll
AICMave <- AICAll
#
SCAll <- matrix(0, L, nclosemodels * 2)
SCTess <- SCAll
SCBo <- SCAll
SCMave <- SCAll

num <- c(1:L)
for (i in 1: L){
  xx <- rbind(num, StatDistModelsAICAll[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  AICAll[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsAICTess[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  AICTess[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsAICBo[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  AICBo[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsAICMave[i,])
  temp  <-  xx[,order(-xx[nrow(xx),])]
  AICMave[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #------------------------------------------
  xx <- rbind(num, StatDistModelsSCAll[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  SCAll[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsSCTess[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  SCTess[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsSCBo[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  SCBo[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
  xx <- rbind(num, StatDistModelsSCMave[i,])
  temp  <-  xx[, order(-xx[nrow(xx),])]
  SCMave[i,] <- cbind(temp[1, 1:3],temp[2, 1:3])
  #
}


dfbestthree <- data.frame(cbind(
        AICTess[,4:6], AICMave[,4:6], AICBo[,4:6], AICAll[,4:6],
        SCTess[,4:6], SCMave[,4:6], SCBo[,4:6], SCAll[,4:6],
        AICTess[,1:3], AICMave[,1:3], AICBo[,1:3], AICAll[,1:3],
        SCTess[,1:3], SCMave[,1:3], SCBo[,1:3], SCAll[,1:3]))
colnames(dfbestthree) <- c(
    rep("P(AICTess)", 3), rep("P(AICMave)", 3), rep("P(AICBo)", 3), rep("P(AICAll)", 3),
    rep("P(SCTess)", 3), rep("P(SCMave)", 3), rep("P(SCBo)", 3), rep("P(SCAll)", 3),
    rep("M(AICTess)", 3), rep("M(AICMave)", 3), rep("M(AICBo)", 3), rep("M(AICAll)", 3),
    rep("M(SCTess)", 3), rep("M(SCMave)", 3), rep("M(SCBo)", 3), rep("M(SCAll)", 3))
#
# end Stationary Probability Distribution
#--------------------------------------------------------------
#--------------------------------------------------------------
# end all
#--------------------------------------------------------------