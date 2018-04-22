################
##
## @description Theoretical reproducibility simulation
##
## @param None
##
## @return None
##
## @lastChange 2018-04-22
##
## @changes
##
################
library(caTools)
library(data.table)
library(permute)
library(matrixStats)
library(MCMCpack)


#############
## PATHS
#############
baseDir <- "."
scriptDir <- paste0(baseDir, "/src/functions")
outputDir <- paste0(baseDir, "/data/modelComparison")


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
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/getModelSelectionConstant.R"))
source(paste0(scriptDir, "/getModelComparison.R"))


###################
## INPUT PARAMETERS
###################
## Number of factors
k = 3

## Sigma (Error variance)
sigmas = c(0.2, 0.5, 0.8)

## Model Selection
mss = c("AIC", "Schwarz")

## Sample size
sampleSize = 100

## Number of independent samples on which model comparison statistic is based
nIter = 1e5

## Generate all models
models = generateModels(k)
models[[1]] = t(as.matrix(models[[1]]))
L = length(models)

## Correlation
correlation = 0.2

## Initialize
xset = generateXSet(sampleSize, k, correlation)

# Get Model Selection Constant for all Models
msConstant = getModelSelectionConstant(models, xset)

for( ms in 1:length(mss) ) {
  for( s in 1:length(sigmas) ) {
    for( i in 1:L ) {
      tModel = models[[i]] # True Model
      ProbMC = getModelComparison(xset, sampleSize, tModel, sigmas[s], models,
          nIter, mss[ms], msConstant)
      
      sigmaTxt = gsub("\\.", "", toString(sigmasTxt[s]))
      save(ProbMC, file = paste0(outputDir,
              '/ModelComparisonProbabilitiesSigma', sigmaTxt, mss[ms], i,
              '.RData', sep=""))
    }
  }
}
