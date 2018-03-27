################
##
## @description Main reproducibility simulation
##
## @param None
##
## @return None
##
## @lastChange 2018-03-22
##
## @changes
##  Terminology adjustment [2018-03-22]
##  Parameters output file [2017-06-17]
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
scriptDir <- paste0(baseDir, "/src/abm/functions")
inputDir <- paste0(baseDir, "/data")
outputDir <- paste0(baseDir, "/data")


#############
## FUNCTIONS
#############
source(paste0(scriptDir, "/analysis.R"))
source(paste0(scriptDir, "/calculateDet.R"))
source(paste0(scriptDir, "/calculateDistance.R"))
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/generateBetas.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/generateXSet.R"))
source(paste0(scriptDir, "/generateY.R"))
source(paste0(scriptDir, "/getBetas.R"))
source(paste0(scriptDir, "/getPredictors.R"))
source(paste0(scriptDir, "/modelSimilarByInteraction.R"))
source(paste0(scriptDir, "/modelSimilarByTerm.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/seedGenerator.R"))
source(paste0(scriptDir, "/simulator.R"))
source(paste0(scriptDir, "/strToModel.R"))


###################
## INPUT PARAMETERS
###################
## Number of replications
replications <- 100

## Length of the simulation
timesteps <- 11000

## Number of factors
k <- 3

## Sigma (Error variance)
sigma <- 0.2

## Sample size
sampleSize <- 100

##
## Generate all possible linear regression models with k number of factors
##
## The linear regression models are represented in a matrix where columns
## represent the factors and each row represent the terms. Cells with value 1
## indicate the factor in each term.
##
## Example: The following matrix represents the linear regression model
##          x1 + x2 + x3 + x1x2 + x1x3
##
##      1   0   0
##      0   1   0
##      0   0   1
##      1   1   0
##      1   0   1
##
models <- generateModels(k)

## Identify the number of predictors of each model
predictors <- getPredictors(models)

## Generate Betas
weights <- generateBetas(models)

## True model
trueModel <- "x1 + x2"
tModel <- strToModel(trueModel, k)

## Correlation
correlation <- 0.2

## Number of replicators
nRey <- 1

## Number of theory testers (i.e., add or remove term)
nTess <- 1

## Number of boundary testers (i.e., add interaction)
nBo <- 1

## Number of mavericks (random select a model)
nMave <- 1

##
## Model comparison
##
## TSTATISTICS    T Statistic
## RSQ            R-Squared
## ARSQ           Adjusted R-Squared
## AIC            Akaike Information Criterion
## BIC            Bayesian Information Criterion
##
modelCompare <- AIC

## Output filename
outputFile <- "output.csv"
paramFile <- "parameters.rds"

## Verbose mode
verbose <- TRUE

## Number of decimal places
ndec <- 4


###################
## SET SEED
###################
seeds <- seedGenerator(replications, paste0(inputDir, "/seeds.csv"))


###################
## SIMULATION
###################
simulator(replications, timesteps, models, k, tModel,
    nRey, nTess, nBo, nMave, weights, sampleSize, correlation, sigma, 
    modelCompare, inputDir, outputDir, outputFile, paramFile,
    verbose, ndec, seeds)
