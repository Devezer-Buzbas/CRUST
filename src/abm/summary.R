################
##
## @description Generate summary of multiple simulations
##
## @param None
##
## @return None
##
## @lastChange 2018-03-22
##
## @changes
##   Terminology adjustment [2018-03-22]
##   Add frequency of replicated No True Model as Global Models [2018-03-06]
##   Add frequency of contiguous sequence sizes [2018-03-006]
##   If no REY, set replicated, replicatedTGM and replicatedNTGM to NA [2018-02-26]
##   Fixed replicatedTGM and replicatedNTGM [2017-06-18]
##   Add Beta Bias [2017-06-18]
##   Fixed changesBeforeFirstTGM [2017-06-18]
##
################
library(caTools)
library(data.table)


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
source(paste0(scriptDir, "/compareModels.R"))
source(paste0(scriptDir, "/constants.R"))
source(paste0(scriptDir, "/convertBinary.R"))
source(paste0(scriptDir, "/generateModels.R"))
source(paste0(scriptDir, "/modelToStr.R"))
source(paste0(scriptDir, "/searchModel.R"))
source(paste0(scriptDir, "/strToModel.R"))


###################
## INPUT PARAMETERS
###################
## Replications
replications <- 100

## Length of the simulation
timesteps <- 11000

## Number of timesteps to skip
skip <- 1000

## Number of factors
k <- 3

## Generate all models
models <- generateModels(k)

## Number of models
numModels <- length(models)

m <- c(2, 7, 14)

sigmas <- 1:3

types <- 1:5

verbose <- TRUE

###################
## SUMMARY
###################
output <- NULL
outputCode <- NULL
frequency <- NULL
frequencyCode <- NULL
freqCont <- NULL
freqContCode <- NULL
freqNTGM <- NULL
freqNTGMCode <- NULL
for(tModelIndex in 1:length(m)){
  for(sigmaIndex in sigmas){
    for(typeIndex in types){
      
      if(verbose){
        print(paste0(tModelIndex, "-", sigmaIndex, "-", typeIndex))
      }
      
      ## True model
      tModel <- modelToStr(models[[m[tModelIndex]]])
      
      ## Sigma
      if(sigmaIndex == 1){
        sigma <- 0.2
      } else if(sigmaIndex == 2){
        sigma <- 0.5
      } else if(sigmaIndex == 3){
        sigma <- 0.8
      }
      
      ## Agent types
      if(typeIndex == 1){
        nRey <- 1
        nTess <- 1
        nBo <- 1
        nMave <- 1
      } else if(typeIndex == 2){
        nRey <- 300
        nTess <- 1
        nBo <- 1
        nMave <- 1
      } else if(typeIndex == 3){
        nRey <- 1
        nTess <- 300
        nBo <- 1
        nMave <- 1
      } else if(typeIndex == 4){
        nRey <- 1
        nTess <- 1
        nBo <- 300
        nMave <- 1
      } else if(typeIndex == 5){
        nRey <- 1
        nTess <- 1
        nBo <- 1
        nMave <- 300
      }
      
      ## Upload simulation data file
      filename <- paste0("output-", tModelIndex, "-", sigmaIndex, "-",
          typeIndex, ".csv")
      data <- fread(paste0(inputDir, "/", filename), sep=";")
      
      data <- data[(skip+1):nrow(data),]
      
      ## Frequency of each model
      freq <- data[which(timestep <= timesteps), list(num=.N),
          by=list(replica, final_global_model)]
      
      for(r in 1:replications){
        aux <- setdiff(1:numModels,
            unique(freq[which(replica == r)]$final_global_model))
        if(length(aux) > 0){
          freq <- rbind(freq, as.data.table(cbind(replica=r,
                      final_global_model=aux, num=rep(0, length(aux)))))
        }
      }
      freq <- freq[order(replica, final_global_model)]
      
      ## Proportion of times True Model was selected as Global Model
      selTrueModel <- data[which((final_global_true_model == 1) &
                  (timestep <= timesteps)), list(nrows=.N / timesteps),
              by=replica]
      
      aux <- setdiff(1:replications, unique(selTrueModel$replica))
      if(length(aux) > 0){
        selTrueModel <- rbind(selTrueModel, as.data.table(cbind(replica=aux,
                    nrows=rep(0, length(aux)))))
      }
      selTrueModel <- selTrueModel[order(replica)]
      
      ## Proportion of times the Global Model departs from the True Model
      ## when the Global Model was the True Model
      numerator <- data[which((initial_global_true_model == 1) &
                  (final_global_true_model == 0) &
                  (timestep <= timesteps)), list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(numerator$replica))
      if(length(aux) > 0){
        numerator <- rbind(numerator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      numerator <- numerator[order(replica)]
      
      denominator <- data[which((initial_global_true_model == 1) &
                  (timestep <= timesteps)), list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(denominator$replica))
      if(length(aux) > 0){
        denominator <- rbind(denominator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      denominator <- denominator[order(replica)]
      
      departTrueModel <- as.data.table(cbind(replica=numerator$replica,
          nrows=ifelse(denominator$nrows == 0, NA,
              1 - (numerator$nrows / denominator$nrows))))
      
      departTrueModel <- departTrueModel[order(replica)]
      
      ## Proportion of times that REY reproduced previous result
      replicated <- data[which((strategy == REY) &
                  (timestep <= timesteps)),
          list(num=mean(replicated)), by=replica]
      
      aux <- setdiff(1:replications, unique(replicated$replica))
      if(length(aux) > 0){
        replicated <- rbind(replicated, as.data.table(cbind(
                    replica=aux, num=rep(NA, length(aux)))))
      }
      replicated <- replicated[order(replica)]
      
      ## First time the Global Model is the True Model
      firstTGM <- data[which((initial_global_true_model == 0) &
                  (final_global_true_model == 1) &
                  (timestep <= timesteps)), .SD, by=replica][,
          list(mints=min(timestep)), by=replica]
      
      aux <- setdiff(1:replications, unique(firstTGM$replica))
      if(length(aux) > 0){
        firstTGM <- rbind(firstTGM, as.data.table(cbind(replica=aux,
                    mints=rep(NA, length(aux)))))
        firstTGM <- firstTGM[order(replica)]
      }
      
      ## Number of switches
      numSwitches <- data[which((((initial_global_true_model == 0) &
                      (final_global_true_model == 1)) |
                    ((initial_global_true_model == 1) &
                      (final_global_true_model == 0))) &
                  (timestep <= timesteps)), list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(numSwitches$replica))
      if(length(aux) > 0){
        numSwitches <- rbind(numSwitches, as.data.table(cbind(replica=aux,
                    nrows=rep(NA, length(aux)))))
      }
      numSwitches <- numSwitches[order(replica)]
      
      ## Average and Standard Deviation of contiguous time steps with the
      ## Global Model as True Model
      auxCont <- vector()
      aux <- data[which(timestep <= timesteps),
          list(final_global_true_model), by=replica]
      runLenTGM <- NULL
      for(i in unique(aux$replica)){
        l <- rle(aux[which(replica == i),]$final_global_true_model)
        runLenTGM <- rbind(runLenTGM, cbind(replica=i,
                avg=ifelse(length(l$lengths[l$values == 1]) > 0,
                    mean(l$lengths[l$values == 1]),
                    NA),
                sd=ifelse(length(l$lengths[l$values == 1]) > 1,
                    sd(l$lengths[l$values == 1]),
                    NA)))
        
        ## Calculate the frequency of continuous blocks
        for(j in l$lengths[l$values == 1]){
          if(is.na(auxCont[j])){
            auxCont[j] <- 1
          } else {
            auxCont[j] <- auxCont[j] + 1
          }
        }
      }
      runLenTGM <- as.data.table(runLenTGM)
      
      aux <- setdiff(1:replications, unique(runLenTGM$replica))
      if(length(aux) > 0){
        runLenTGM <- rbind(runLenTGM, as.data.table(cbind(replica=aux,
                    avg=rep(NA, length(aux)),
                    sd=rep(NA, length(aux)))))
      }
      runLenTGM <- runLenTGM[order(replica)]
      
      ## Number of changes before the True Model is selected as
      ## Global Model for the first time
      changesBeforeFirstTGM <- NULL
      for(r in 1:replications){
        changesBeforeFirstTGM <- rbind(changesBeforeFirstTGM,
            cbind(r, data[
            which((initial_global_model != final_global_model) &
                    (replica == r) &
                    (timestep <= firstTGM[which(replica == r),]$mints)),
            .N]))
      }
      changesBeforeFirstTGM <- data.table(changesBeforeFirstTGM)
      names(changesBeforeFirstTGM) <- c("replica", "num")
      
      aux <- setdiff(1:replications, unique(changesBeforeFirstTGM$replica))
      if(length(aux) > 0){
        changesBeforeFirstTGM <- rbind(changesBeforeFirstTGM,
            as.data.table(cbind(replica=aux, num=rep(NA, length(aux)))))
      }
      changesBeforeFirstTGM <- changesBeforeFirstTGM[order(replica)]
      
      ## Level of reproducibility when the True Model is the Global Model
      numerator <- data[which((initial_global_true_model == 1) &
                  (strategy == REY) & (replicated == 1) &
                  (timestep <= timesteps)),
          list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(numerator$replica))
      if(length(aux) > 0){
        numerator <- rbind(numerator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      numerator <- numerator[order(replica)]
      
      denominator <- data[which((initial_global_true_model == 1) &
                  (strategy == REY) &
                  (timestep <= timesteps)), list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(denominator$replica))
      if(length(aux) > 0){
        denominator <- rbind(denominator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      denominator <- denominator[order(replica)]
      
      replicatedTGM <- as.data.table(cbind(replica=numerator$replica,
              nrows=ifelse(denominator$nrows > 0,
                  numerator$nrows / denominator$nrows, NA)))
      
      aux <- setdiff(1:replications, unique(replicatedTGM$replica))
      if(length(aux) > 0){
        replicatedTGM <- rbind(replicatedTGM,
            as.data.table(cbind(replica=aux, nrows=rep(NA, length(aux)))))
      }
      replicatedTGM <- replicatedTGM[order(replica)]
      
      ## Level of reproducibility when the True Model is not
      ## the Global Model
      numerator <- data[which((initial_global_true_model == 0) &
                  (strategy == REY) & (replicated == 1) &
                  (timestep <= timesteps)),
          list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(numerator$replica))
      if(length(aux) > 0){
        numerator <- rbind(numerator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      numerator <- numerator[order(replica)]
      
      denominator <- data[which((initial_global_true_model == 0) &
                  (strategy == REY) &
                  (timestep <= timesteps)), list(nrows=.N), by=replica]
      
      aux <- setdiff(1:replications, unique(denominator$replica))
      if(length(aux) > 0){
        denominator <- rbind(denominator, as.data.table(cbind(
                    replica=aux, nrows=rep(0, length(aux)))))
      }
      denominator <- denominator[order(replica)]
      
      replicatedNTGM <- as.data.table(cbind(replica=numerator$replica,
              nrows=ifelse(denominator$nrows > 0,
                  numerator$nrows / denominator$nrows, NA)))
      
      aux <- setdiff(1:replications, unique(replicatedNTGM$replica))
      if(length(aux) > 0){
        replicatedNTGM <- rbind(replicatedNTGM,
            as.data.table(cbind(replica=aux, nrows=rep(NA, length(aux)))))
      }
      replicatedNTGM <- replicatedNTGM[order(replica)]
      
      ## Frequency of reproducibility of each model when the model is
      ## not the True Model
      frequencyNTGM <- data[which((initial_global_true_model == 0) &
                  (strategy == REY) & (replicated == 1) &
                  (initial_global_model != tModelIndex) &
                  (timestep <= timesteps)),
          list(nrows=.N), by=initial_global_model]
      
      frequencyNTGM <- frequencyNTGM[order(initial_global_model)]
      
      ## Bias Betas
      biasBeta <- data[which(timestep <= timesteps),
                       list(avg=mean(beta1_true - beta1_estimate)),
                       by=replica]
      
      aux <- setdiff(1:replications, unique(biasBeta$replica))
      if(length(aux) > 0){
        biasBeta <- rbind(biasBeta, as.data.table(cbind(
          replica=aux, nrows=rep(NA, length(aux)))))
        biasBeta <- biasBeta[order(replica)]
      }
      
      ## Output
      output <- rbind(output, cbind(1:replications,
              tModel, sigma, nRey, nTess, nBo, nMave,
              selTrueModel$nrows, departTrueModel$nrows, replicated$num,
              firstTGM$mints, numSwitches$nrows, runLenTGM$avg,
              runLenTGM$sd, changesBeforeFirstTGM$num,
              replicatedTGM$nrows, replicatedNTGM$nrows, biasBeta$avg))
      
      outputCode <- rbind(outputCode, cbind(1:replications, tModelIndex,
              sigmaIndex, typeIndex,
              selTrueModel$nrows, departTrueModel$nrows, replicated$num,
              firstTGM$mints, numSwitches$nrows, runLenTGM$avg,
              runLenTGM$sd, changesBeforeFirstTGM$num,
              replicatedTGM$nrows, replicatedNTGM$nrows, biasBeta$avg))
      
      for(r in 1:replications){
        frequency <- rbind(frequency, cbind(r, tModel, sigma,
                nRey, nTess, nBo, nMave,
                array(freq$final_global_model, dim=c(length(models),
                        replications))[,r],
                array(freq$num, dim=c(length(models), replications))[,r]))
        
        frequencyCode <- rbind(frequencyCode, cbind(r, tModelIndex,
                sigmaIndex, typeIndex,
                array(freq$final_global_model, dim=c(length(models),
                        replications))[,r],
                array(freq$num, dim=c(length(models), replications))[,r]))
      }
      
      for(j in 1:length(auxCont)){
        freqCont <- rbind(freqCont, cbind(tModel, sigma,
                nRey, nTess, nBo, nMave, j, auxCont[j]))
        
        freqContCode <- rbind(freqContCode, cbind(tModelIndex,
                sigmaIndex, typeIndex, j, auxCont[j]))
      }
      
      for(j in 1:nrow(frequencyNTGM)){
        freqNTGM <- rbind(freqNTGM, cbind(tModel, sigma,
                nRey, nTess, nBo, nMave,
                as.integer(frequencyNTGM[j, 1]),
                as.integer(frequencyNTGM[j, 2])))
        
        freqNTGMCode <- rbind(freqNTGMCode, cbind(tModelIndex,
                sigmaIndex, typeIndex,
                as.integer(frequencyNTGM[j, 1]),
                as.integer(frequencyNTGM[j, 2])))
      }
    }
  }
}

## Summary data
output <- data.table(output)
names(output) <- c("replica", "tModel", "sigmaIndex",
    "nRey", "nTess", "nBo", "nMave",
    "selTrueModel", "departTrueModel", "replicated",
    "firstTGM", "numSwitches", "avgContTGM", "sdContTGM",
    "changesBeforeFirstTGM", "replicatedTGM", "replicatedNTGM", "biasBeta1")

write.table(output, file=paste0(outputDir, "/summary.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

outputCode <- data.table(outputCode)
names(outputCode) <- c("replica", "tModelIndex", "sigmaIndex", "typeIndex",
    "selTrueModel", "departTrueModel", "replicated",
    "firstTGM", "numSwitches", "avgContTGM", "sdContTGM",
    "changesBeforeFirstTGM", "replicatedTGM", "replicatedNTGM", "biasBeta1")

write.table(outputCode, file=paste0(outputDir, "/summaryCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

## Frequency of the models
frequency <- data.table(frequency)
names(frequency) <- c("replica", "tModel", "sigmaIndex",
    "nRey", "nTess", "nBo", "nMave",
    "final_global_model", "frequency")

write.table(frequency, file=paste0(outputDir, "/frequency.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

frequencyCode <- data.table(frequencyCode)
names(frequencyCode) <- c("replica", "tModelIndex", "sigmaIndex", "typeIndex",
    "final_global_model", "frequency")

write.table(frequencyCode, file=paste0(outputDir, "/frequencyCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

## Frequency contiguity lentgh
freqCont <- data.table(freqCont)
names(freqCont) <- c("tModel", "sigmaIndex", "nRey", "nTess", "nBo", "nMave",
    "contiguityLength", "frequency")

write.table(freqCont, file=paste0(outputDir, "/freqCont.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

freqContCode <- data.table(freqContCode)
names(freqContCode) <- c("tModelIndex", "sigmaIndex", "typeIndex",
    "contiguityLength", "frequency")

write.table(freqContCode, file=paste0(outputDir, "/freqContCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

## Frequency models replicated NTGM
freqNTGM <- data.table(freqNTGM)
names(freqNTGM) <- c("tModel", "sigmaIndex", "nRey", "nTess", "nBo", "nMave",
    "model", "frequency")

write.table(freqNTGM, file=paste0(outputDir, "/freqNTGM.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")

freqNTGMCode <- data.table(freqNTGMCode)
names(freqNTGMCode) <- c("tModelIndex", "sigmaIndex", "typeIndex",
    "model", "frequency")

write.table(freqNTGMCode, file=paste0(outputDir, "/freqNTGMCode.csv"),
    append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=";")
