################
##
## @description Generate similar model(s) adding an interaction
##
## @param model          Model in matrix format
## @param models         List of all possible models
## @param mode           "all" returns all possible similar models
##                       "random" returns one similar model
## @param modelSelection "soft" select model(s) from the set of all models
##                       "hard" select model(s) from the set of possible models
##
## @return Similar model(s) adding an interaction
##
## @lastChange 2017-09-18
##
## @changes
##
################
modelSimilarByInteraction <- function(model, models, mode=c("all", "random"),
                                      modelSelection=c("hard", "soft")){

  mode <- match.arg(mode)
  modelSelection <- match.arg(modelSelection)

  if(mode == "all"){
    if(!is.matrix(model)){
      similarModel <- t(as.matrix(model))
      model <- t(as.matrix(model))
    } else {
      similarModel <- model
    }

    k <- length(similarModel[1,])

    # Get the index of the single terms and
    # identify existing interactions
    f <- 10^((k-1):0)
    nTerms <- 0
    terms <- c()
    inters <- c()
    for(r in 1:nrow(model)){
      if(sum(model[r,]) == 1){
        nTerms <- nTerms + 1
        terms <- c(terms, which(model[r,] == 1))
      } else {
        inters <- c(inters, sum(model[r,] * f))
      }
    }

    # Generate set of possible restricted interactions
    inter <- NULL
    if(length(terms) > 1){
      for(i in 2:length(terms)){
        comb <- combs(terms, i)
        for(c in 1:nrow(comb)){
          aux <- rep(0, k)
          aux[comb[c,]] <- 1
          if(!(sum(aux * f) %in% inters)){
            inter <- rbind(inter, aux)
          }
        }
      }
    }

    # Generate all models
    similarModels <- list()
    if(!is.null(nrow(inter))){
      rownames(inter) <- c()

      for(r in 1:nrow(inter)){
        newModel <- similarModel
        addInter <- inter[r,]
        # Identify the index of the terms in the interaction
        newTerms <- which(addInter == 1)

        # Add all necessary terms
        for(i in 1:length(newTerms)){
          if(!(newTerms[i] %in% terms)){
            newTerm <- rep(0, k)
            newTerm[newTerms[i]] <- 1
            newModel <- rbind(newModel, newTerm)
            rownames(newModel) <- c()
          }
        }

        # Add all necessary interactions
        newModel <- rbind(newModel, addInter)
        rownames(newModel) <- c()
        newInters <- c(inters, sum(addInter * f))
        for(i in 2:length(newTerms)){
          comb <- combs(newTerms, i)
          for(c in 1:nrow(comb)){
            aux <- rep(0, k)
            aux[comb[c,]] <- 1
            if(!(sum((aux) * f) %in% newInters)){
              newModel <- rbind(newModel, aux)
              rownames(newModel) <- c()
            }
          }
        }

        similarModels[length(similarModels) + 1] <- list(newModel)
      }
    }

    result <- array(0, dim=c(1,length(models)))
    for(m in similarModels){
      mIndex <- searchModel(m, models)
      if(mIndex > 0){
        result[mIndex] <- 1
      }
    }

    return(result)
  } else if(mode == "random"){
    aux <- modelSimilarByInteraction(model, models, "all", "hard")
    nModels <- sum(aux)
    if(modelSelection == "soft"){
      aux <- aux * (1 / (nModels + 1))
      aux[which(aux == 0)] <- (1 - sum(aux)) / (length(aux) - nModels)
    } else if(modelSelection == "hard"){
      if(nModels > 0){
        aux <- aux * (1 / nModels)
      } else {
        return(model)
      }
    }

    threshold <- runif(1)
    index <- 0
    value <- 0
    repeat{
      index <- index + 1
      value <- value + aux[index]
      if(value > threshold){
        break;
      }
    }

    similarModel <- models[[index]]
    return(similarModel)
  }
}