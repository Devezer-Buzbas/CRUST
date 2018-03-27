################
##
## @description Generate similar model(s) adding an interaction
##
## @param model  Model in matrix format
## @param models List of all possible models
## @param mode   "all" returns all possible similar models
##               "random" returns one similar model
##
## @return Similar model(s) adding an interaction
##
## @lastChange 2017-04-25
##
## @changes
##
################
modelSimilarByInteraction <- function(model, models, mode=c("all", "random")){
  
  mode <- match.arg(mode)
  
  if (mode == "all"){
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
    
    # Generate all possible interactions
    inter <- NULL
    if(length(terms) > 1){
      for(i in 2:length(terms)){
        comb <- combs(terms, i)
        for(c in 1:nrow(comb)){
          aux <- rep(0, k)
          aux[comb[c,]] <- 1
          if (!(sum((aux) * f) %in% inters)){
            inter <- rbind(inter, aux)
          }
        }
      }
    }
    
    # Select an interaction randomly
    similarModels <- list()
    if(!is.null(nrow(inter))){
      rownames(inter) <- c()
      
      for(r in 1:nrow(inter)){
        newModel <- similarModel
        addInter <- inter[r,]
        # Identify the index of the terms in the interaction
        terms <- which(addInter == 1)
        
        # Add all necessary interactions
        newModel <- rbind(newModel, addInter)
        rownames(newModel) <- c()
        newInters <- c(inters, sum(addInter * f))
        for(i in 2:length(terms)){
          comb <- combs(terms, i)
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
  } else if (mode == "random"){
    if(!is.matrix(model)){
      similarModel <- t(as.matrix(model))
      model <- t(as.matrix(model))
    } else {
      similarModel <- model
    }
    
    k <- length(model[1,])
    
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
    
    # Generate all possible interactions
    inter <- NULL
    if(length(terms) > 1){
      for(i in 2:length(terms)){
        comb <- combs(terms, i)
        for(c in 1:nrow(comb)){
          aux <- rep(0, k)
          aux[comb[c,]] <- 1
          if (!(sum((aux) * f) %in% inters)){
            inter <- rbind(inter, aux)
          }
        }
      }
    }
    
    # Select an interaction randomly
    if(!is.null(nrow(inter))){
      rownames(inter) <- c()
      addInter <- inter[as.integer(runif(1, min=1, max=nrow(inter)+1)),]
      
      # Identify the index of the terms in the interaction
      terms <- which(addInter == 1)
      
      # Add all necessary interactions
      similarModel <- rbind(similarModel, addInter)
      rownames(similarModel) <- c()
      inters <- c(inters, sum(addInter * f))
      for(i in 2:length(terms)){
        comb <- combs(terms, i)
        for(c in 1:nrow(comb)){
          aux <- rep(0, k)
          aux[comb[c,]] <- 1
          if(!(sum((aux) * f) %in% inters)){
            similarModel <- rbind(similarModel, aux)
            rownames(similarModel) <- c()
          }
        }
      }
    }
    
    return(similarModel)
  }
}