################
##
## @description Generate all possible models with k number of factors
##
## @param k Number of factors
##
## @return List of all possible models
##
## @lastChange 2016-12-28
##
## @changes
##
################
generateModels <- function(k){
  
  models <- NULL
  
  if(k < 1){
    return(list(c(0)))
  } else if(k == 1){
    return(list(c(1)))
  }
  
  ## Generate all possible combination of factors
  rows <- NULL
  for(i in 1:k){
    comb <- combs(1:k, i)
    
    for(j in 1:nrow(comb)){
      rows <- rbind(rows, convertBinary(comb[j,], k))
    }
  }
  
  ## Generate all possible models
  factor1 <- c(1, rep(0, k-1))
  
  pModels <- list()
  for(i in 2:nrow(rows)){
    comb <- combs(2:nrow(rows), i-1)
    
    for(j in 1:nrow(comb)){
      
      model <- matrix(0, nrow=length(comb[j,]) + 1, ncol=k)
      
      model[1,] <- factor1
      index <- 2
      for(r in comb[j,]){
        model[index,] <- rows[r,]
        index <- index + 1
      }
      
      pModels[length(pModels) + 1] <- list(model)
    }
  }
  
  ## Validate the models
  f <- 10^((k-1):0)
  models <- list(c(1, rep(0, k-1)))
  for(model in pModels){
    
    ## Get the index of single terms and identify existing interactions
    nTerms <- 0
    terms <- c()
    inters <- c()
    interTerms <- c()
    indInters <- NULL
    for(r in 1:nrow(model)){
      if(sum(model[r,]) == 1){
        nTerms <- nTerms + 1
        terms <- c(terms, which(model[r,] == 1))
      } else {
        inters <- c(inters, sum(model[r,] * f))
        indInters <- rbind(indInters, model[r,])
        aux <- which(model[r,] == 1)
        for(z in aux){
          if(!(z %in% interTerms)){
            interTerms <- c(interTerms, z)
          }
        }
      }
    }
    
    ## Check if factors in the interaction exist in isolation
    if(is.null(inters)){
      if(searchModel(model, models) == 0){
        models[length(models) + 1] <- list(model)
      }
    } else if((setequal(interTerms, terms)) ||
              (all(is.element(setdiff(interTerms, terms), terms)))){
      
      ## Add all necessary interactions
      newModel <- model
      if(length(interTerms) > 1){
        for(inter in 1:nrow(indInters)){
          indexInter <- which(indInters[inter,] == 1)
          for(i in 2:length(indexInter)){
            comb <- combs(indexInter, i)
            for(c in 1:nrow(comb)){
              aux <- rep(0, k)
              aux[comb[c,]] <- 1
              if(!(sum((aux) * f) %in% inters)){
                newModel <- rbind(newModel, aux)
                inters <- c(inters, sum((aux) * f))
              }
            }
          }
        }
      }
      rownames(newModel) <- c()
      
      if(searchModel(newModel, models) == 0){
        models[length(models) + 1] <- list(newModel)
      }
    }
  }
  
  return(models)
}
