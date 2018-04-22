################
##
## @description Upload seeds from a file or generate randomly
##
## @param N        Number of seeds to upload or generate
## @param filename Name of the filename
##
## @return Vector of seeds
##
## @lastChange 2016-12-28
##
## @changes
##
################
seedGenerator <- function(N, filename = "seeds.csv"){
  
  if(file.exists(filename)){
    s <- fread(filename, sep=";", header=FALSE)
    seeds <- array(0, dim=nrow(s))
    for(i in 1:nrow(s)){
      seeds[i] <- s[i,]
    }
  } else {
    seeds <- array(0, dim=N)
    s <- 1
    while(s <= N){
      seed <- as.integer(runif(1, min=0, max=(N*10)))
      seeds[s] <- seed
      s <- s + 1
    }
    
    s <- data.table(seeds)
    fname <- write.table(s, file=filename, quote=FALSE, sep=";",
        row.names=FALSE, col.names=FALSE)
  }
  
  return(unlist(seeds))
}
