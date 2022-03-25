library(dplyr)
library(tidyverse)

pathToFolder <- "outputTabs/iPSC/outRDS/"
df <- paste0(pathToFolder,list.files(path=pathToFolder,pattern = ".RDS")) %>%
  map(readRDS) 
names(df) <- sapply(strsplit(list.files(path=pathToFolder,pattern = ".RDS"),"\\."), function(x) x[2])

stopifnot(all(sapply(df, function(x) dim(x)[1])==19653))


matrixBuilder <- function(varCategory="All"){
  
  varCat <- sapply(df, function(x) {
    
    vec <- x[,varCategory]
    names(vec) <- x[,"geneSymbol"]
    vec
    
  }, simplify=F)
  
  varCat <- do.call("cbind", varCat)
  return(varCat)
  
}

mutBurden_All <- matrixBuilder("All")
mutBurden_missPatho <- matrixBuilder("missPatho")
mutBurden_ptv <- matrixBuilder("ptv")
mutBurden_del <- matrixBuilder("del")
mutBurden_missnonPatho <- matrixBuilder("missense")
mutBurden_synonymous <- matrixBuilder("synonymous")
mutBurden_other <- matrixBuilder("other")

stopifnot(sum(sapply(df, function(x) sum(x["All"])))==sum(rowSums(mutBurden_All)))


newPathToFolder <- "outputTabs/iPSC/"
saveRDS(mutBurden_All,
        file=paste0(newPathToFolder, "mutBurden_All.RDS"))
saveRDS(mutBurden_missPatho,
        file=paste0(newPathToFolder, "mutBurden_missPatho.RDS"))
saveRDS(mutBurden_ptv,
        file=paste0(newPathToFolder, "mutBurden_ptv.RDS"))
saveRDS(mutBurden_del,
        file=paste0(newPathToFolder, "mutBurden_del.RDS"))
saveRDS(mutBurden_missnonPatho,
        file=paste0(newPathToFolder, "mutBurden_missnonPatho.RDS"))
saveRDS(mutBurden_synonymous,
        file=paste0(newPathToFolder, "mutBurden_synonymous.RDS"))
saveRDS(mutBurden_other,
        file=paste0(newPathToFolder, "mutBurden_other.RDS"))

