
mutTab_logReg <- readRDS("outputTabs/mutTab_logReg2.RDS")
metadata_table <- read.table("inputTabs/hipsci.qc1_sample_info.20170927.tsv",
                             header=TRUE, sep="\t")
mutTab_logReg$donor <- metadata_table[match(mutTab_logReg$cline, metadata_table$name),]$donor
mutTab_logReg_sub <- subset(mutTab_logReg, !is.na(NeuroSeq_pred))
listRep <- split(mutTab_logReg_sub$cline, mutTab_logReg_sub$donor)

namedCline <- mutTab_logReg_sub$NeuroSeq_pred
names(namedCline) <- mutTab_logReg_sub$cline

listRep <- sapply(listRep, function(x){
  namedCline[x]
})

listRep <- listRep[elementNROWS(listRep)>1]

outcomePred <- sapply(listRep, function(x){
  
  tmp <- unique(x)
  
  if (length(tmp)==1){
    
    if (tmp=="Failed"){
      "Concordant_Failed"
    } else if (tmp=="Successful"){
      "Concordant_Successful"
    }
    
  } else if (length(tmp)==2){
    "Discordant"
  }
  
})

##########################
##########################

listDisc <- listRep[names(listRep) %in% names(outcomePred[outcomePred=="Discordant"])]
#first failed, second successful replicate line (from the same donor)

listDisc[elementNROWS(listDisc)>2]

set.seed(123)
listDisc[elementNROWS(listDisc)>2] <- sapply(listDisc[elementNROWS(listDisc)>2], function(x){
  c(sample(x[x %in% "Failed"],1), sample(x[x %in% "Successful"],1))
}, simplify=F)


listDisc <- sapply(listDisc, function(x) sort(x), simplify=F)

saveRDS(listDisc, "outputTabs/iPSC/discordant/summTabs.RDS")

discordantPairsTab <- do.call("rbind", sapply(listDisc, function(x) names(x), simplify=F))
rownames(discordantPairsTab) <- NULL

filesVcf <- read.table("listOFfilesVEP.txt")

discordantPairsTab <- sapply(1:nrow(discordantPairsTab), function(x){
  
  tmp <- discordantPairsTab[x,]
  tmpFileName <- gsub(".wes.+", "", gsub("pathToWESdataVEPannot/","",filesVcf$V1))
  tmp <- t(as.data.frame(unname(sapply(tmp, function(y) filesVcf$V1[tmpFileName==y]))))
  colnames(tmp) <- NULL
  rownames(tmp) <- NULL
  tmp
  
}, simplify=F)

discordantPairsTab <- do.call("rbind", discordantPairsTab)


fileDiscordant <- "outputTabs/iPSC/discordant/discordantPairsTab.txt"
write.table(discordantPairsTab, file=fileDiscordant,
            quote=F,
            col.names = F,
            row.names = F,
            sep="\t")




