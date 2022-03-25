library(hdf5r)
library(Seurat)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(e1071)
library(rhdf5)
library(scales)

pathToOutput <- "../singleCellProcessing/output/"
file <- paste0(pathToOutput,"allpools.scanpy.harmonyPCA.clust.wmeta.notNorm.h5ad")
h5_data <- H5File$new(file, mode = "r")
h5ls(file)

metadata <- data.frame(cellid=h5_data[["obs/_index"]][],
                       batch=h5_data[["obs/batch"]][],
                       donor_id=h5_data[["obs/donor_id"]][],
                       leiden_id=h5_data[["obs/leiden_id"]][],
                       pool_id=h5_data[["obs/pool_id"]][],
                       sample_id=h5_data[["obs/sample_id"]][],
                       tp_treat=h5_data[["obs/tp_treat"]][])
kk=h5read(file,"/obs/__categories")

print(object.size(metadata), units = "auto")

relabelColumns <- function(metadata, kk){
  corr <-sapply(kk, function(x) {
    tt <- seq(0,length(x),1)
    tt <- tt[-length(tt)]
    tt2 <- x
    names(tt2) <- tt
    tt2
  })
  colsToMatch <- colnames(metadata)[colnames(metadata) %in% names(corr)]
  for (col in colsToMatch){
    mtch_corr <- match(col, names(corr))
    tmp <- corr[mtch_corr][[1]]
    metadata[,match(col, colnames(metadata))] <- tmp[as.character(metadata[,match(col, colnames(metadata))])]
  }
  return(metadata)
}

metadata <- relabelColumns(metadata, kk)
batchNames <- read.table("../singleCellProcessing/pools_to_merge_strictTP.csv", sep=",", header = T)
batchNames <- batchNames[batchNames$sample_id %in% metadata$sample_id,]
batchNames$batchInfo <- paste0(batchNames$pool_id,"-",batchNames$time_point,
                               "-BIO", batchNames$bioRep, "-TEC", batchNames$techRep,"-TENX", batchNames$tenXRep)
metadata$batchInfo <- batchNames$batchInfo[match(metadata$sample_id, batchNames$sample_id)]
batchNames$batch <- metadata[match(batchNames$sample_id, metadata$sample_id),]$batch
batchNames <- batchNames[!is.na(batchNames$batch),]

koCorr <- c("WT","ASXL3","SNCA","CTNNB1","TCF4","WT","CHD2","SET","GATAD2B","TBL1XR1")
names(koCorr) <- paste0("pool",c(10:17,20,21))


### treat kolf2 KO as different lines

metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$donor_id <- paste0(metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$donor_id, "/",
                                                                    unname(koCorr[metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$pool_id]))


allbatches <- unique(batchNames$batchInfo)

### Biological replicates ###
bioRep <- allbatches[grep("BIO2", allbatches)]
vecRep <- sapply(1:length(bioRep), function(x){
  tmp1 <- gsub("-TEC.+","",bioRep[x])
  tmp2 <- gsub("BIO2","BIO1",tmp1)
  if (length(grep(tmp2, allbatches))){
    c(allbatches[grep(tmp1, allbatches)], allbatches[grep(tmp2, allbatches)])
  } else {
    NA
  }
},simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
vecRep <- vecRep[!duplicated(vecRep)]
stopifnot(all(sapply(vecRep, function(x) length(unique(x))>1)))
batchNames$bioRepLabel <- NA_character_
vecUpd <- as.character(1:length(unique(gsub("-BIO.+","",unlist(vecRep)))))
names(vecUpd) <- unique(gsub("-BIO.+","",unlist(vecRep)))


for (i in 1:length(vecRep)){
  
  fullName <- batchNames[batchNames$batchInfo %in% vecRep[[i]],]$batchInfo
  rootName <- unname(vecUpd[gsub("-BIO.+","", fullName)])
  shortName <- gsub("BIO","",sapply(strsplit(batchNames[batchNames$batchInfo %in%vecRep[[i]],]$batchInfo, "-"), function(x) x[3]))
  
  batchNames[batchNames$batchInfo %in%vecRep[[i]],]$bioRepLabel <-
    paste0(rootName,"-",shortName)
}


### Technical replicates ###
techRep <- allbatches[grep("TEC2", allbatches)]
vecRep <- sapply(1:length(techRep), function(x){
  tmp1 <- gsub("-TENX.+","",techRep[x])
  tmp2 <- gsub("TEC2","TEC1",tmp1)
  if (length(grep(tmp2, allbatches))){
    c(allbatches[grep(tmp1, allbatches)], allbatches[grep(tmp2, allbatches)])
  } else {
    NA
  }
},simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
vecRep <- vecRep[!duplicated(vecRep)]

stopifnot(all(sapply(vecRep, function(x) length(unique(x))>1)))
batchNames$techRepLabel <- NA_character_
vecUpd <- as.character(1:length(unique(gsub("-TEC.+","",unlist(vecRep)))))
names(vecUpd) <- unique(gsub("-TEC.+","",unlist(vecRep)))


for (i in 1:length(vecRep)){
  fullName <- batchNames[batchNames$batchInfo %in%vecRep[[i]],]$batchInfo
  rootName <- unname(vecUpd[gsub("-TEC.+","", fullName)])
  shortName <- gsub("TEC","",sapply(strsplit(batchNames[batchNames$batchInfo %in%vecRep[[i]],]$batchInfo, "-"), function(x) x[4]))
  batchNames[batchNames$batchInfo %in%vecRep[[i]],]$techRepLabel <-
    paste0(rootName,"-",shortName)
}


### TenX replicates ###
tenXRep <- allbatches[grep("TENX2", allbatches)]
vecRep <- sapply(1:length(tenXRep), function(x){
  tmp1 <- tenXRep[x]
  tmp2 <- gsub("TENX2","TENX1",tmp1)
  if (length(grep(tmp2, allbatches))){
    c(allbatches[grep(tmp1, allbatches)], allbatches[grep(tmp2, allbatches)])
  } else {
    NA
  }
},simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
vecRep <- vecRep[!duplicated(vecRep)]

stopifnot(all(sapply(vecRep, function(x) length(unique(x))>1)))
batchNames$tenXRepLabel <- NA_character_
vecUpd <- as.character(1:length(unique(gsub("-TENX.+","",unlist(vecRep)))))
names(vecUpd) <- unique(gsub("-TENX.+","",unlist(vecRep)))

for (i in 1:length(vecRep)){
  fullName <- batchNames[batchNames$batchInfo %in%vecRep[[i]],]$batchInfo
  rootName <- unname(vecUpd[gsub("-TENX.+","", fullName)])
  shortName <- gsub("TENX","",sapply(strsplit(batchNames[batchNames$batchInfo %in%vecRep[[i]],]$batchInfo, "-"), function(x) x[5]))
  batchNames[batchNames$batchInfo %in%vecRep[[i]],]$tenXRepLabel <-
    paste0(rootName,"-",shortName)
}

metadata <- cbind(metadata, batchNames[match(metadata$sample_id, batchNames$sample_id),][,c("bioRepLabel","techRepLabel","tenXRepLabel")])
metadata$donorRepLabel <- NA
rownames(metadata) <- NULL


metadata$donor_extended <- paste0(metadata$donor_id, "/", metadata$pool_id)

# > length(unique(metadata$donor_extended))
# [1] 281
# > length(unique(metadata$donor_id))
# [1] 236

annotLeiden <- c("DA","Astro","FPP-1","Sert-like","Epend-1","proFPP-1","FPP-2","proFPP-2","Unk-1","FPP-3","Pro.Sert-like","Unk-2")
names(annotLeiden) <- as.character(c(0:11))
metadata$annot <- unname(annotLeiden[metadata$leiden_id])
metadata$tp <- gsub("-no","", metadata$tp_treat)

summData <- sapply(unique(metadata$tp)[order(unique(metadata$tp))], function(x){
  
  tmp <- subset(metadata, tp==x)
  
  summData <- sapply(unique(tmp$donor_extended), function(y){
    
    tmp2 <- subset(tmp, donor_extended==y)
    
    
    summData <- sapply(unname(annotLeiden), function(z){
      
      tmp3 <- subset(tmp2, annot==z)
      
      if (dim(tmp3)[1]){
        
        summData <- data.frame(donor_id=unique(tmp3$donor_id),
                               donor_extended=unique(tmp3$donor_extended),
                               nBatches=length(unique(tmp2$batchInfo)),
                               tp=unique(tmp3$tp),
                               annot=unique(tmp3$annot),
                               nCells=dim(tmp3)[1],
                               nTotalCells=dim(tmp2)[1])
        summData$cfrac <- summData$nCells/summData$nTotalCells
        summData
        
      } else {
        
        summData <- data.frame(donor_id=unique(tmp2$donor_id),
                               donor_extended=unique(tmp2$donor_extended),
                               nBatches=length(unique(tmp2$batchInfo)),
                               tp=unique(tmp2$tp),
                               annot=z,
                               nCells=dim(tmp3)[1],
                               nTotalCells=dim(tmp2)[1])
        summData$cfrac <- summData$nCells/summData$nTotalCells
        summData
        
      }
      
      
      
    }, simplify=F)
    
    summData <- do.call("rbind", summData)
    rownames(summData) <- NULL
    summData
    
  }, simplify=F)
  
  summData <- do.call("rbind", summData)
  rownames(summData) <- NULL
  summData
  
}, simplify=F)

summData <- do.call("rbind", summData)
rownames(summData) <- NULL


labelDivisionsSumm <- function(summData, numDiv=10){
  
  summDataQ <- sapply(unique(summData$tp)[order(unique(summData$tp))], function(x){
    
    tmpsummDataQ <- subset(summData, tp==x)
    tmpsummDataQ <- tmpsummDataQ[order(tmpsummDataQ$nTotalCells),]
    cuts <- seq(0,1, by = 1/numDiv)
    cuts <- quantile(tmpsummDataQ$nTotalCells, cuts)
    print(unname(cuts[2]))
    tmpsummDataQ$quantile <- cut(tmpsummDataQ$nTotalCells, cuts, include.lowest = T)
    vec <- 1:length(levels(tmpsummDataQ$quantile))
    names(vec) <- levels(tmpsummDataQ$quantile)
    tmpsummDataQ$quantile <- unname(vec[tmpsummDataQ$quantile])
    tmpsummDataQ
    
  }, simplify=F)
  
  summDataQ <- do.call("rbind", summDataQ)
  rownames(summDataQ) <- NULL
  
  
  return(summDataQ)
  
}

summDataQ20 <- labelDivisionsSumm(summData, numDiv=20)

saveRDS(metadata, file="outputTabs/suppData1.RDS")
saveRDS(summDataQ20, file="outputTabs/summData_donorCtype_Cells.RDS")



