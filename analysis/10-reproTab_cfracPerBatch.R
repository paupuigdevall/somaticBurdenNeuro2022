
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

## save metadata?

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


#######################################
#### Raw cell counts per 10x inlet ####
#######################################

## expected cell-lines per pool
dirfiles <- "../demuxlet/deconvolution/"
files <- list.files(dirfiles)[grep("^pool[0-9].+.txt", list.files(dirfiles))]
linesPerPool <- sapply(files, function(x){
  tabFile <- read.table(paste0(dirfiles,x))
  colnames(tabFile) <- "cell_lines"
  tabFile
}, simplify=F)
names(linesPerPool) <- gsub("_.+","",names(linesPerPool))

## There are no independent pools (with exact same donors)
stopifnot(all(rowSums(sapply(linesPerPool, function(x){
  sapply(linesPerPool, function(y){
    all(x$cell_lines %in% y$cell_lines)
  })
}))==1))

summTab <- sapply(unique(metadata$batch), function(x) {
  #print(x)
  tmp <- subset(metadata, batch==x)
  pool <- unique(tmp$pool_id)
  
  df <- data.frame(cell_lines=linesPerPool[pool][[1]]$cell_lines,
                   batch=x,
                   pool=pool,
                   numCells=0,
                   tp=gsub("-.+","",unique(tmp$tp_treat)))
  vec <- table(tmp$donor_id)
  
  if (length(grep("kolf_2", names(vec)))){
    original <- names(vec[grep("kolf_2", names(vec))])
    mod <- gsub("/.+","", names(vec[grep("kolf_2", names(vec))]))
    names(vec)[grep("kolf_2", names(vec))] <- mod
  }
  
  df[match(names(vec), df$cell_lines),]$numCells <- unname(vec)
  df$norm_numCells <- signif(df$numCells/sum(df$numCells),3)
  if (length(grep("kolf_2", df$cell_lines))){
    df[grep("kolf_2", df$cell_lines),]$cell_lines <- original
  }
  
  df
}, simplify=F)

summTab <- do.call("rbind", summTab)
rownames(summTab) <- NULL

## everything together with mutational burden
summTab$cell_lines2 <- gsub("/.+","",summTab$cell_lines)
exomes_allinfo <- read.table("../suppTabs/suppTable1.txt",
                             header=T)
summTab <- cbind(summTab,
                 exomes_allinfo[match(summTab$cell_lines2, exomes_allinfo$donor_iPSC), c("all","synonymous","other","deleterious","missPatho","ptv")])
rownames(summTab) <- NULL

summTab <- cbind(summTab, batchNames[match(summTab$batch, batchNames$batch),c("bioRepLabel","techRepLabel","tenXRepLabel")])
rownames(summTab) <- NULL


##Re-label the replicates (per donor)
summTab <- sapply(unique(summTab$cell_lines), function(x){
  
  tmp <- subset(summTab, cell_lines==x)
  
  if (length(unique(tmp$pool))>1){
    vec <- 1:length(unique(tmp$pool))
    names(vec) <- unique(tmp$pool)
    tmp$donorRep <- unname(vec[tmp$pool])
  } else {
    
    tmp$donorRep <- 0
  }
  
  
  ##RepLabel
  if (sum(!is.na(tmp$techRepLabel))){
    newvec <- as.character(1:length(unique(gsub("-.+","",tmp$techRepLabel))))
    names(newvec) <- unique(gsub("-.+","",tmp$techRepLabel))
    tmp$techRepLabel <- paste0(unname(newvec[gsub("-.+","",tmp$techRepLabel)]),"-",gsub(".+-","",tmp$techRepLabel))
    tmp$techRepLabel[grep("NA", tmp$techRepLabel)] <- "0"
  }
  
  ##BioLabel
  if (sum(!is.na(tmp$bioRepLabel))){
    newvec <- as.character(1:length(unique(gsub("-.+","",tmp$bioRepLabel))))
    names(newvec) <- unique(gsub("-.+","",tmp$bioRepLabel))
    tmp$bioRepLabel <- paste0(unname(newvec[gsub("-.+","",tmp$bioRepLabel)]),"-",gsub(".+-","",tmp$bioRepLabel))
    tmp$bioRepLabel[grep("NA", tmp$bioRepLabel)] <- "0"
  }
  
  ##TenxLabel
  if (sum(!is.na(tmp$tenXRepLabel))){
    newvec <- as.character(1:length(unique(gsub("-.+","",tmp$tenXRepLabel))))
    names(newvec) <- unique(gsub("-.+","",tmp$tenXRepLabel))
    tmp$tenXRepLabel <- paste0(unname(newvec[gsub("-.+","",tmp$tenXRepLabel)]),"-",gsub(".+-","",tmp$tenXRepLabel))
    tmp$tenXRepLabel[grep("NA", tmp$tenXRepLabel)] <- "0"
  }
  
  
  
  tmp
  
}, simplify=F)

summTab <- do.call("rbind", summTab)
rownames(summTab) <- NULL

summTab[is.na(summTab$techRepLabel),]$techRepLabel <- 0
summTab[is.na(summTab$bioRepLabel),]$bioRepLabel <- 0
summTab[is.na(summTab$tenXRepLabel),]$tenXRepLabel <- 0

saveRDS(summTab, "outputTabs/reproTab_cfracPerBatch.RDS")

