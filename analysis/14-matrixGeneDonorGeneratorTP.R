library(hdf5r)
library(Seurat)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)
library(e1071)
library(rhdf5)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library(qvalue)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
tp <- args[1]

file <- paste0("../singleCellProcessing/output/allpools.scanpy.",tp,".wMetaClustUmapGraph.exprLogNormNotScaled.h5ad")
pathToPGS="outputTabs/DEsinglecell/"

h5_data <- H5File$new(file, mode = "r")
h5ls(file)

feature_matrix <- Matrix::sparseMatrix(
  i = h5_data[['X/indices']][],
  p = h5_data[['X/indptr']][],
  x = h5_data[['X/data']][],
  dimnames = list(
    h5_data[["/var/_index"]][],
    h5_data[["obs/_index"]][]),
  dims = as.integer(c(length(h5_data[["/var/_index"]][]),
                      length(h5_data[["obs/_index"]][]))),
  index1 = FALSE
)

seurat <- CreateSeuratObject(feature_matrix, min.cells = 0.01*dim(feature_matrix)[2])
rm(feature_matrix)

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

## Identify 10x replicates (technical, biological, 10x)
## Note that donor added in different 10x replicates are replicates on donor level
donorRep <- rowSums(table(metadata$donor_id, metadata$pool_id)>0)[rowSums(table(metadata$donor_id, metadata$pool_id)>0)>1]
#table(donorRep)

allbatches <- unique(batchNames$batchInfo)

#bioReplicate
bioRep <- allbatches[grep("BIO2", allbatches)]
vecRep <- sapply(1:length(bioRep), function(x){
  tmp1 <- bioRep[x]
  tmp2 <- gsub("BIO2","BIO1",bioRep[x])
  if (tmp2 %in% allbatches){
    c(tmp1,tmp2)
  } else {
    NA
  }
},simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
if (length(vecRep)){
  stopifnot(all(sapply(vecRep, function(x) length(unique(x)))==2))
  batchNames$bioRepLabel <- NA_character_
  for (i in 1:length(vecRep)){
    batchNames[batchNames$batchInfo %in%vecRep[[i]],]$bioRepLabel <- as.character(i)
  }
}


#techReplicate
techRep <- allbatches[grep("TEC2", allbatches)]
vecRep <- sapply(1:length(techRep), function(x){
  tmp1 <- techRep[x]
  tmp2 <- gsub("TEC2","TEC1",techRep[x])
  if (tmp2 %in% allbatches){
    c(tmp1,tmp2)
  } else {
    NA
  }
}, simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
stopifnot(all(sapply(vecRep, function(x) length(unique(x)))==2))
batchNames$techRepLabel <- NA
for (i in 1:length(vecRep)){
  batchNames[batchNames$batchInfo %in%vecRep[[i]],]$techRepLabel <- as.character(i)
}


## 10x replicate
tenXRep <- allbatches[grep("TENX2", allbatches)]
vecRep <- sapply(1:length(tenXRep), function(x){
  tmp1 <- tenXRep[x]
  tmp2 <- gsub("TENX2","TENX1",tenXRep[x])
  if (tmp2 %in% allbatches){
    c(tmp1,tmp2)
  } else {
    NA
  }
}, simplify=F)
vecRep <- vecRep[!is.na(vecRep)]
stopifnot(all(sapply(vecRep, function(x) length(unique(x)))==2))
batchNames$tenXRepLabel <- NA
for (i in 1:length(vecRep)){
  batchNames[batchNames$batchInfo %in%vecRep[[i]],]$tenXRepLabel <- as.character(i)
}

metadata <- cbind(metadata, batchNames[match(metadata$sample_id, batchNames$sample_id),][,c("techRepLabel","tenXRepLabel")])
metadata$donorRepLabel <- NA
rownames(metadata) <- NULL

## donor replicate
donorRep_uniq <- unique(names(donorRep))
for (dn in donorRep_uniq){
  reppools <- unique(subset(metadata, donor_id==dn)$batchInfo)
  vecRep <- sapply(1:length(reppools), function(x){
    tmp1 <- reppools[x]
    tmp1_sub <- strsplit(tmp1,"-")[[1]][2]
    otherbatches <- reppools[-x]
    otherbatches_sub <- sapply(strsplit(reppools[-x],"-"), function(x) x[2])
    
    if (any(otherbatches_sub %in% tmp1_sub)){
      c(tmp1,otherbatches[otherbatches_sub %in% tmp1_sub])
    } else {
      NA
    }
  }, simplify=F)
  vecRep <- vecRep[!is.na(vecRep)]
  vecRep <-sapply(vecRep, function(x) sort(x), simplify=F)
  vecRep <- vecRep[!duplicated(vecRep)]
  
  for (element in vecRep){
    
    
    
  }
  
  count=1
  
  if (length(vecRep)){
    print(dn)
    for (element in 1:length(vecRep)){
      mask_donor <- metadata$donor_id==dn
      mask_vecRep <- metadata$batchInfo %in% vecRep[[element]]
      metadata[mask_donor & mask_vecRep,]$donorRepLabel <- count
      count <- count+1
      
    }
    
  }
  
}


stopifnot(length(which(rowSums(table(metadata$donor_id, metadata$donorRepLabel)>0)!=0))==length(donorRep))

mask_na <- !is.na(metadata$bioRepLabel)
if (sum(mask_na)){
  metadata[mask_na,]$bioRepLabel <- paste0(metadata[mask_na,]$bioRepLabel,"-", sapply(strsplit(metadata[mask_na,]$batchInfo, "-"), function(x) x[3]))
}

mask_na <- !is.na(metadata$techRepLabel)
if (sum(mask_na)){
  metadata[mask_na,]$techRepLabel <- paste0(metadata[mask_na,]$techRepLabel,"-", sapply(strsplit(metadata[mask_na,]$batchInfo, "-"), function(x) x[4]))
}

mask_na <- !is.na(metadata$tenXRepLabel)
if (sum(mask_na)){
  metadata[mask_na,]$tenXRepLabel <- paste0(metadata[mask_na,]$tenXRepLabel,"-", sapply(strsplit(metadata[mask_na,]$batchInfo, "-"), function(x) x[5]))
}

mask_na <- !is.na(metadata$donorRepLabel)
if (sum(mask_na)){
  metadata[mask_na,]$donorRepLabel <- paste0(metadata[mask_na,]$donorRepLabel,"-",
                                             gsub("pool", "p", metadata[mask_na,]$pool_id),"-",
                                             gsub(".+-","",metadata[mask_na,]$donor_id))
}



seurat@meta.data <- metadata



### Generate Seurat expression object per donor (instead of per cell)
# Donor across several batches are averaged (The final object is a list of the 12 clusters per TP)
# kolf-2 line is labelled according to its CTRL or KO status. Each KO gene is treated as independent donor.
# Log Norm UMI counts


metadata$donorIdExpanded <- metadata$donor_id
genes_state_KO <- c("wt","ASXL3_hom","SNCA_hom","CTNNB1","TCF4_het","CHD2_het","SET_het","GATAD2B_het","TBL1XR1_het")
names(genes_state_KO) <- c("pool10","pool11","pool12","pool13","pool14","pool16","pool17","pool20","pool21")
poolsWithKO <- names(genes_state_KO)[names(genes_state_KO) %in% unique(metadata$pool_id)]

koGenes <- unname(genes_state_KO[match(poolsWithKO, names(genes_state_KO))])
koGenes <- gsub("_.+","", koGenes[!koGenes %in% "wt"])
stopifnot(all(koGenes %in% rownames(seurat)))



for (pool in poolsWithKO){
  geneKO <- unname(genes_state_KO[names(genes_state_KO)==pool])
  mask_KO <- grepl("kolf_2", metadata$donorIdExpanded) & (metadata$pool_id==pool)
  stopifnot(any(mask_KO))
  metadata[mask_KO,]$donorIdExpanded <- paste0(metadata[mask_KO,]$donorIdExpanded,"/",geneKO)
}


## We need to create a matrix per cell-type
listOfSeuratsObjects <- vector("list", length = length(as.character(sort(as.numeric(unique(metadata$leiden_id))))))
names(listOfSeuratsObjects) <- paste0("D30-Cluster",as.character(sort(as.numeric(unique(metadata$leiden_id)))))
count=1

unique(metadata$donorIdExpanded)

for (cluster in as.character(sort(as.numeric(unique(metadata$leiden_id))))){
  
  print(cluster)
  newLogNormCountsPerDonor <- matrix(NA, nrow=dim(seurat)[1], ncol=length(unique(metadata$donorIdExpanded)))
  rownames(newLogNormCountsPerDonor) <- rownames(seurat)
  colnames(newLogNormCountsPerDonor) <- as.character(unique(metadata$donorIdExpanded))
  newMetadataTable <- data.frame()
  
  for (don in unique(metadata$donorIdExpanded)){
    
    print(don)
    matchDon <- match(don, colnames(newLogNormCountsPerDonor))
    mask_cells <- colnames(seurat) %in% (subset(metadata, donorIdExpanded==don & leiden_id==cluster)$cellid)
    #rowMeans(seurat@assays$RNA[,mask_cells])
    #names(newLogNormCountsPerDonor[,matchDon])
    stopifnot(all(names(newLogNormCountsPerDonor[,matchDon])==names(rowMeans(as.matrix(seurat@assays$RNA[,mask_cells])))))
    
    if (sum(mask_cells)>=10){
      newLogNormCountsPerDonor[,matchDon] <- unname(rowMeans(as.matrix(seurat@assays$RNA[,mask_cells])))
      tmp_newMetadataTable <- data.frame(donorIdExpanded=don,
                                         leiden_id=cluster,
                                         nCells=dim(subset(metadata, donorIdExpanded==don & leiden_id==cluster))[1],
                                         batches=paste0(unique(subset(metadata, donorIdExpanded==don)$batch),
                                                        collapse=","))
      
      newMetadataTable <- rbind(newMetadataTable, tmp_newMetadataTable)
      
    } else {
      
      tmp_newMetadataTable <- data.frame(donorIdExpanded=don,
                                         leiden_id=cluster,
                                         nCells=dim(subset(metadata, donorIdExpanded==don & leiden_id==cluster))[1],
                                         batches=paste0(unique(subset(metadata, donorIdExpanded==don)$batch),
                                                        collapse=","))
      
      newMetadataTable <- rbind(newMetadataTable, tmp_newMetadataTable)
      
    }
  }
  
  newSeurat <- CreateSeuratObject(newLogNormCountsPerDonor)
  newSeurat@meta.data <- newMetadataTable
  listOfSeuratsObjects[[count]] <- newSeurat
  count=count+1
}
saveRDS(listOfSeuratsObjects, file=paste0(pathToPGS,tp,"NeuroExpPerClust_logNormCts_List.RDS"))




