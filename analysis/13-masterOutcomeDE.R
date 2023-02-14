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
pathToDE="outputTabs/DEsinglecell/"

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
seurat@meta.data <- metadata


metadata$donorIdExpanded <- metadata$donor_id
genes_state_KO <- c("wt","ASXL3_hom","SNCA_hom","CTNNB1","TCF4_het","wt","CHD2_het","SET_het","GATAD2B_het","TBL1XR1_het")
names(genes_state_KO) <- c("pool10","pool11","pool12","pool13","pool14","pool15","pool16","pool17","pool20","pool21")
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

seurat@meta.data$donorIdExpanded <- metadata$donorIdExpanded

summData <- readRDS("outputTabs/DEsinglecell/CTfractionPerLinePerTPCurated2.RDS")
seurat@meta.data$outcomeCurated <- summData[match(seurat@meta.data$donor_id, summData$donor_id),]$outcome

metadata <- seurat@meta.data
mask_na <- !is.na(seurat@meta.data$outcomeCurated)

seurat <- seurat[,mask_na]
seurat@meta.data <- metadata[mask_na,]

annotLeiden <- c("DA","Astro","FPP1","Sertlike","Epend1","proFPP1","FPP2","proFPP2","Unk1","FPP3","ProSertlike","Unk2")
names(annotLeiden) <- as.character(c(0:11))
seurat@meta.data$annot <- unname(annotLeiden[seurat@meta.data$leiden_id])

metadata <- seurat@meta.data

## make sure there are at least 10 cells from each outcome for each cell-type
## Otherwise, the cell-type is not considered for further analysis.

comp1_fail <- as.matrix(table(subset(metadata, outcomeCurated=="Failed")$donorIdExpanded,
                      subset(metadata, outcomeCurated=="Failed")$annot))
clusters_fail <- names(which(colSums(comp1_fail)>10))
comp1_succ <- as.matrix(table(subset(metadata, outcomeCurated=="Successful")$donorIdExpanded,
                              subset(metadata, outcomeCurated=="Successful")$annot))
clusters_succ <- names(which(colSums(comp1_succ)>10))
clusters <- intersect(clusters_fail, clusters_succ)
                      
dfConverter <- function(de.response, y, tp, test){
  
  if (dim(de.response)[1]>0){
    
    de.response <- de.response[de.response$p_val_adj<0.05,]
    numDE <- dim(de.response)[1]
    
    if (numDE>0){
      
      metaData <- test@meta.data
    
      de.response$geneDE <- rownames(de.response)
      rownames(de.response) <- NULL
      de.response$ident.1 <- "Failed"
      de.response$ident.2 <- "Successful"
      de.response$annot <- y
      de.response$timepoint <- tp
      de.response <- de.response[,c(6:dim(de.response)[2],1:5)]
      de.response$numDE <- numDE
      de.response$numFailedDonors <- length(unique(subset(metaData, outcomeCurated=="Failed")$donorIdExpanded))
      de.response$numSuccessfulDonors <- length(unique(subset(metaData, outcomeCurated=="Successful")$donorIdExpanded))
      de.response$numFailedCells <- dim(subset(metaData, outcomeCurated=="Failed"))[1]
      de.response$numSuccessfulCells <- dim(subset(metaData, outcomeCurated=="Successful"))[1]
      de.response$avgFailedCellsperDonor <- round(de.response$numFailedCells/de.response$numFailedDonors,2)
      de.response$avgSuccessfulCellsperDonor <- round(de.response$numSuccessfulCells/de.response$numSuccessfulDonors,2)
      
      
      return(de.response)
      
    } else if (numDE==0){
      
      de.response <- data.frame(geneDE=NA,
                                ident.1="Failed",
                                ident.2="Successful",
                                annot=y,
                                timepoint=tp,
                                p_val=NA,
                                avg_logFC=NA,
                                pct.1=NA,
                                pct.2=NA,
                                p_val_adj=NA,
                                numDE=0,
                                de.response$numFailedDonors <- length(unique(subset(metaData, outcomeCurated=="Failed")$donorIdExpanded)),
                                de.response$numSuccessfulDonors <- length(unique(subset(metaData, outcomeCurated=="Successful")$donorIdExpanded)),
                                de.response$numFailedCells <- dim(subset(metaData, outcomeCurated=="Failed"))[1],
                                de.response$numSuccessfulCells <- dim(subset(metaData, outcomeCurated=="Successful"))[1],
                                de.response$avgFailedCellsperDonor <- round(de.response$numFailedCells/de.response$numFailedDonors,2),
                                de.response$avgSuccessfulCellsperDonor <- round(de.response$numSuccessfulCells/de.response$numSuccessfulDonors,2))
                            
      return(de.response)
      
    }
    
    
  } else {
    
    de.response <- data.frame(geneDE=NA,
                              ident.1="Failed",
                              ident.2="Successful",
                              annot=y,
                              timepoint=tp,
                              p_val=NA,
                              avg_logFC=NA,
                              pct.1=NA,
                              pct.2=NA,
                              p_val_adj=NA,
                              numDE=0,
                              de.response$numFailedDonors <- length(unique(subset(metaData, outcomeCurated=="Failed")$donorIdExpanded)),
                              de.response$numSuccessfulDonors <- length(unique(subset(metaData, outcomeCurated=="Successful")$donorIdExpanded)),
                              de.response$numFailedCells <- dim(subset(metaData, outcomeCurated=="Failed"))[1],
                              de.response$numSuccessfulCells <- dim(subset(metaData, outcomeCurated=="Successful"))[1],
                              de.response$avgFailedCellsperDonor <- round(de.response$numFailedCells/de.response$numFailedDonors,2),
                              de.response$avgSuccessfulCellsperDonor <- round(de.response$numSuccessfulCells/de.response$numSuccessfulDonors,2))
    
    return(de.response)
    
  }
}


outputCluster <- sapply(clusters, function(y){
  
  print(paste0("Processing cluster ", y))
  Idents(seurat) <- "annot"
  test = seurat[,Idents(seurat)==y]
  metData <- seurat@meta.data
  test@meta.data <- metData[metData$annot==y,]
  
  ## filter-out those cells below 10 cells
  both <- sum(test@meta.data$outcomeCurated=="Successful")>10 & sum(test@meta.data$outcomeCurated=="Failed")>10
  
  if (both){
    
    ##run DE analysis
    Idents(test) <- "outcomeCurated"
    de.response <- FindMarkers(test, ident.1 = "Failed", ident.2 = "Successful",
                               verbose = FALSE, logfc.threshold = log10(1.5)/log10(10))
    de.response <- dfConverter(de.response, y, tp, test)
    return(de.response)  
    
  } else {
    
    return(NA)
  }
  
}, simplify=F)

outputCluster <- outputCluster[!sapply(outputCluster, function(x) all(is.na(x)))]
outputCluster <- do.call("rbind", outputCluster)
rownames(outputCluster) <- NULL

saveRDS(outputCluster, 
        file=paste0(pathToDE,"resultsDEFailVsSucc2_",tp,".RDS"))  

  








