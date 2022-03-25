

## EXAMPLE: KMT2D in dopaminergic neurons and FPP-1 at day 11


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
library(GenomicRanges)


fileToPath <- "../singleCellProcessing/output/"
file <- paste0("allpools.scanpy.D11.wMetaClustUmapGraph.exprLogNormNotScaled.h5ad")

h5_data <- H5File$new(paste0(fileToPath,file), mode = "r")
#h5ls(file)

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

kk=h5read(paste0(fileToPath,file),"/obs/__categories")

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
metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$donor_id <- paste0(metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$donor_id, "/",
                                                                    unname(koCorr[metadata[metadata$donor_id=="HPSI0114i-kolf_2",]$pool_id]))
allbatches <- unique(batchNames$batchInfo)
metadata$donor_extended <- paste0(metadata$donor_id, "/", metadata$pool_id)

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
                               whichBatches=paste(unique(tmp2$batch), collapse=","),
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
                               whichBatches=paste(unique(tmp2$batch), collapse=","),
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
# [1] 69
summData <- subset(summDataQ20, quantile!=1)

alldonorsPerPool <- unique(metadata$donor_extended)
names(alldonorsPerPool) <- gsub("/pool.+","",alldonorsPerPool)
poolRep <- split(unname(alldonorsPerPool), names(alldonorsPerPool))
poolRep <- poolRep[elementNROWS(poolRep)>1]



geneOfSelection <- "KMT2D"
tmp_seurat <- t(as.matrix(seurat[rownames(seurat) %in% geneOfSelection,]@assays$RNA[]))

linesData <- sapply(unique(summData$donor_extended), function(y){
  #print(y)
  tmp_metadata <- subset(metadata, donor_extended==y)
  batches <- unique(unlist(sapply(strsplit(subset(summData, donor_extended==y)$whichBatches,","), function(x) x, simplify=F)))
  tmp_metadata <- tmp_metadata[tmp_metadata$batch %in% batches,]
  tmp_seurat2 <- tmp_seurat[match(tmp_metadata$cellid, rownames(tmp_seurat)),]
  lineData <- sapply(unique(subset(summData, donor_extended==y)$annot), function(z){
    
    tmp_metadata2 <- subset(tmp_metadata, annot==z)
    
    if (dim(tmp_metadata2)[1]>0){
      
      tmp_seurat3 <- tmp_seurat2[match(tmp_metadata2$cellid, names(tmp_seurat2))]
      stopifnot(dim(tmp_metadata2)[1]==length(tmp_seurat3))
      stopifnot(dim(tmp_metadata2)[1]==subset(summData, donor_extended==y & annot==z)$nCells)
      
      
      data.frame(donor_extended=y,
                 gene=geneOfSelection,
                 annot=z,
                 tp=timePoint,
                 meanExp=mean(tmp_seurat3),
                 cfrac=subset(summData, donor_extended==y & annot==z)$cfrac,
                 nCellsLineAnnot=subset(summData, donor_extended==y & annot==z)$nCells,
                 nCellsLine=unique(subset(summData, donor_extended==y)$nTotalCells))
      
    } else {
      
      
      data.frame(donor_extended=y,
                 gene=geneOfSelection,
                 annot=z,
                 tp=timePoint,
                 meanExp=NA,
                 cfrac=0,
                 nCellsLineAnnot=0,
                 nCellsLine=unique(subset(summData, donor_extended==y)$nTotalCells))
      
      
      
    }
    
    
  }, simplify=F)
  
  lineData <- do.call("rbind", lineData)
  rownames(lineData) <- NULL
  lineData
  
}, simplify=F)

linesData <- do.call("rbind", linesData)
rownames(linesData) <- NULL


clustersToSelect <- c("DA","FPP-1")
linesData <- linesData[linesData$annot %in% clustersToSelect,]

saveRDS(linesData, file="outputTabs/outAnalysis/exampleCorr_KMT2D.RDS")

