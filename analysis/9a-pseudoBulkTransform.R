

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
gc()

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

batchNames <- read.table("../../pools_to_merge_strictTP.csv", sep=",", header = T)
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

annotLeiden <- c("DA","Astro","FPP-1","Sert-like","Epend-1","proFPP-1","FPP-2","proFPP-2","Unk-1","FPP-3","Pro.Sert-like","Unk-2")
names(annotLeiden) <- as.character(c(0:11))
seurat@meta.data$annot <- unname(annotLeiden[seurat@meta.data$leiden_id])


###
###

##discard KO

metData <- seurat@meta.data
seurat <- seurat[,!grepl("/",seurat$donorIdExpanded)]
seurat@meta.data <- metData[!grepl("/",metData$donorIdExpanded),]

summData <- readRDS("outputTabs/DEsinglecell/CTfractionPerLinePerTPCurated2.RDS")
seurat@meta.data$DNoutcome <- summData[match(seurat@meta.data$donor_id, summData$donor_id),]$outcome

## discard missing DN outcome
geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
mutTab_logReg <- readRDS("outputTabs/mutTab_logReg2.RDS")
pathToFile <- "outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))

###Comparisons of interest
seurat@meta.data$LoF_BCOR <- mutTab_logReg[match(seurat@meta.data$donor_id, mutTab_logReg$cline),]$BCORpos
stopifnot(!any(unique(seurat@meta.data[is.na(seurat@meta.data$LoF_BCOR),]$donor_id) %in% mutTab_logReg$cline))

metData <- seurat@meta.data
seurat <- seurat[,!is.na(seurat@meta.data$LoF_BCOR)]
seurat@meta.data <- metData[!is.na(metData$LoF_BCOR),]

###Comparisons of interest
bcorExprSeurat <- seurat["BCOR",]
metDataBCOR <- seurat@meta.data
bcorExprSeurat@meta.data <- metDataBCOR

##pseudoBulk per line

exprSummary <- sapply(unique(bcorExprSeurat@meta.data$donor_id), function(x){
  
  print(x)
  print(paste0("Line number ", which(unique(bcorExprSeurat@meta.data$donor_id) %in% x) ,
               " out of ", length(unique(bcorExprSeurat@meta.data$donor_id))))
  
  tmp <- bcorExprSeurat[, bcorExprSeurat@meta.data$donor_id==x]
  metData <- bcorExprSeurat@meta.data
  tmp@meta.data <- metData[bcorExprSeurat@meta.data$donor_id==x,]
  
  if (dim(tmp)[2]>10){
    data.frame(cline=x,
               numCells=dim(tmp)[2],
               meanExpPseudoBulk=mean(as.numeric(tmp@assays$RNA[])),
               DNoutcome=unique(tmp$DNoutcome),
               BCORpositive=unique(tmp$LoF_BCOR))
  } else {
    data.frame(cline=x,
               numCells=NA,
               meanExpPseudoBulk=NA,
               DNoutcome=unique(tmp$DNoutcome),
               BCORpositive=unique(tmp$LoF_BCOR))
    
  }
  
}, simplify=F)


exprSummary <- do.call("rbind", exprSummary)
rownames(exprSummary) <- NULL

exprSummary <- exprSummary[!is.na(exprSummary$meanExpPseudoBulk),]
exprSummary <- exprSummary[!is.na(exprSummary$DNoutcome),]

saveRDS(exprSummary, file=paste0("outputTabs/expr_",tp,"_pseudoBulkPerLine.RDS"))





