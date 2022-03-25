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

tmp <- sapply(c("D11","D30","D52"), function(x){
  pgs <- readRDS(paste0("outputTabs/DEsinglecell/",x,"NeuroExpPerClust_logNormCts_List.RDS"))
  pgs <- rownames(pgs[[1]]@assays$RNA)
  data.frame(symbol=pgs,
             tp=x)
}, simplify=F)


tmp <- do.call("rbind", tmp)
rownames(tmp) <- NULL

geneUniverse <- data.frame(symbol=unique(tmp$symbol),
                           D11=unique(tmp$symbol) %in% subset(tmp,tp=="D11")$symbol,
                           D30=unique(tmp$symbol) %in% subset(tmp,tp=="D30")$symbol,
                           D52=unique(tmp$symbol) %in% subset(tmp,tp=="D52")$symbol)

genelist <- read.table("../singleCellProcessing/gene_list.csv", sep="\t", header=F)
geneUniverse$ensembl <- genelist[match(geneUniverse$symbol, genelist$V2),]$V1

mask_na <- is.na(geneUniverse$ensembl)

vec <- sapply(geneUniverse[mask_na,]$symbol, function(x){
  
  tmp <- sapply(strsplit(x, "-"), function(x) x[1])
  genelist[grepl(tmp, genelist$V2),]$V1[2]

})

geneUniverse[match(names(vec), geneUniverse$symbol),]$ensembl <- unname(vec)
stopifnot(sum(is.na(geneUniverse$ensembl))==0)

## Download Cosmic database, version 90 (GRCh37)
cat( "Annotating COSMIC Tier 1 gene-level information \n")
path_to_cosmic_db <- "inputTabs/CosmicMutantExportCensus.tsv"
cosmic_db <- read.table(path_to_cosmic_db, sep="\t", header=TRUE)
cosmic_db <- cosmic_db[cosmic_db$Tier==1,]
cosmic_db <- as.character(unique(cosmic_db$Gene.name))

## Download DDG2P genes (24.1.2020)
cat( "Annotating DDD gene-level curated list \n")
path_to_ddd_genes <- "inputTabs/DDG2P_24_1_2020.tsv"
ddd_curated <- read.csv(path_to_ddd_genes, sep="\t")
ddd_dominantMOI <- as.character(unique(ddd_curated[ddd_curated$allelic.requirement=="monoallelic",]$gene.symbol))
ddd_curated <- as.character(unique(ddd_curated$gene.symbol))


geneUniverse$ddd <- geneUniverse$symbol %in% ddd_curated
geneUniverse$cosmic <- geneUniverse$symbol %in% cosmic_db
geneUniverse$ddd_dominantMOI <- geneUniverse$symbol %in% ddd_dominantMOI

saveRDS(geneUniverse, "outputTabs/DEsinglecell/geneUniverse_seurat.RDS")



