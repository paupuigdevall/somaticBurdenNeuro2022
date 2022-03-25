library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library(ggpubr)
library(Category)
library(GSEABase)
library(GOstats)
library(org.Hs.eg.db)
library(knitr)
library(limma)
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(ggupset)
library(wesanderson)
library(tidyverse)

source("../functionsToImport.R")

## outlier analysis
summData <- readRDS("outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)

summData <- sapply(unique(summData$tp), function(x){
  
  tmp <- subset(summData, tp==x)
  annotTmp <- sapply(unique(tmp$annot), function(y){
    
    tmp2 <- subset(tmp, annot==y)
    tmp2$zsco <- round((tmp2$cfrac-mean(tmp2$cfrac))/sd(tmp2$cfrac),2)
    tmp2$outlier <- abs(tmp2$zsco)>2
    tmp2
    
  }, simplify=F)
  
  annotTmp <- do.call("rbind", annotTmp)
  rownames(annotTmp) <- NULL
  annotTmp
  
}, simplify=F)

summData <- do.call("rbind", summData)
rownames(summData) <- NULL


summData$pool <- sapply(strsplit(summData$donor_extended, "/"), function(x) x[length(x)])
outlier_analysis <- data.frame(donor_extended=unique(summData$donor_extended),
                               outlier_tp11_z2=NA,
                               outlier_tp30_z2=NA,
                               outlier_tp52_z2=NA,
                               outlier_alltp_z2=NA)
avail <- sapply(unique(summData$donor_extended), function(x){ paste(unique(subset(summData, donor_extended==x)$tp), collapse="-")}, simplify=T)
outlier_analysis$avail <- unname(avail[match(outlier_analysis$donor_extended, names(avail))])

#all outliers (all tp)
outliers_z2 <- unique(summData[summData$outlier==TRUE,]$donor_extended)
outlier_analysis$outlier_alltp_z2 <- outlier_analysis$donor_extended %in% outliers_z2
outlier_analysis$donor_id <- sapply(strsplit(outlier_analysis$donor_extended,"/"), function(x) paste(x[-length(x)],collapse="/"))

tpOutlierSet <- function(summData, timePoint="D11", colname="outlier_tp11_z2"){
  outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE)$donor_extended)
  colIndex <- match(colname, colnames(outlier_analysis))
  outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
  return(outlier_analysis)
}

outlier_analysis <- tpOutlierSet(summData, timePoint="D11", colname="outlier_tp11_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D30", colname="outlier_tp30_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D52", colname="outlier_tp52_z2")

outlierLong <- outlier_analysis[,c("donor_extended","outlier_alltp_z2","outlier_tp11_z2","outlier_tp30_z2","outlier_tp52_z2")]
colnames(outlierLong)[2:5] <- c("agg_allTP", "Day11", "Day30", "Day52")


vec4 <- c("Outlier","Non-outlier")
names(vec4) <- c("TRUE","FALSE")
outlierLong$agg_allTP <- unname(vec4[as.character(outlierLong$agg_allTP)])
outlierLong$Day11 <- unname(vec4[as.character(outlierLong$Day11)])
outlierLong$Day30 <- unname(vec4[as.character(outlierLong$Day30)])
outlierLong$Day52 <- unname(vec4[as.character(outlierLong$Day52)])


##########

metadata <- readRDS("outputTabs/suppData1.RDS")

valid_3tp <- sapply(unique(metadata$donor_extended), function(x){
  
  tmp <- subset(metadata, donor_extended==x)
  
  if (length(unique(tmp$tp))==3){
    x
  } else {
    NA
  }
  
}, simplify=T)

valid_3tp <- unname(valid_3tp[!is.na(valid_3tp)])
metadata <- metadata[metadata$donor_extended %in% valid_3tp,]

## read lines from RDS
dirFiles <- "../demuxlet/deconvolution/"
Files <- dir(dirFiles)[grepl("_sample_list.txt",dir(dirFiles))]
datalist = lapply(paste0(dirFiles,Files), function(x) read.table(x, header=F))
names(datalist) <- gsub("_sample.+","", Files)

genes_state_KO <- c("ASXL3","SNCA","CTNNB1","TCF4","CHD2","SET","GATAD2B","TBL1XR1")
names(genes_state_KO) <- c("pool11","pool12","pool13","pool14","pool16","pool17","pool20","pool21")

tmpPerPool <- sapply(unique(metadata$pool_id), function(x){
  
  tmpPool <- subset(metadata, pool_id==x)
  tmpclines <- datalist[[x]]$V1
  
  if (any(x %in% names(genes_state_KO))){
    
    mask <- grepl("kolf_2", tmpclines)
    tmpko <- unname(genes_state_KO[match(x, names(genes_state_KO))])
    tmpclines[mask] <- paste0(tmpclines[mask],"/", tmpko)
    
  }
  
  tmp_tp <- sapply(sort(c("D0",unique(tmpPool$tp))), function(y){
    
    if (y=="D0"){
      
      tmp <- data.frame(cline=tmpclines,
                        pool=x,
                        timePoint=y,
                        cfrac=signif(1/length(tmpclines),3))
      
      
    } else {
      
      tmpPool2 <- subset(tmpPool, tp==y)
      tabNum <- table(tmpPool2$donor_id)
      
      tmp <- data.frame(cline=tmpclines,
                        pool=x,
                        timePoint=y,
                        cfrac=NA)
      tmp$cfrac <- signif(tabNum[match(tmp$cline, names(tabNum))]/sum(tabNum),3)
      
      
    }
    
    
    tmp
    
    
    
  }, simplify=F)
  
  tmp_tp <- do.call("rbind", tmp_tp)
  rownames(tmp_tp) <- NULL
  tmp_tp
  
}, simplify=F)

tmpPerPool <- do.call("rbind", tmpPerPool)
rownames(tmpPerPool) <- NULL

tmpPerPool$pool <- firstup(tmpPerPool$pool)
tmpPerPool$pool <- as.factor(tmpPerPool$pool)
tmpPerPool$pool <- factor(tmpPerPool$pool,
                          levels=levels(tmpPerPool$pool)[order(as.numeric(gsub("Pool","",levels(tmpPerPool$pool))))])


tmpPerPool$cline_expanded <- paste0(tmpPerPool$cline,"/", tmpPerPool$pool)
tmpPerPool$cfrac_log1p <- signif(log1p(tmpPerPool$cfrac),3)
ratiosTab <- sapply(unique(tmpPerPool$cline_expanded), function(z){
  
  tmp <- subset(tmpPerPool, cline_expanded==z)
  
  data.frame(cline_expanded=z,
             pool=as.character(unique(tmp$pool)),
             ratio_d11_d0=signif(subset(tmp, timePoint=="D11")$cfrac/subset(tmp, timePoint=="D0")$cfrac,3),
             ratio_d30_d0=signif(subset(tmp, timePoint=="D30")$cfrac/subset(tmp, timePoint=="D0")$cfrac,3),
             ratio_d52_d0=signif(subset(tmp, timePoint=="D52")$cfrac/subset(tmp, timePoint=="D0")$cfrac,3))
  
  
}, simplify=F)

ratiosTab <- do.call("rbind", ratiosTab)
rownames(ratiosTab) <- NULL

ratiosTab$pool <- firstup(ratiosTab$pool)
outlierLong$donor_extended <- gsub("/pool","/Pool",outlierLong$donor_extended)

## linear regression (using cell-fraction at day 52)


summData$donor_extended <- gsub("pool","Pool", summData$donor_extended)
pathToFile <- "outputTabs/iPSC/"
ptvIPSC <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
delIPSC <- readRDS(paste0(pathToFile, "mutBurden_del.RDS"))

tmpTab <- subset(summData, tp=="D52")
mergedTab <- merge(tmpTab, ratiosTab, by.x="donor_extended", by.y="cline_expanded")

allClust <- sapply(rownames(delIPSC), function(y){

  print(y)
  mergedTab$gene <- unname(delIPSC[y,][match(mergedTab$donor_id, names(delIPSC[y,]))])
  mergedTab$geneStatus <- mergedTab$gene>0
  stopifnot(all(table(mergedTab$donor_extended)==12))
  
  annotClust <- sapply(unique(mergedTab$annot), function(x){
    
    tmpTab2 <- subset(mergedTab, annot==x & !is.na(geneStatus))
    res <- glm(formula=geneStatus ~ cfrac, data=tmpTab2, family="binomial") 
    
    
    if (length(subset(tmpTab2, geneStatus==FALSE)$cfrac)>0 & length(subset(tmpTab2, geneStatus==TRUE)$cfrac)>0){
      res2 <- wilcox.test(subset(tmpTab2, geneStatus==FALSE)$cfrac,
                          subset(tmpTab2, geneStatus==TRUE)$cfrac)$p.value
      dirMut <- which.max(c(mean(subset(tmpTab2, geneStatus==FALSE)$cfrac),mean(subset(tmpTab2, geneStatus==TRUE)$cfrac)))
      vecdirMut <- c("unMut","Mut")
      names(vecdirMut) <- as.character(c("1","2"))
      dirMut <- unname(vecdirMut[as.character(dirMut)])
      
    } else{
      res2 <- NA
      dirMut <- NA
    }
    
    
    if (is.na(res$coefficients["cfrac"])){
      
      data.frame(annotId=x,
                 logRegPval=NA,
                 wilcoxPval=NA,
                 dirWilcox=NA)
      
    } else {
      
      
      
      data.frame(annotId=x,
                 logRegPval=summary(res)$coef[2,4],
                 wilcoxPval=res2,
                 dirWilcox=dirMut)
      
    }
    
  }, simplify=F)
  
  annotClust <- do.call("rbind", annotClust)
  rownames(annotClust) <- NULL
  annotClust$gene <- y
  
  return(annotClust)
  
}, simplify=F)

allClust <- do.call("rbind", allClust)
rownames(allClust) <- NULL

saveRDS(allClust, file="outputTabs/corr_cfracBurdenprolif.RDS")

