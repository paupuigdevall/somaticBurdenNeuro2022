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
library(VennDiagram)
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(RColorBrewer)
library(MASS)

source("../functionsToImport.R")

summData <- readRDS("outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)

### define successful and failed lines based on our data (to include KO, basically!)
summDataTP <- subset(summData, tp=="D52")
summDataTP <- subset(summDataTP, annot=="DA" | annot=="Sert-like")


summDataDiff <- sapply(unique(summDataTP$donor_extended), function(x) {
  
  tmp <- subset(summDataTP, x==donor_extended)
  data.frame(donor_id=unique(tmp$donor_id),
             donor_extended=unique(tmp$donor_extended),
             nBatches=paste0(tmp$nBatches, collapse=","),
             cfrac=sum(tmp$nCells)/unique(tmp$nTotalCells))
  
}, simplify=F)

summDataDiff <- do.call("rbind", summDataDiff)
rownames(summDataDiff) <- NULL

summDataDiff$outcome <- NA
mask_succ <- summDataDiff$cfrac>=0.2
summDataDiff$outcome[mask_succ] <- "Successful"
summDataDiff$outcome[!mask_succ] <- "Failed"
summDataDiff$outcome2 <- NA
## outcome2 set to outcome (for those without pool replicates)
summDataDiff[match(names(which(table(summDataDiff$donor_id)==1)), summDataDiff$donor_id),]$outcome2 <-
  summDataDiff[match(names(which(table(summDataDiff$donor_id)==1)), summDataDiff$donor_id),]$outcome
## outcome2 definition for those cell-lines with pool replicates (either both succ, both fail, or discordant)
ids_dup <- names(which(table(summDataDiff$donor_id)>1))

outDisamb <- sapply(ids_dup, function(x){
  
  tmp <- subset(summDataDiff, donor_id==x)
  
  if (length(unique(tmp$outcome))==1){
    return(unique(tmp$outcome))
  } else {
    return("Discordant")
  }
  
  
}, simplify=T)

summDataDiff[is.na(summDataDiff$outcome2),]$outcome2 <-
  unname(outDisamb[match(summDataDiff[is.na(summDataDiff$outcome2),]$donor_id, (names(outDisamb)))])
summDataDiff <- subset(summDataDiff, outcome2!="Discordant")

pathToFile <- "outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
delMutBurden <- readRDS(paste0(pathToFile, "mutBurden_del.RDS"))
synMutBurden <- readRDS(paste0(pathToFile, "mutBurden_synonymous.RDS"))

geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
mutTab_logReg <- readRDS("outputTabs/mutTab_logReg2.RDS")

pvalConverter <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.05)){
    symbolVec[vectorPval<0.05] <- "*"
  } 
  
  if (any(vectorPval<0.01)){
    symbolVec[vectorPval<0.01] <- "**"
  }
  
  if (any(vectorPval<0.001)){
    symbolVec[vectorPval<0.001] <- "***"
  }
  return(symbolVec)
}

inputGSEA_tabs <- function(ptvMutBurden, mutTab_logReg, dataSet="NeuroSeq_obs", varCategory="ptv"){
  
  match_col <- match(dataSet, colnames(mutTab_logReg))
  
  regSummary <- sapply(rownames(ptvMutBurden), function(x){
    mutTab_logReg$gene <- unname(ptvMutBurden[x,][match(mutTab_logReg$cline, names(ptvMutBurden[x,]))])
    
    df <- data.frame(gene=x,
                     pval=NA,
                     effectSize=NA,
                     nFailed=length(subset(mutTab_logReg, get(dataSet)=="Failed")$gene),
                     nSucc=length(subset(mutTab_logReg, get(dataSet)=="Successful")$gene),
                     mutFailed=sum(subset(mutTab_logReg, get(dataSet)=="Failed")$gene),
                     mutSucc=sum(subset(mutTab_logReg, get(dataSet)=="Successful")$gene),
                     genomicLength=geneUniverse[match(x, geneUniverse$symbol),]$genomic_length,
                     totBurden=sum(subset(mutTab_logReg, get(dataSet)=="Failed" | get(dataSet)=="Successful")$gene),
                     mutState=NA,
                     norm_mutFailed=NA,
                     norm_mutSucc=NA,
                     log2FC=NA)
    df$norm_mutFailed <- df$mutFailed*1000/(df$genomicLength*df$nFailed)
    df$norm_mutSucc <- df$mutSucc*1000/(df$genomicLength*df$nSucc)
    
    mask <- df$totBurden!=0
    
    if (mask){
      tmpDf <- subset(mutTab_logReg, !is.na(get(dataSet)))
      tmpDf[,dataSet] <- as.factor(tmpDf[,dataSet])
      
      second_mask <- length(unique(tmpDf$gene))==1
      
      if (!second_mask){
        df$pval <- coef(summary(glm(get(dataSet)~gene, family="binomial", tmpDf)))[,"Pr(>|z|)"][2]
        df$effectSize <- coef(summary(glm(get(dataSet)~gene, family="binomial", tmpDf)))["gene","Estimate"]
      }
      
      }
    
    df
    
  }, simplify=F)
  
  regSummary <- do.call("rbind", regSummary)
  rownames(regSummary) <- NULL
  regSummary$adjPval <- NA
  mask <- !is.na(regSummary$pval)
  regSummary$adjPval[mask] <- p.adjust(regSummary$pval[mask],"BH")
  regSummary$signif <- pvalConverter(regSummary$pval)
  regSummary$signifAdj <- "ns"
  regSummary$signifAdj[!is.na(regSummary$adjPval)] <- pvalConverter(regSummary$adjPval[!is.na(regSummary$adjPval)])

  pseudocount <- 0.001
  regSummary$log2FC <- log2((regSummary$norm_mutFailed+pseudocount)/(regSummary$norm_mutSucc+pseudocount))
  regSummary <- subset(regSummary, !is.na(log2FC))
  regSummary

  return(regSummary)
  
}

runs <- list(c("delMutBurden","NeuroSeq_obs","del"),
             c("delMutBurden","NeuroSeq_pred","del"))
names(runs) <- c("observed_del","predicted_del")

inputGSEA_res <- sapply(runs, function(y){
  inputGSEA_tabs(get(y[1]), mutTab_logReg, dataSet=y[2], varCategory=y[3])
}, simplify=F)


geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$gene, "ENTREZID", "ENSEMBL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$ENSEMBL),]
geneUniverse$entrezid <- tabCorr[match(geneUniverse$gene,tabCorr$ENSEMBL),]$ENTREZID
geneUniverse <- subset(geneUniverse, !is.na(entrezid))

inputGSEA_res <- sapply(inputGSEA_res, function(x){
  x$entrezid <- geneUniverse[match(x$gene, geneUniverse$symbol),]$entrezid
  x <- subset(x, !is.na(entrezid))
  x
}, simplify=F)


selectGenes <- function(inputGSEA_res, selection="observed_ptv"){
  
  match_idx <- match(selection, names(inputGSEA_res))
  tmp <- inputGSEA_res[[match_idx]]
  list(successfulGenes=tmp[tmp$log2FC<quantile(tmp$log2FC, c(0.1)),]$entrezid,
       failedGenes=tmp[tmp$log2FC>quantile(tmp$log2FC, c(0.9)),]$entrezid)
}


runGO <- list(c("observed_del","failedGenes", "Failed_Observed_Del"),
              c("predicted_del","failedGenes", "Failed_Predicted_Del"),
              c("observed_del","successfulGenes", "Succ_Observed_Del"),
              c("predicted_del","successfulGenes", "Succ_Predicted_Del")) 


GOresults <- sapply(runGO, function(y){
  GOenrichmentAndReport(selectGenes(inputGSEA_res, selection=y[1])[[y[2]]], unique(geneUniverse$entrezid),
                        p.value=0.05, minCount =10, minSize=20, label=y[3])
}, simplify=F)

names(GOresults) <- sapply(runGO, function(y) y[3])


suppTable3 <- do.call("rbind", GOresults)
rownames(suppTable3) <- NULL

write.table(suppTable3, file="../suppTabs/suppTable3.txt",
            sep="\t", quote=F, col.names=T, row.names = F)

