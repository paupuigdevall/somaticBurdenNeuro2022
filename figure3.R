library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(dplyr)
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
#library(VennDiagram)
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(MASS)

source("functionsToImport.R")


##############
### Fig 3A ###
##############

summTab <- readRDS("analysis/outputTabs/reproTab_cfracPerBatch.RDS")

reproTab <- rbind(functionReproTab(summTab, repro="bioRep"),
                  functionReproTab(summTab, repro="techRep"),
                  functionReproTab(summTab, repro="tenXRep"),
                  functionReproTab(summTab, repro="donorRep"))

newlab <- c("bioRep"="Biological replicates",
            "donorRep"="Pool replicates",
            "techRep"="Technical replicates",
            "tenXRep"="10x replicate")

reproTab$type <- factor(reproTab$type, levels=c("donorRep","bioRep","techRep","tenXRep"))
reproPlot <- ggplot(reproTab, aes(x=rep.x, y=rep.y))+
  geom_point(alpha=0.5, size=2)+
  facet_wrap(~type, labeller = labeller(type=newlab))+
  xlab("Mean cell line proportion for Replicate 1")+
  ylab("Mean cell line proportion for Replicate 2")+
  #ggtitle("Reproducibility on line proportion")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        strip.text = element_text(size = 14))+
  geom_smooth(method="lm", col="grey40", alpha=0.25, size=0.5)


corr <- sapply(unique(reproTab$type), function(x){
  
  tmp <- subset(reproTab, type==x)
  fit11 <- lm(rep.y ~ rep.x, data = tmp)
  
  data.frame(type=x,
             adj.r.squared=signif(summary(fit11)$adj.r.squared, 3),
             p=signif(summary(fit11)$coef[2,4], 3))
  
}, simplify=F)

corr <- do.call("rbind", corr)
rownames(corr) <- NULL
rownames(corr) <- corr$type

reproPlot <- reproPlot +
  geom_text(data=corr,
            aes(label=paste0("italic(R) ^ 2 == ", adj.r.squared)), x=0.1, y=0.5, col="black",size=4, parse= T)+
  geom_text(data=corr,
            aes(label=paste0("italic(p) == ", p)), x=0.1, y=0.45, col="black",size=4, parse= T)

pdf(file=paste0("figures/mainFigs/figure3A.pdf"))
plot(reproPlot)
dev.off()


##############
### Fig 3B ###
##############


##raw data
metadata <- readRDS("analysis/outputTabs/suppData1.RDS")
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

dirFiles <- "demuxlet/deconvolution/"
Files <- dir(dirFiles)[grepl("_sample_list.txt",dir(dirFiles))]
datalist = lapply(paste0(dirFiles,Files), function(x)read.table(x, header=F))
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



poolComposition <- ggplot(tmpPerPool, aes(x=timePoint, y=cfrac, col=cline, group=cline))+
  geom_point()+
  geom_line()+
  facet_wrap(~pool)+
  theme_bw()+
  theme(legend.position="none",
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        strip.text = element_text(size = 14))+
  xlab("Differentiation timepoint")+
  ylab("Cell-line proportion")

pdf(file=paste0("figures/mainFigs/figure3B.pdf"))
plot(poolComposition)
dev.off()


##############
### Fig 3C ###
##############

pathToFile <- "analysis/outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
bcorPositive <- names(which(ptvMutBurden["BCOR",]>0))

##raw data
metadata <- readRDS("analysis/outputTabs/suppData1.RDS")
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

dirFiles <- "demuxlet/deconvolution/"
Files <- dir(dirFiles)[grepl("_sample_list.txt",dir(dirFiles))]
datalist = lapply(paste0(dirFiles,Files), function(x)read.table(x, header=F))
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

cols <- colnames(ratiosTab)[grepl("ratio",colnames(ratiosTab))]

annotData <- sapply(cols, function(y){
  
  vecValues <- ratiosTab[,y]
  a <- mean(vecValues)
  s <- sd(vecValues)
  n <- length(vecValues)
  error <- qt(0.975,df=n-1)*s/sqrt(n)
  left <- a-error
  right <- a+error
  
  annotData <- data.frame(ratioComparison=rep(y,3),
                          value=c(signif(a,3),
                                  signif(left,3),
                                  signif(right,3)),
                          annot=c("mean","95CI_low","95CI_high"))
  
}, simplify=F)

annotData <- do.call("rbind", annotData)
rownames(annotData) <- NULL

ratiosTablong <- as.data.frame(ratiosTab %>% pivot_longer(-c(1,2), names_to="ratioComparison", values_to="ratioVal"))

summData <- readRDS("analysis/outputTabs/summData_donorCtype_Cells.RDS")
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

ratiosTablong$cline <- sapply(strsplit(ratiosTablong$cline_expanded, "/"), function(x) paste0(x[-length(x)],collapse="/"))
ratiosTablong$obsAndPredDA <- summDataDiff[match(ratiosTablong$cline, summDataDiff$donor_id),]$outcome

my_comparisons <- list(c("Failed","Successful"))
vecRato <- c("D11/D0", "D30/D0", "D52/D0")
names(vecRato) <- unique(ratiosTablong$ratioComparison)

ratiosTablong$ratioComparison <- vecRato[ratiosTablong$ratioComparison]
ratiosTablong$newLabel <- ratiosTablong$obsAndPredDA
ratiosTablong$newLabel[ratiosTablong$cline %in% bcorPositive] <- "Failed/BCOR+"
ratiosTablong[ratiosTablong$newLabel=="Failed",]$newLabel <- "Failed/BCOR-"

my_comparisonsTest <- list(c("Failed/BCOR-","Successful"), c("Failed/BCOR+","Successful"),c("Failed/BCOR+","Failed/BCOR-"))

ratiosTablong$newLabel <- factor(ratiosTablong$newLabel, levels=c("Failed/BCOR+","Failed/BCOR-","Successful"))
bcorGrowthRate <- ggplot(subset(ratiosTablong, !is.na(newLabel)),
                         aes(x=newLabel, y=log10(ratioVal)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_jitter(width = 0.2), alpha=0.2, size=2)+
  theme_bw()+
  xlab("")+
  facet_wrap(~ratioComparison)+
  #scale_fill_manual(name="Differentiation outcome",values=c('#E69F00', '#56B4E9','grey60'))+
  theme(axis.text.x=element_text(angle=60, size=12, vjust=0.5),
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        plot.title=element_text(hjust=0.5, face="bold", size=13),
        legend.title=element_text(face="bold"),
        strip.text = element_text(size = 14),
        legend.position="top")+
  #stat_compare_means(label = "p.signif", size=4, method = "wilcox.test", comparisons=my_comparisonsTest)+
  stat_compare_means(
    comparisons = my_comparisonsTest,
    label = "wilcox.test", step.increase = 0.07, size=4
  )+
  ylab("Cell-line proliferation rate (in log10)")

pdf(file=paste0("figures/mainFigs/figure3C.pdf"))
plot(bcorGrowthRate)
dev.off()


##############
### Fig 3D ###
##############

## outlier analysis
summData <- readRDS("analysis/outputTabs/summData_donorCtype_Cells.RDS")
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
summData$donor_extended <- gsub("pool","Pool", summData$donor_extended)

#####

allClust <- readRDS("analysis/outputTabs/corr_cfracBurdenprolif.RDS")
adjList <- sapply(unique(allClust$annotId), function(x){
  
  tmp <- subset(allClust, annotId==x)
  tmp <- tmp[order(tmp$logRegPval, decreasing = F),]
  tmp$adjPval <- p.adjust(tmp$wilcoxPval, "BH")
  tmp
  
}, simplify=F)


adjList_filt <- sapply(adjList, function(x){
    subset(x, adjPval<0.05)
  
}, simplify=F)

adjList_filt <- do.call("rbind", adjList_filt)
rownames(adjList_filt) <- NULL

retain <- table(metadata$annot, metadata$tp)*100/colSums(table(metadata$annot, metadata$tp))>2
retain <- retain[,"D52"]
retain <- names(retain[retain])

adjList_filt2 <- adjList_filt[adjList_filt$annotId %in% retain,]

proliferationAndCellFraction <- function(ratio="ratio_d11_d0", timePoint="D11"){
  
  tmpCfrac <- subset(summData, tp==timePoint)
  tmpCfrac$prolifRate <- ratiosTab[match(tmpCfrac$donor_extended, ratiosTab$cline_expanded),ratio]
  tmpCfrac$ratioType <- ratio
  
  majorTypes <-names(which(sapply(unique(tmpCfrac$annot), function(x){
    mean(subset(summData, tp==timePoint & annot==x)$cfrac)>0.02
  })))
  
  annotCorr <- sapply(majorTypes, function(x){
    
    tmpCfrac_annot <- subset(tmpCfrac, annot==x)
    res <- lm(cfrac~prolifRate, data=tmpCfrac_annot)
    data.frame(annot=x,
               corrPearson=unname(cor.test(tmpCfrac_annot$prolifRate, tmpCfrac_annot$cfrac, family="pearson")$estimate),
               pval_lm=summary(res)$coef["prolifRate","Pr(>|t|)"],
               effect_size_lm=summary(res)$coef["prolifRate","Estimate"],
               tp=timePoint,
               prolifRate=ratio)
    
  }, simplify=F)
  
  annotCorr <- do.call("rbind", annotCorr)
  rownames(annotCorr) <- NULL
  annotCorr$pAdj <- p.adjust(annotCorr$pval_lm, "BH")
  annotCorr$signif <- pvalConverter(annotCorr$pAdj)
  
  return(annotCorr)
  
}

prolifTable <- proliferationAndCellFraction(ratio="ratio_d52_d0", timePoint="D52")

signif.floor <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- floor(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

plotD52 <- ggplot(data=prolifTable, aes(x=annot, y=tp, fill=corrPearson))+
  geom_tile(col="black")+
  geom_point(data=prolifTable, aes(x=annot, y=tp, size=signif, shape=signif))+
  theme_bw()+
  ylab("")+xlab("")+
  scale_y_discrete(position = "right")+
  theme(legend.position="top",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5, size=14),
        plot.title=element_text(hjust=0.5, size=13, face="bold"),
        legend.text=element_text(size=12))+
  #ggtitle("Cell-line proliferation ~ Cell-type fraction (Day 52)")+
  scale_size_manual(name="LinearReg P.adj", values=c(4,4,6,8),breaks=c("ns","*","**","***"))+
  scale_shape_manual(name="LinearReg P.adj",values=c(1,16,16,16), breaks=c("ns","*","**","***"))+
  scale_fill_gradientn(name="Pearson correlation",
                       colours = c("orange",
                                   "white", "blue"),
                       limits=c(signif.floor(min(prolifTable$corrPearson),1),
                                signif.ceiling(max(prolifTable$corrPearson),1)),
                       labels=seq(signif.floor(min(prolifTable$corrPearson),1),signif.ceiling(max(prolifTable$corrPearson),1),0.2),
                       breaks=seq(signif.floor(min(prolifTable$corrPearson),1),signif.ceiling(max(prolifTable$corrPearson),1),0.2),
                       guide = guide_colourbar(barwidth = 10, nbin = 5))+
  guides(size=guide_legend(nrow=3,byrow=TRUE))

pdf(file=paste0("figures/mainFigs/figure3D.pdf"))
plot(plotD52)
dev.off()





