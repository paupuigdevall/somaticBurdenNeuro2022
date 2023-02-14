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
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(MASS)

source("functionsToImport.R")

###############
### SFig 3A ###
###############

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

ratio_distr <- ggplot(ratiosTablong, aes(x=log10(ratioVal)))+
  geom_histogram(col="black", fill="white")+
  facet_wrap(~ratioComparison, ncol=1,
             labeller = labeller(ratioComparison=c("ratio_d11_d0"="D11/D0",
                                                   "ratio_d30_d0"="D30/D0",
                                                   "ratio_d52_d0"="D52/D0")))+
  #ggtitle("Cell-line proportion changes across differentiation")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        axis.text=element_text(size=20),
        strip.text.x = element_text(size = 22))+
  xlab("")+ylab("")+
  #ylab("Number of cell-lines (within a pool)")+
  #xlab("Cell-line proliferation rate (log10)")+
  geom_vline(data=subset(annotData,annot=="mean"), aes(xintercept=log10(value)), col="red", size=1)+
  geom_vline(data=subset(annotData,annot!="mean"), aes(xintercept=log10(value)), col="grey40", size=0.25, linetype="dashed")

pdf(file=paste0("figures/suppFigs/suppfig3A.pdf"))
plot(ratio_distr)
dev.off()


###############
### SFig 3B ###
###############

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


###


ratiosTablong$cline <- sapply(strsplit(ratiosTablong$cline_expanded, "/"), function(x) paste0(x[-length(x)],collapse="/"))
ratiosTablong$obsAndPredDA <- summDataDiff[match(ratiosTablong$cline, summDataDiff$donor_id),]$outcome

my_comparisons <- list(c("Failed","Successful"))


vecRato <- c("D11/D0", "D30/D0", "D52/D0")
names(vecRato) <- unique(ratiosTablong$ratioComparison)

ratiosTablong$ratioComparison <- vecRato[ratiosTablong$ratioComparison]

failedprolif <- ggplot(subset(ratiosTablong, !is.na(obsAndPredDA)),
                       aes(x=ratioComparison, y=log10(ratioVal), fill=obsAndPredDA))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha=0.2, size=3)+
  theme_bw()+
  xlab("")+ylab("")+
  scale_fill_manual(name="Differentiation outcome",values=c('#E69F00', '#56B4E9'))+
  theme(axis.text.x=element_text(angle=90, size=20, vjust=0.5, hjust=0.5),
        axis.text.y=element_text(size=20),
        legend.title=element_text(face="bold", size=22),
        legend.text=element_text(size=18),
        legend.position="top")+
  stat_compare_means(label = "p.format", size=6, method = "wilcox.test")
  #ylab("Proliferation rate per time point (in log10)")


pdf(file=paste0("figures/suppFigs/suppfig3B.pdf"))
plot(failedprolif)
dev.off()


###############
### SFig 3D ###
###############

allClust <- readRDS("analysis/outputTabs/corr_cfracBurdenprolif.RDS")

adjList <- sapply(unique(allClust$annotId), function(x){
  
  tmp <- subset(allClust, annotId==x)
  tmp <- tmp[order(tmp$logRegPval, decreasing = F),]
  tmp$adjPval <- p.adjust(tmp$wilcoxPval, "BH")
  tmp
  
}, simplify=F)


adjList_filt <- sapply(adjList, function(x){
  
  #x$adjPval <- p.adjust(x$wilcoxPval, "BH")
  subset(x, adjPval<0.05)
  
}, simplify=F)

adjList_filt <- do.call("rbind", adjList_filt)
rownames(adjList_filt) <- NULL

retain <- table(metadata$annot, metadata$tp)*100/colSums(table(metadata$annot, metadata$tp))>2
retain <- retain[,"D52"]
retain <- names(retain[retain])

adjList_filt2 <- adjList_filt[adjList_filt$annotId %in% retain,]

genesOfInterest <- sapply(adjList[retain], function(x){
  
  x[x$gene %in% unique(adjList_filt2$gene),]
  
}, simplify=F)
genesOfInterest <- do.call("rbind", genesOfInterest)
rownames(genesOfInterest) <- NULL


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

genesOfInterest$signif <- pvalConverter(genesOfInterest$adjPval)
genesOfInterest[genesOfInterest$signif=="ns",]$dirWilcox <- "none"

plotBCOR <- ggplot(data=genesOfInterest, aes(x=annotId, y=gene, fill=dirWilcox))+
  geom_tile(col="black")+
  geom_point(data=genesOfInterest, aes(x=annotId, y=gene, size=signif, shape=signif))+
  theme_bw()+
  ylab("")+xlab("")+
  scale_y_discrete(position = "right")+
  theme(legend.position="top",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.title=element_text(size=22, face="bold"),
        legend.text = element_text(size=18),
        axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5, size=22),
        plot.title=element_text(hjust=0.5, size=13, face="bold"))+
  #ggtitle("Cell-type fraction differences given deleterious mutational status (per gene)")+
  #scale_fill_manual(name="Cell-type proportion\nin BCOR+ vs BCOR- lines\n(Deleterious mutations)",
  scale_fill_manual(name="",
                    values=c("blue","orange","white"),
                    breaks=c("Mut","unMut","none"),
                    label=c("BCOR+ >> BCOR-","BCOR+ << BCOR-","No sign difference"))+
  scale_size_manual(name="Wilcox P.adj", values=c(4,4,6,8),breaks=c("ns","*","**","***"))+
  scale_shape_manual(name="Wilcox P.adj",values=c(1,16,16,16), breaks=c("ns","*","**","***"))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE),
         size=guide_legend(nrow=3,byrow=TRUE))


pdf(file=paste0("figures/suppFigs/suppfig3D.pdf"))
plot(plotBCOR)
dev.off()




