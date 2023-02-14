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
library(RColorBrewer)
library(MASS)

source("functionsToImport.R")
summData <- readRDS("analysis/outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)
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

### successful or failed
dopaminergic_neurons <- "analysis/inputTabs/dopNeurons/predDopaminergic.csv"
dopaminergic_neurons <- read.csv(dopaminergic_neurons)

# predicted differentiation efficiency
dopaminergic_neurons$pred_diff_efficiency <- NA
dopaminergic_neurons[dopaminergic_neurons$model_score>=0.2,]$pred_diff_efficiency <- "Successful"
dopaminergic_neurons[dopaminergic_neurons$model_score<0.2,]$pred_diff_efficiency <- "Failed"

summDataDiff$JerberObserved <-
  dopaminergic_neurons[match(summDataDiff$donor_id, dopaminergic_neurons$donor_id),]$pred_diff_efficiency
vec <- c("Successful","NA","NotAssessed","Failed")
names(vec) <- unique(summDataDiff$JerberObserved)
summDataDiff$JerberObserved <- unname(vec[summDataDiff$JerberObserved])

summDataDiff$JerberPredicted <-
  dopaminergic_neurons[match(summDataDiff$donor_id, dopaminergic_neurons$donor_id),]$pred_diff_efficiency

table(summDataDiff$outcome2, summDataDiff$JerberPredicted)
# Failed Successful
# Discordant      0          6
# Failed         51          2
# Successful      1        158

summDataDiff$duplicatedPool <- FALSE
summDataDiff[summDataDiff$donor_id %in% names(outDisamb),]$duplicatedPool <- TRUE
table(sapply(unique(summDataDiff[summDataDiff$duplicatedPool==TRUE,]$donor_id),
             function(x) unique(subset(summDataDiff, donor_id==x)$outcome2), simplify=T))
# Discordant     Failed Successful 
# 3          6         27

summData$outcome <- summDataDiff[match(summData$donor_extended, summDataDiff$donor_extended),]$outcome

summData <- summData[!summData$donor_id %in% unique(summDataDiff[summDataDiff$outcome2=="Discordant",]$donor_id),]

outcomeDetail <- sapply(split(summData$outcome, summData$donor_id), function(x) unique(x))
toAddPools <- sapply(outcomeDetail[elementNROWS(outcomeDetail)>1], function(x){
  x[!is.na(x) & length(x)>1]
})

for (i in 1:length(toAddPools)){
  summData[summData$donor_id %in% names(toAddPools)[i],]$outcome <- unname(toAddPools[i])
}

summData$measuredOut <- "observed"

## There is a list of 18 donors/pool that outcome was not measured as they were not present in D52.
## To take advantage of those samples (and avoid removing them), we use the JerberPredicted data to annotate the outcome.
## Note here that the agreement between those predicitions and our experimental data is very high for those with exp. data available.


outcomeDetail <- sapply(split(summData$outcome, summData$donor_id), function(x) unique(x))
missing <- outcomeDetail[is.na(outcomeDetail)]

mask_na <- is.na(summData$outcome)
summData[mask_na,]$outcome <- dopaminergic_neurons[match(summData[mask_na,]$donor_id, dopaminergic_neurons$donor_id),]$pred_diff_efficiency
summData[mask_na,]$measuredOut <- "predicted"

outcomeDetail <- sapply(split(summData$outcome, summData$donor_id), function(x) unique(x))
missing <- outcomeDetail[is.na(outcomeDetail)]
summData <- subset(summData, !is.na(outcome))

saveRDS(summData, file="analysis/outputTabs/DEsinglecell/CTfractionPerLinePerTPCurated2.RDS")

#############
### fig4A ###
#############

## We finally remove those annotations from which we cannot derive the outcome
pvaltp <- sapply(unique(summData$tp), function(x){
  
  tmp_summData <- subset(summData, tp==x)
  
  pvalAnnot <- sapply(unique(tmp_summData$annot), function(y){
    
    tmp_summData2 <- subset(tmp_summData, annot==y)
    model.nb = glm.nb(nCells ~  outcome + offset(log(nTotalCells)), data=tmp_summData2)
    res = coef(summary(model.nb))
    
    data.frame(annot=y,
               tp=x,
               pval.nb=res["outcomeSuccessful","Pr(>|z|)"],
               pval.wilcox=wilcox.test(subset(tmp_summData2, outcome=="Failed")$cfrac,subset(tmp_summData2, outcome=="Successful")$cfrac)$p.value)
    
    
  }, simplify=F)
  
  pvalAnnot <- do.call("rbind", pvalAnnot)
  rownames(pvalAnnot) <- NULL
  pvalAnnot
  
}, simplify=F)

pvaltp <- do.call("rbind", pvaltp)
rownames(pvaltp) <- NULL

pvaltp$pvalAdj.nb <- p.adjust(pvaltp$pval.nb, "BH")
pvaltp$signif.nb <- pvalConverter(pvaltp$pvalAdj.nb)
pvaltp$comb <- paste0(pvaltp$annot,"-", pvaltp$tp)

cfracTab <- as.data.frame(summData %>%
                            group_by(tp, annot) %>%
                            summarise_at(.vars = c("cfrac"), .funs = mean))
cfracTab$rareClust <- FALSE
cfracTab[cfracTab$cfrac<0.02,]$rareClust <- TRUE
cfracTab$comb <- paste0(cfracTab$annot,"-", cfracTab$tp)

pvaltp$rareClust <- cfracTab[match(pvaltp$comb, cfracTab$comb),]$rareClust

fig4a <- ggplot(summData, aes(x=annot, y=cfrac, fill=factor(outcome)))+
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~tp, nrow = 3)+
  geom_jitter(position=position_jitterdodge(jitter.width =0.2),size=0.25, alpha=0.25)+
  theme_bw()+
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=30, vjust=0.5, hjust=0.5, size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  xlab("")+
  ylab("Cell-type fraction")+
  #ggtitle("Cell-type composition")+
  # geom_text(data = pvaltp, aes(x = annot, y = 0.9,
  #                              label = paste0("p=",formatC(pvalAdj.nb, format = "e", digits = 1)),
  #                              face="bold"),
  #            size=3, inherit.aes=F)+
  geom_text(data = pvaltp, aes(x = annot, y = 0.9, label = signif.nb, color=rareClust),
            size=4, inherit.aes=F)+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_manual(name="Differentiation outcome", values=c("Failed"='#E69F00', "Successful"='#56B4E9'))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  scale_color_manual(name=c("Rare Cluster"),values=c("FALSE"="black","TRUE"="red"))+
  theme(legend.title=element_text(face="bold"))+
  guides(col=FALSE)


pdf(file="figures/mainFigs/figure4A.pdf")
plot(fig4a)
dev.off()


##############
### Fig 4B ###
##############

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


## do that for the three Ratios

cols <- colnames(ratiosTab)[grepl("ratio",colnames(ratiosTab))]


## outlier analysis
summData <- readRDS("tabs/summData_donorCtype_Cells.RDS")
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

pdf(file=paste0("figures/mainFigs/figure4B.pdf"))
plot(plotD52)
dev.off()


#############
### fig4C ###
#############

## DE analysis (add to scripts)
pathTODE <- "analysis/outputTabs/DEsinglecell/"
summResults <- paste0(pathTODE, list.files(path=pathTODE,
                                           pattern = "resultsDEFailVsSucc2.+.RDS")) %>%
  map(readRDS) %>% 
  bind_rows()

heatmap_df <- summResults
heatmap_df <- heatmap_df[,match(c("annot","timepoint","numDE","numFailedDonors","numSuccessfulDonors","numFailedCells","numSuccessfulCells"), colnames(heatmap_df))]
heatmap_df <- heatmap_df[!duplicated(heatmap_df),]

# creation of geneUniverse (add to scripts)
geneUniverse <- readRDS("analysis/outputTabs/DEsinglecell/geneUniverse_seurat.RDS")
pval_df <- rbind(enrichmentCalc(summResults, geneUniverse, geneSet="ddd"),
                 enrichmentCalc(summResults, geneUniverse, geneSet="cosmic"),
                 enrichmentCalc(summResults, geneUniverse, geneSet="ddd_dominantMOI"))

pval_df$pvalAdj <- p.adjust(pval_df$pval, "BH")
pval_df$signif <- pvalConverter(pval_df$pvalAdj)
pval_df$signif <- factor(pval_df$signif, levels=c("ns","*","**","***"))

genVec <- c("DDD","Cosmic-T1","DDD-Dominant")
names(genVec) <- unique(pval_df$geneSet)

pval_df$geneSet <- unname(genVec[pval_df$geneSet])

heatmap_df <- subset(heatmap_df, annot!="Unk2")
heatmap_df$annot <- as.factor(heatmap_df$annot)
heatmap_df$annot <- factor(heatmap_df$annot, levels=rev(levels(heatmap_df$annot)))
pval_df <- subset(pval_df, annot!="Unk2")
pval_df$annot <- as.factor(pval_df$annot)
pval_df$annot <- factor(pval_df$annot, levels=levels(pval_df$annot))

heatmap_df$annot <- as.factor(heatmap_df$annot)
heatmap_df$annot <- factor(heatmap_df$annot, levels=rev(levels(pval_df$annot)))

size_values=c("ns"=0,
              "*"=4,
              "**"=6,
              "***"=8)

myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
my_breaks <- c(0,100,200,300,400,500)

fig4c <- ggplot(data=heatmap_df, aes(x=timepoint, y=annot, fill=numDE))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  scale_fill_gradientn(name="DE genes",
                       colours=myPalette(100),
                       limits=c(0,550),
                       labels=my_breaks,
                       breaks=my_breaks)+
  xlab("Time-point")+
  ylab("Cell-types")+
  #ggtitle("Number of DE genes / Gene set enrichment")+
  geom_point(data=pval_df, aes(x=xpos, y=annot, size=signif, col=geneSet), inherit.aes = F)+
  scale_color_manual(name="Gene-Set", values=c("DDD"="cyan1","Cosmic-T1"="plum1","DDD-Dominant"="cornflowerblue"))+
  theme(legend.title=element_text(face="bold"))+
  scale_size_manual(name="Adj.Pval", values=size_values)+
  guides(shape=F)

pdf(file="figures/mainFigs/figure4C.pdf")
plot(fig4c)
dev.off()

#############
### fig4D ###
#############

stopifnot(all(!is.na(match(summResults$geneDE, geneUniverse$symbol))))

tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$ensembl, "ENTREZID", "ENSEMBL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$ENSEMBL),]

geneUniverse$entrezid <- tabCorr[match(geneUniverse$ensembl,tabCorr$ENSEMBL),]$ENTREZID
geneUniverse <- subset(geneUniverse, !is.na(entrezid))

summResults$entrezid <- geneUniverse[match(summResults$geneDE, geneUniverse$symbol),]$entrezid
summResults$comb <- paste0(summResults$timepoint,"-", summResults$annot)
summResults <- subset(summResults, !is.na(entrezid))
pval_df$comb <- paste0(pval_df$timepoint,"-", pval_df$annot)

## all DE considering all clusters together
allclusters <- unique(summResults$entrezid)
geneUniverse_allclusters <- unique(geneUniverse$entrezid)

## all DE considering only signif clusters (either ddd or cosmic)
mask_ddd <- pval_df$geneSet=="DDD"
mask_cosmic <- pval_df$geneSet=="Cosmic-T1"

allsignif_ddd <- pval_df[mask_ddd & pval_df$signif!="ns",]$comb
#allsignif_cosmic <- pval_df[mask_cosmic & pval_df$signif!="ns",]$comb

allsignif_ddd <- unique(summResults[summResults$comb %in% allsignif_ddd,]$entrezid)
#allsignif_cosmic <- unique(summResults[summResults$comb %in% allsignif_cosmic,]$entrezid)
geneUniverse_allsignif <- geneUniverse_allclusters


## list of DE for each signif cluster (ddd)
list_clustr_ddd <- sapply(pval_df[mask_ddd & pval_df$signif!="ns",]$comb, function(x){
  subset(summResults, comb==x)$entrezid
}, simplify=F)

list_geneUniverse_ddd <- gsub("-.+","",pval_df[mask_ddd & pval_df$signif!="ns",]$comb)
list_geneUniverse_ddd <- sapply(list_geneUniverse_ddd, function(x){
  colMatch <- match(x, colnames(geneUniverse))
  unique(geneUniverse[geneUniverse[,colMatch]==TRUE,]$entrezid)
}, simplify=F)

## list of DE for each signif cluster (ddd)
list_clustr_cosmic <- sapply(pval_df[mask_cosmic & pval_df$signif!="ns",]$comb, function(x){
  subset(summResults, comb==x)$entrezid
}, simplify=F)

list_geneUniverse_cosmic <- gsub("-.+","",pval_df[mask_cosmic & pval_df$signif!="ns",]$comb)
list_geneUniverse_cosmic <- sapply(list_geneUniverse_cosmic, function(x){
  colMatch <- match(x, colnames(geneUniverse))
  unique(geneUniverse[geneUniverse[,colMatch]==TRUE,]$entrezid)
}, simplify=F)


## GOenrichment tests
report_allclusters <- GOenrichmentAndReport(allclusters, geneUniverse_allclusters, minSize=30, maxSize=200, minCount=20, p.value=0.001, label="allDE")
##allsignif-ddd
report_allsignif_ddd <- GOenrichmentAndReport(allsignif_ddd, geneUniverse_allsignif, minSize=30, maxSize=200, minCount=20, p.value=0.001, label="allsignif")
report_allclusters$GeneSyms <- NULL
report_allsignif_ddd$GeneSyms <- NULL

##ddd
report_testClust <- sapply(1:length(list_clustr_ddd), function(x){
  print(x)
  labelClust= names(list_clustr_ddd)[x]
  report_testClust <- GOenrichmentAndReport(list_clustr_ddd[[x]], list_geneUniverse_ddd[[x]], minSize=10, maxSize=200, minCount=7, p.value=0.001, label=labelClust)
  report_testClust$GeneSyms <- NULL
  report_testClust
}, simplify=F)
#saveRDS(report_testClust, file="analysis/outputTabs/DEsinglecell/report_testClust_ddd2.RDS")
#report_testClust <- readRDS("analysis/outputTabs/DEsinglecell/report_testClust_ddd2.RDS")

suppTable7 <- rbind(report_allclusters, report_allsignif_ddd, do.call("rbind", report_testClust))
write.table(suppTable7, file="suppTabs/suppTable7.txt",
            quote=F, sep="\t", col.names=T, row.names=F)


addPosition <- function(report_allclusters){
  report_allclusters$position <- 1:dim(report_allclusters)[1]
  return(report_allclusters)
}

report_allclusters <- addPosition(report_allclusters)
report_allsignif_ddd <- addPosition(report_allsignif_ddd)
report_testClust <- sapply(report_testClust, function(x){
  addPosition(x)
}, simplify=F)
report_testClust <- do.call("rbind", report_testClust)
reportInfo <- rbind(report_allclusters, report_allsignif_ddd, report_testClust)

vecMatch <- c("axon","neuron","glial","brain",
              "hindbrain","forebrain","midbrain","synapse","chromatin",
              "cerebellum", "neural","cortex","neurogenesis",
              "axonogenesis","nervous","hippocampus","neurotransmitter",
              "dopaminergic", "axenome", "action potential","synaptic")

vecLogic <- sapply(vecMatch, function(x){
  grepl(x,reportInfo$Term)
}, simplify=T)

reportInfo$neuroRelated <- rowSums(vecLogic)>0

reportInfo$top1 <- FALSE
reportInfo[match(unique(reportInfo$label), reportInfo$label),]$top1 <- TRUE

reportInfo$top2 <- FALSE
reportInfo[sort(c(match(unique(reportInfo$label), reportInfo$label), (match(unique(reportInfo$label), reportInfo$label))+1)),]$top2 <- TRUE

reportInfo$mostShared <- !is.na(match(reportInfo$Term,names(sort(table(reportInfo$Term), decreasing=T)[1:20])))

### Neuro-related ###

neuroTab <- subset(reportInfo, neuroRelated==TRUE)
neuroTab <- rbind(neuroTab, fillNA(neuroTab, list_clustr_ddd, highlight="neuroRelated"))
maxOddsNeuro <- ceiling(max(subset(reportInfo, neuroRelated==TRUE)$OddsRatio))
fig4d <- ggplot(neuroTab,
                aes(y=Term, x=label, fill=OddsRatio))+
  geom_tile(colour = "black")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold", size=13))+
  xlab("")+
  ylab("")+
  #ggtitle("GO:BP enrichment in critical cell-types")+
  scale_fill_gradientn(name="OddsRatio",
                       colours=myPalette(100),
                       limits=c(0,maxOddsNeuro),
                       labels=seq(0,maxOddsNeuro,2),
                       breaks=seq(0,maxOddsNeuro,2),
                       na.value = 'grey90')+
  geom_text(aes(label=position))

pdf(file="figures/mainFigs/figure4D.pdf", width = 7, height = 4)
plot(fig4d)
dev.off()

