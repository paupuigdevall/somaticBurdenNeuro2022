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


#############
### fig4B ###
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

fig4b <- ggplot(data=heatmap_df, aes(x=timepoint, y=annot, fill=numDE))+
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

pdf(file="figures/mainFigs/figure4B.pdf")
plot(fig4b)
dev.off()

#############
### fig4C ###
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

suppTable5 <- rbind(report_allclusters, report_allsignif_ddd, do.call("rbind", report_testClust))
write.table(suppTable5, file="suppTabs/suppTable5.txt",
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
fig4c <- ggplot(neuroTab,
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

pdf(file="figures/mainFigs/figure4C.pdf", width = 7, height = 4)
plot(fig4c)
dev.off()

