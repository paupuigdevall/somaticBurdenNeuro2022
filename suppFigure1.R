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
library(VennDiagram)
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(MASS)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)



###############
### SFig 1A ###
###############

source("functionsToImport.R")
diffmetrics <- read.table("suppTabs/suppTable1.txt",
                          header=T)

colnames_of_interest <- c("donor_iPSC","all","synonymous","other","deleterious","missPatho","ptv")
ctypes <- c("sensory_neurons","macrophages","dopaminergic_neurons_exp","dopaminergic_neurons_pred")

outcomes <- sapply(ctypes, function(x){
  
  tmp <- diffmetrics[,match(c(x,colnames_of_interest), colnames(diffmetrics))]
  matchCol <- match(x, colnames(tmp))
  colnames(tmp)[matchCol] <- "outcome"
  order <- c("donor_iPSC","outcome",colnames_of_interest[2:length(colnames_of_interest)])
  tmp <- tmp[,match(colnames(tmp), order)]
  tmp <- as.data.frame(tmp %>% pivot_longer(-c("donor_iPSC","outcome"), names_to="varCategory", values_to="mutBurden"))
  tmp <- tmp[!is.na(tmp$mutBurden) & !is.na(tmp$outcome),]
  tmp$differentiation <- x
  tmp
  
}, simplify=F)

outcomes <- do.call("rbind", outcomes)
rownames(outcomes) <- NULL

newvec <- c("SN","M","DAexp","DApred")
names(newvec) <- unique(outcomes$differentiation)
outcomes$differentiation2 <- unname(newvec[outcomes$differentiation])
outcomes$differentiation2 <- paste0(outcomes$differentiation2, "_", sapply(strsplit(outcomes$outcome,""), function(x) x[1]))

sn <- table(subset(outcomes, differentiation=="sensory_neurons" & varCategory=="all")$outcome)
ma <- table(subset(outcomes, differentiation=="macrophages" & varCategory=="all")$outcome)
da_exp <- table(subset(outcomes, differentiation=="dopaminergic_neurons_exp" & varCategory=="all")$outcome)
da_pred <- table(subset(outcomes, differentiation=="dopaminergic_neurons_pred" & varCategory=="all")$outcome)

vec <- c(paste0("SN_S", paste0("\n(n=",sn[["Successful"]],")")),
         paste0("SN_F", paste0("\n(n=",sn[["Failed"]],")")),
         paste0("M_S", paste0("\n(n=",ma[["Successful"]],")")),
         paste0("M_F", paste0("\n(n=",ma[["Failed"]],")")),
         paste0("DAobs_S", paste0("\n(n=",da_exp[["Successful"]],")")),
         paste0("DAobs_F", paste0("\n(n=",da_exp[["Failed"]],")")),
         paste0("DApred_S", paste0("\n(n=",da_pred[["Successful"]],")")),
         paste0("DApred_F", paste0("\n(n=",da_pred[["Failed"]],")")))
names(vec) <- unique(outcomes$differentiation2)
outcomes$differentiation3 <- unname(vec[outcomes$differentiation2])

vec2 <- c("All variants","Synonymous variants",
          "Other","Deleterious variants",
          "Missense Pathogenic","Loss-of-function")
names(vec2) <- unique(outcomes$varCategory)
outcomes$varCategory <- unname(vec2[outcomes$varCategory])

my_comparisons <- list(c(unname(vec[1]), unname(vec)[2]),
                       c(unname(vec[3]), unname(vec[4])),
                       c(unname(vec[5]), unname(vec[6])),
                       c(unname(vec[7]), unname(vec[8])))

outcomes$varCategory <- factor(outcomes$varCategory, levels=unique(outcomes$varCategory))

outplot <- ggplot(outcomes, aes(x=differentiation3, y=mutBurden, col=outcome))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha=0.2, size=1)+
  facet_wrap(~varCategory, scales="free_y", ncol = 2)+
  scale_color_manual(name="Outcome",values=c('#E69F00', '#56B4E9'))+
  xlab("")+
  #xlab("DA=Dopaminergic neurons (obs=observed, pred=predicted), M=macrophages, SN=sensory neurons, n=Number of cell-lines")+
  ylab("Somatic mutational burden per cell-line")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
        axis.text.y=element_text(size=12, colour="black"),
        axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1, colour="black"),
        legend.position="none",
        strip.text.x = element_text(size = 12))



outplot <- outplot + 
  geom_signif(test="wilcox.test",
              comparisons=my_comparisons,
              map_signif_level = F,
              textsize=5,
              col="black",
              vjust=2)

pdf(file=paste0("figures/suppFigs/suppfig1A.pdf"))
outplot
dev.off()


###############
### SFig 1B ###
###############


runLogisticReg <- function(exomes_allinfo, cellType="sensory_neurons"){
  
  print(cellType)
  columns_to_test <- colnames(exomes_allinfo)[2:15]
  
  res <- sapply(columns_to_test, function(x){
    
    #print(x)
    mtrix <- exomes_allinfo[,match(c(x,cellType), colnames(exomes_allinfo))]
    mtrix <- mtrix[!(is.na(mtrix[,1]) | is.na(mtrix[,2])),]
    logistic <- glm(as.factor(mtrix[,2]) ~ mtrix[,1], data=mtrix, family="binomial")
    df <- data.frame(variantCat=x,
                     cType=cellType,
                     failed=sum(mtrix[,2]=="Failed"),
                     successful=sum(mtrix[,2]=="Successful"),
                     intercept=unname(summary(logistic)$coefficients[,"Estimate"][1]),
                     logOddsRatio=unname(summary(logistic)$coefficients[,"Estimate"][2]),
                     pVal=summary(logistic)$coefficients[,"Pr(>|z|)"][2])
    
  }, simplify=F)
  
  res <- do.call("rbind", res)
  rownames(res) <- NULL
  res
  
}


resLogistic <- rbind(runLogisticReg(diffmetrics, cellType="sensory_neurons"),
                     runLogisticReg(diffmetrics, cellType="macrophages"),
                     runLogisticReg(diffmetrics, cellType="dopaminergic_neurons_exp"),
                     runLogisticReg(diffmetrics, cellType="dopaminergic_neurons_pred"))




newvec <- c("SN","M","DAobs","DApred")
names(newvec) <- unique(resLogistic$cType)
resLogistic$cType <- unname(newvec[resLogistic$cType])

newvec <- c("All","Synonymous","Missense",
            "Deleterious","Missense Pathogenic","LoF",
            "Other","Dinucleotides","Only SNVs",
            "Indels","LoF Indels","Number of CNV",
            "CNV Length (different, Mb)","CNV Length (shared, Mb)")
names(newvec) <- unique(resLogistic$variantCat)
resLogistic$variantCat <- unname(newvec[resLogistic$variantCat])


resLogistic$variantCat <- factor(resLogistic$variantCat, levels=rev(unname(newvec)))
my_breaks <- c(0.001,0.01,0.05,0.1,1)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
logReg_heatmap <- ggplot(resLogistic, aes(x=cType, y=variantCat, fill=pVal)) + 
  geom_tile() +
  theme_bw()+
  theme(axis.text.x=element_text(vjust=0.5, size=12, colour="black"),
        axis.text.y=element_text(size=12, colour="black"),
        axis.title=element_text(size=14, colour="black"),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=11, colour="black"))+
  scale_fill_gradientn(name="Raw P-value",
                       colours=myPalette(100),
                       trans = "log10",
                       limits=c(0.001,1),
                       labels=my_breaks,
                       breaks=my_breaks)+
  xlab("iPSC-derived cell-type")+
  ylab("Variant category")

pdf(file=paste0("figures/suppFigs/suppfig1B.pdf"), width = 6, height = 4)
plot(logReg_heatmap)
dev.off()

###############
### SFig 1C ###
###############

allvariants <- readRDS("outputTabs/allvariants_v2.RDS")

### UCSC ###
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cdsbytx <- cdsBy(txdb, by="tx")
txbygene <- transcriptsBy(txdb, by="gene")
txGene <- data.frame(gene=names(unlist(txbygene)),
                     tx=unlist(txbygene)$tx_id)
cds_length <- sapply(width(cdsbytx), function(x) sum(x))
txGene$cds_tx_length <- NA
mask <- txGene$tx %in% names(cds_length)
txGene$cds_tx_length[mask] <- unname(cds_length[match(txGene$tx[mask], names(cds_length))])
maxCdsxGene <- split(txGene$cds_tx_length, txGene$gene)
maxCdsxGene <- maxCdsxGene[!sapply(maxCdsxGene, function(x) all(is.na(x)))]
maxCdsxGene <- sapply(maxCdsxGene, function(x) max(x[!is.na(x)]))
maxCdsxGene <- data.frame(gene=names(maxCdsxGene),
                          max_cds_length=unname(maxCdsxGene))
symbolEntrez <- AnnotationDbi::select(org.Hs.eg.db, keys = maxCdsxGene$gene, columns = "SYMBOL", keytype = "ENTREZID")
symbolEntrez <- symbolEntrez[complete.cases(symbolEntrez),]
maxCdsxGene$symbol <- symbolEntrez[match(maxCdsxGene$gene, symbolEntrez$ENTREZID),]$SYMBOL

## Gene universe using 87 gff3 version of the GRCh37 human genome
test <- import("inputTabs/Homo_sapiens.GRCh37.87.gff3.gz")
cdsbytx <- test[test$type=="CDS",]
cdsbytx$transcript_id <- gsub(".+\\:","", cdsbytx$Parent)
names(cdsbytx) <- cdsbytx$transcript_id
cdsbytx <- split(cdsbytx, names(cdsbytx))

txbygene <- test[test$type=="transcript" | test$type=="mRNA"]
txbygene <- txbygene[txbygene$biotype=="protein_coding"]
txbygene$gene_id <- gsub(".+\\:","", txbygene$Parent)
names(txbygene) <- txbygene$gene_id
txbygene <- split(txbygene, names(txbygene))

##genomic coordinates
genes <- test[test$type=="gene"]
txGene <- data.frame(gene=rep(names(txbygene), unname(elementNROWS(txbygene))),
                     tx=unlist(txbygene)$transcript_id)
cds_length <- sapply(width(cdsbytx), function(x) sum(x))
txGene$cds_tx_length <- NA
mask <- txGene$tx %in% names(cds_length)
txGene$cds_tx_length[mask] <- unname(cds_length[match(txGene$tx[mask], names(cds_length))])
maxCdsxGene <- split(txGene$cds_tx_length, txGene$gene)
maxCdsxGene <- maxCdsxGene[!sapply(maxCdsxGene, function(x) all(is.na(x)))]
maxCdsxGene <- sapply(maxCdsxGene, function(x) max(x[!is.na(x)]))
maxCdsxGene <- data.frame(gene=names(maxCdsxGene),
                          max_cds_length=unname(maxCdsxGene))
maxCdsxGene <- maxCdsxGene[!is.na(match(maxCdsxGene$gene, genes$gene_id)),]
maxCdsxGene$symbol <- genes[match(maxCdsxGene$gene, genes$gene_id),]$Name
txMaxCds <- sapply(maxCdsxGene$gene, function(x) subset(txGene, txGene$gene==x)[which.max(subset(txGene, txGene$gene==x)$cds_tx_length),]$tx)
stopifnot(all(maxCdsxGene$gene==names(txMaxCds)))
maxCdsxGene$tx <- unname(txMaxCds)
cdsbytx <- cdsbytx[match(maxCdsxGene$tx, names(cdsbytx))]
genes <- genes[genes$biotype=="protein_coding"]

### The variants from the genes that are not accounted by the used annotation in downstream analysis are filtered-out
##Remove variants in genes not present in the annotation

allvariants <- allvariants[allvariants$gene_name %in% genes$Name,]
##Remove variants that do not fall within the genes of the annotation
tmpwes <- GRanges(allvariants$chrom, IRanges(allvariants$pos, width=nchar(allvariants$ALT)), strand="*")
hits2 <- findOverlaps(tmpwes, genes)
stopifnot(length(unique(queryHits(hits2)))==length(tmpwes))


##hg19 RepeatMasker

##Requires bedops installation in bin path
system("export PATH=/soft/bedops/bin:$PATH")

## Most up-to-date version
system('wget -qO- ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/rmsk.txt.gz | gunzip -c - | awk -v OFS="\t" "{ print $6, $7, $8, $12, $11, $10 }" - | sort-bed - > inputTabs/rmsk.bed')

## Repetitive elements
repMasker <- read.table("inputTabs/rmsk.bed")
colnames(repMasker) <- c("genoName","genoStart","genoEnd","repClass","repName","strand")

discarded <- c("DNA?","LINE?","LTR?","RC","RC?","SINE?","Unknown","Unknown?","Other")
repMasker <- repMasker[!repMasker$repClass %in% discarded,]

tmpwes <- GRanges(seqnames=allvariants$chrom, IRanges(allvariants$pos, width=nchar(allvariants$ALT)), strand="*")
seqlevelsStyle(tmpwes) <- "UCSC"
configs <- levels(seqnames(tmpwes))

repMasker <- repMasker[repMasker$genoName %in% configs,]
repMaskGR <- GRanges(seqnames=repMasker$genoName, IRanges(start=repMasker$genoStart, end=repMasker$genoEnd), strand="*")
repMaskGR$repClass <- repMasker$repClass
seqlevelsStyle(repMaskGR) <- "NCBI"

##anyHit (proximity to repetitive element)

repElementGenerator <- function(genes, repMaskGR, typeElement="All Rep Elements"){
  
  seqlevelsStyle(tmpwes) <- "NCBI"
  hitsGenes <- findOverlaps(tmpwes, genes)
  genesMut <- genes[unique(subjectHits(hitsGenes))]
  genesUnMut <- genes[-unique(subjectHits(hitsGenes))]
  
  if (typeElement=="All Rep Elements"){
    
    nearestDistMut <- distanceToNearest(genesMut,repMaskGR)
    nearestDistUnMut <- distanceToNearest(genesUnMut,repMaskGR)
    
    ## mutated
    dfNearestMut <- data.frame(genes=genesMut$Name,
                               distance=NA,
                               element=NA,
                               geneType="Mutated")
    dfNearestMut$distance[queryHits(nearestDistMut)] <- mcols(nearestDistMut)$distance
    dfNearestMut$element[queryHits(nearestDistMut)] <- repMaskGR[subjectHits(nearestDistMut),]$repClass
    
    ##unmutated
    dfNearestUnMut <- data.frame(genes=genesUnMut$Name,
                                 distance=NA,
                                 element=NA,
                                 geneType="Unmutated")
    dfNearestUnMut$distance[queryHits(nearestDistUnMut)] <- mcols(nearestDistUnMut)$distance
    dfNearestUnMut$element[queryHits(nearestDistUnMut)] <- repMaskGR[subjectHits(nearestDistUnMut),]$repClass
    
    dfNearest <- rbind(dfNearestMut, dfNearestUnMut)
    
    return(dfNearest)
    
  } else {
    
    repMaskGR_sub <- subset(repMaskGR, repClass==typeElement)
    
    nearestDistMut <- distanceToNearest(genesMut,repMaskGR_sub)
    nearestDistUnMut <- distanceToNearest(genesUnMut,repMaskGR_sub)
    
    ## mutated
    dfNearestMut <- data.frame(genes=genesMut$Name,
                               distance=NA,
                               element=NA,
                               geneType="Mutated")
    dfNearestMut$distance[queryHits(nearestDistMut)] <- mcols(nearestDistMut)$distance
    dfNearestMut$element[queryHits(nearestDistMut)] <- repMaskGR_sub[subjectHits(nearestDistMut),]$repClass
    
    ##unmutated
    dfNearestUnMut <- data.frame(genes=genesUnMut$Name,
                                 distance=NA,
                                 element=NA,
                                 geneType="Unmutated")
    dfNearestUnMut$distance[queryHits(nearestDistUnMut)] <- mcols(nearestDistUnMut)$distance
    dfNearestUnMut$element[queryHits(nearestDistUnMut)] <- repMaskGR_sub[subjectHits(nearestDistUnMut),]$repClass
    
    dfNearest <- rbind(dfNearestMut, dfNearestUnMut)
    
    return(dfNearest)
    
  }
}


repelements <- c("All Rep Elements",unique(repMaskGR$repClass)[-length(unique(repMaskGR$repClass))])

dfNearest <- sapply(repelements, function(y){
  repElementGenerator(genes, repMaskGR, typeElement=y)
}, simplify=F)

## summary statistics

dfSumm <- sapply(dfNearest, function(x){
  
  x <- subset(x, !is.na(distance))
  pvaltest <- t.test(subset(x, geneType=="Mutated")$distance,
                     subset(x, geneType=="Unmutated")$distance)$p.value
  
  data.frame(pval=pvaltest,
             mean.dist.mut=format(mean(subset(x, geneType=="Mutated")$distance), scientific=F),
             mean.dist.unmut=format(mean(subset(x, geneType=="Unmutated")$distance),scientific=F),
             perc.zero.mut=signif(dim(subset(x, geneType=="Mutated" & distance==0))[1]*100/dim(subset(x, geneType=="Mutated"))[1],2),
             perc.zero.unmut=signif(dim(subset(x, geneType=="Unmutated" & distance==0))[1]*100/dim(subset(x, geneType=="Unmutated"))[1],2))
  
}, simplify=F)


dfSumm <- do.call("rbind", dfSumm)
dfSumm$test <- rownames(dfSumm)
rownames(dfSumm) <- NULL

dfSumm <- dfSumm[,colnames(dfSumm)[c(length(colnames(dfSumm)),1: (length(colnames(dfSumm))-1))]]
dfSumm$pvalAdj <- p.adjust(dfSumm$pval, "BH")

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

dfSumm$signif <- pvalConverter(dfSumm$pvalAdj)

dfNearestAll <- do.call("rbind", dfNearest)
dfNearestAll$type <- gsub("\\..+","",rownames(dfNearestAll))
rownames(dfNearestAll) <- NULL
dfNearestAll <- dfNearestAll[!is.na(dfNearestAll$distance),]
dfNearestAll <- dfNearestAll[dfNearestAll$element!="RNA",]


dfNearestAll$type <- factor(dfNearestAll$type, levels=rev(c("All Rep Elements","Satellite","LTR","LINE","SINE","Simple_repeat",
                                                            "Low_complexity","DNA","rRNA","snRNA","scRNA","srpRNA","tRNA")))

repPlot <- ggplot(dfNearestAll, aes(y=type, x=distance, fill=geneType))+
  theme_bw()+
  geom_boxplot(outlier.size=0.05)+
  scale_fill_brewer(name="T-test", palette="Set2")+
  xlab("Proximity to repetitive element (in bp)")+
  ylab("Repetitive element")+
  scale_x_continuous(trans="log1p",
                     breaks=c(0,1,10,100,1000,10000,100000,1000000,25000000),
                     labels = scales::label_comma(accuracy = 1),
                     expand = expansion(add = c(0, 0.5)))+
  geom_text(data = dfSumm, aes(y = test, x = 32000000, label = signif),
            size=5.5, inherit.aes=F)+
  theme(axis.text.x=element_text(size=10, colour="black"),
        axis.text.y=element_text(size=13, colour="black"),
        axis.title.y=element_text(size=16, colour="black"),
        axis.title.x=element_text(size=16, colour="black", hjust=0.5, vjust=0),
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=13, face="bold"),
        legend.position="top")+
  stat_summary(fun=mean, geom="point", shape=18, size=3, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))

pdf(file=paste0("figures/suppFigs/suppfig1C.pdf"), width = 8, height = 4)
logReg_heatmap
dev.off()

###############
### SFig 1D ###
###############

seqlevelsStyle(tmpwes) <- "NCBI"
hitsGenes <- findOverlaps(tmpwes, genes)
genesMut <- genes[unique(subjectHits(hitsGenes))]
genesUnMut <- genes[-unique(subjectHits(hitsGenes))]

allGenes <- rbind(data.frame(genes=genesMut$Name,
                             genLength=width(genesMut),
                             type="Mutated"),
                  data.frame(genes=genesUnMut$Name,
                             genLength=width(genesUnMut),
                             type="Unmutated"))


lenPlot <- ggplot(allGenes, aes(x=type, y=genLength))+
  theme_bw()+
  geom_violin()+
  geom_boxplot(outlier.size=0.75, width=0.4)+
  scale_fill_brewer(name="T-test", palette="Set2")+
  ylab("Gene length (in bp)")+
  xlab("")+
  scale_y_continuous(trans="log1p",
                     breaks=c(0,1,10,100,1000,10000,100000,200000,500000,1000000,2000000),
                     labels = scales::label_comma(accuracy = 1),
                     expand = expansion(add = c(0, 1)))+
  theme(axis.text.x=element_text(size=16, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.title.y=element_text(size=16, colour="black", vjust=-1.5))+
  geom_signif(comparisons = list(c("Unmutated", "Mutated")),
              map_signif_level = TRUE,
              textsize = 8)+
  stat_summary(fun=mean, geom="point", shape=18, size=4, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))


pdf(file=paste0("figures/suppFigs/suppfig1D.pdf"), width = 4, height = 6)
plot(lenPlot)
dev.off()


###############
### SFig 1E ###
###############


directoryRes <- "outputTabs/"
filesIN <- dir(directoryRes)[grepl("_PseudoBulkMeanExpPerLine.RDS", dir(directoryRes))]
dfPseudoBulk <- sapply(filesIN, function(x){
  
  tmp <- readRDS(paste0(directoryRes,x))
  tmp$timePoint <- sapply(strsplit(x, "_"), function(y) y[1])
  tmp$category <- NA
  tmp[tmp$genes %in% genesMut$Name,]$category <- "Mutated"
  tmp[tmp$genes %in% genesUnMut$Name,]$category <- "Unmutated"
  tmp <- tmp[!is.na(tmp$category),]
  tmp
  
}, simplify=F)

dfPseudoBulk <- do.call("rbind", dfPseudoBulk)
rownames(dfPseudoBulk) <- NULL


exprPlot <- ggplot(dfPseudoBulk, aes(x=category, y=AvgExp_PseudoBulkMeanExpPerLine))+
  theme_bw()+
  geom_violin()+
  geom_boxplot(outlier.size=0.5, width=0.05)+
  scale_fill_brewer(name="T-test", palette="Set2")+
  ylab("Average pseudobulk expression per gene")+
  xlab("")+
  facet_wrap(~timePoint)+
  scale_y_continuous(trans="log1p",
                     breaks=c(0,0.1,0.2,0.5,1,2,5),
                     expand = expansion(add = c(0.05, 0.2)))+
  theme(axis.text.x=element_text(size=12, colour="black", angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=12, colour="black"),
        axis.title=element_text(size=16, colour="black"))+
  geom_signif(comparisons = list(c("Unmutated", "Mutated")),
              map_signif_level = TRUE,
              textsize = 8,
              test = "t.test")+
  stat_summary(fun=mean, geom="point", shape=18, size=4, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))

pdf(file=paste0("figures/suppFigs/suppfig1E.pdf"), width = 4, height = 6)
plot(exprPlot)
dev.off()



###############
### SFig 1F ###
###############

new_genesMut <- readRDS("outputTabs/genesMutAnnotRepli.RDS")
new_genesUnMut <- readRDS("outputTabs/genesUnMutAnnotRepli.RDS")

mutState <- rbind(data.frame(gene=new_genesMut$Name,
                             category="Mutated",
                             avgRepliSeq=new_genesMut$avgRepliSeq),
                  data.frame(gene=new_genesUnMut$Name,
                             category="Unmutated",
                             avgRepliSeq=new_genesUnMut$avgRepliSeq))


repliPlot <- ggplot(mutState, aes(x=category, y=avgRepliSeq))+
  theme_bw()+
  geom_violin()+
  geom_boxplot(outlier.size=0.75, width=0.2)+
  scale_fill_brewer(name="T-test", palette="Set2")+
  ylab("ENCODE Replication time (Repli-Seq)")+
  xlab("")+
  theme(axis.text.x=element_text(size=16, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.title.y=element_text(size=16, colour="black", vjust=1.5))+
  geom_signif(comparisons = list(c("Unmutated", "Mutated")),
              map_signif_level = FALSE,
              textsize = 8)+
  stat_summary(fun=mean, geom="point", shape=18, size=4, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  scale_y_continuous(expand=expansion(add = c(2, 8)),
                     breaks=seq(0,100,20))

## suppFig1F
pdf(file=paste0("figures/suppFigs/suppfig1F.pdf"), width = 4, height = 6)
plot(repliPlot)
dev.off()


###############
### SFig 1G ###
###############

set.seed(123)
mutTab_logReg <- readRDS("analysis/outputTabs/mutTab_logReg2.RDS")
pathToFile <- "analysis/outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
synMutBurden <- readRDS(paste0(pathToFile, "mutBurden_synonymous.RDS"))
delMutBurden <- readRDS(paste0(pathToFile, "mutBurden_del.RDS"))

mutTab_logReg <- mutTab_logReg[,-match(c("sensoryNeurons_outcome","macrophages_outcome", "BCORpos","phenotype"), colnames(mutTab_logReg))]
mutTab_logReg$bcorPtv <- unname(ptvMutBurden["BCOR",])[match(mutTab_logReg$cline, names(ptvMutBurden["BCOR",]))]
mutTab_logReg$bcorSyn <- unname(synMutBurden["BCOR",])[match(mutTab_logReg$cline, names(synMutBurden["BCOR",]))]
mutTab_logReg$bcorDel <- unname(delMutBurden["BCOR",])[match(mutTab_logReg$cline, names(delMutBurden["BCOR",]))]
mutTab_logReg_long <- as.data.frame(mutTab_logReg %>% pivot_longer(-c(1:3), names_to="varClass", values_to="nMut"))
mutTab_logReg_long <- as.data.frame(mutTab_logReg_long %>% pivot_longer(-c(1,4,5), names_to="dataset", values_to="Outcome"))

mutTab_logReg_pred <- subset(mutTab_logReg_long, dataset=="NeuroSeq_pred")
mutTab_logReg_pred <- subset(mutTab_logReg_pred, !is.na(Outcome))
vecVarClass <- c("LoF variants", "Synonymous", "Deleterious")
names(vecVarClass) <- unique(mutTab_logReg_pred$varClass)
mutTab_logReg_pred$varClass <- unname(vecVarClass[mutTab_logReg_pred$varClass])

sizeGroups <- table(mutTab_logReg_pred$Outcome)/unique(table(mutTab_logReg_pred$cline))
vecOutcome <- paste0(names(sizeGroups),"\n","n=", unname(sizeGroups))
names(vecOutcome) <- rev(unique(mutTab_logReg_pred$Outcome))
mutTab_logReg_pred$Outcome <- unname(vecOutcome[mutTab_logReg_pred$Outcome])

my_comparison <- list(unique(mutTab_logReg_pred$Outcome))

predOut <- ggplot(mutTab_logReg_pred, aes(x=Outcome, y=nMut))+
  geom_jitter(height = 0.15, width=0.35, alpha=0.25, size=5)+
  facet_wrap(~varClass)+
  theme_bw()+
  ggtitle("Dopaminergic Neurons: Predicted outcome")+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        axis.text.x=element_text(size=14, colour="black"),
        axis.text.y=element_text(size=14, colour="black"),
        strip.text.x = element_text(size = 14))+
  xlab("")+ylab("")+
  stat_compare_means(
    comparisons = my_comparison,
    label = "wilcox.test", step.increase = 0.07, size=5)+
  coord_cartesian(ylim=c(0,3.5))

mutTab_logReg_obs <- subset(mutTab_logReg_long, dataset=="NeuroSeq_obs")
mutTab_logReg_obs <- subset(mutTab_logReg_obs, !is.na(Outcome))
vecVarClass <- c("LoF variants", "Synonymous", "Deleterious")
names(vecVarClass) <- unique(mutTab_logReg_obs$varClass)
mutTab_logReg_obs$varClass <- unname(vecVarClass[mutTab_logReg_obs$varClass])

sizeGroups <- table(mutTab_logReg_obs$Outcome)/unique(table(mutTab_logReg_obs$cline))
vecOutcome <- paste0(names(sizeGroups),"\n","n=", unname(sizeGroups))
names(vecOutcome) <- rev(unique(mutTab_logReg_obs$Outcome))
mutTab_logReg_obs$Outcome <- unname(vecOutcome[mutTab_logReg_obs$Outcome])

my_comparison <- list(unique(mutTab_logReg_obs$Outcome))

obsOut <- ggplot(mutTab_logReg_obs, aes(x=Outcome, y=nMut))+
  geom_jitter(height = 0.15, width=0.35, alpha=0.25, size=5)+
  facet_wrap(~varClass)+
  theme_bw()+
  ggtitle("Dopaminergic Neurons: Observed outcome")+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        axis.text.x=element_text(size=14, colour="black"),
        axis.text.y=element_text(size=14, colour="black"),
        strip.text.x = element_text(size = 14))+
  xlab("")+ylab("")+
  stat_compare_means(
    comparisons = my_comparison,
    label = "wilcox.test", step.increase = 0.07, size=5)+
  coord_cartesian(ylim=c(0,3.5))


bcorPlot <- ggarrange(predOut, obsOut,
                      ncol=1, nrow=2)

bcorPlot <- annotate_figure(bcorPlot,
                            left=text_grob("Number of BCOR variants per iPSC cell-line", rot=90, size=11, face="bold"))
bcorPlot <- annotate_figure(bcorPlot)

pdf(file=paste0("figures/suppFigs/suppfig1G.pdf"))
plot(bcorPlot)
dev.off()


###############
### SFig 1H ###
###############

directoryRes <- "outputTabs/"

filesIN <- dir(directoryRes)[grepl("_pseudoBulkPerLine.RDS", dir(directoryRes))]

dfPseudoBulk <- sapply(filesIN, function(x){
  
  tmp <- readRDS(paste0(directoryRes,x))
  tmp$timePoint <- sapply(strsplit(x, "_"), function(y) y[1])
  tmp$label1 <- ifelse(tmp$BCORpositive==TRUE & tmp$DNoutcome=="Failed", "Failed/BCOR+","Unmutated")
  tmp$label2 <- ifelse(tmp$BCORpositive==TRUE & tmp$DNoutcome=="Failed", "Failed/BCOR+",NA)
  tmp$label2[tmp$BCORpositive==FALSE & tmp$DNoutcome=="Failed"] <- "Failed/BCOR-"
  tmp$label2[tmp$BCORpositive==FALSE & tmp$DNoutcome=="Successful"] <- "Successful"
  
  tmp
  
}, simplify=F)

dfPseudoBulk <- do.call("rbind", dfPseudoBulk)
rownames(dfPseudoBulk) <- NULL


my_comparisonsTest <- list(c("Failed/BCOR+","Unmutated"))

exprLevel <- ggplot(dfPseudoBulk,
                    aes(x=label1, y=log10(meanExpPseudoBulk)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_jitter(width = 0.2), alpha=0.2, size=2)+
  theme_bw()+
  xlab("")+
  facet_wrap(~timePoint)+
  theme(axis.text.x=element_text(angle=60, size=12, vjust=0.5),
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        plot.title=element_text(hjust=0.5, face="bold", size=13),
        legend.title=element_text(face="bold"),
        strip.text = element_text(size = 14),
        legend.position="top")+
  ggtitle("Pseudobulk expression per iPSC line and time point")+
  stat_compare_means(
    comparisons = my_comparisonsTest,
    label = "wilcox.test", step.increase = 0.07, size=4
  )+
  ylab("LogNormalised BCOR expression (in log10)")


pdf(file="figures/suppFigs/suppfig1H.pdf")
plot(exprLevel)
dev.off()





