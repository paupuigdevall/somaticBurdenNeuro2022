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

source("functionsToImport.R")

##############
### Fig 4A ###
##############

pathTODE <- "analysis/outputTabs/DEsinglecell/"
summResults <- paste0(pathTODE, list.files(path=pathTODE,
                                           pattern = "resultsDEFailVsSucc2.+.RDS")) %>%
  map(readRDS) %>% 
  bind_rows()

heatmap_df <- summResults
heatmap_df <- heatmap_df[,match(c("annot","timepoint","numDE","numFailedDonors","numSuccessfulDonors","numFailedCells","numSuccessfulCells"), colnames(heatmap_df))]
heatmap_df <- heatmap_df[!duplicated(heatmap_df),]

sfig4a <- ggplot(data=heatmap_df, aes(x=numDE, y=numFailedCells, col=timepoint))+
  geom_point(alpha=0.5, size=3)+
  xlab("Number of DE genes")+
  ylab("Number of cells")+
  scale_y_continuous(breaks=seq(0,80000, 20000), labels = comma)+
  geom_smooth(method="lm", col="grey40", alpha=0.25, size=0.5)+
  stat_cor(inherit.aes=F,
           aes(x=numDE, y=numFailedCells),
           label.y = 75000, label.x=300, method="pearson", digits=2, size=4.25, hjust=0.5)+
  theme_bw()+
  ggtitle("Failed")+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=15),
        axis.title=element_text(size=13),
        axis.text=element_text(size=12),
        legend.position="top",
        legend.title=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(name="Timepoint",
                     values=wes_palette(n = length(unique(heatmap_df$timepoint)), "GrandBudapest1"))+
  coord_cartesian(ylim=c(0,80000))+
  xlab("")+ylab("")

sfig4b <- ggplot(data=heatmap_df, aes(x=numDE, y=numSuccessfulCells, col=timepoint))+
  geom_point(alpha=0.5, size=3)+
  xlab("Number of DE genes")+
  ylab("Number of cells")+
  scale_y_continuous(breaks=seq(0,80000, 20000), labels = comma)+
  geom_smooth(method="lm", col="grey40", alpha=0.25, size=0.5)+
  stat_cor(inherit.aes=F,
           aes(x=numDE, y=numSuccessfulCells),
           label.y = 75000, label.x=300, method="pearson", digits=2, size=4.25, hjust=0.5)+
  theme_bw()+
  ggtitle("Successful")+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=15),
        axis.title=element_text(size=13),
        axis.text=element_text(size=12))+
  scale_color_manual(name="Timepoint",
                     values=wes_palette(n = length(unique(heatmap_df$timepoint)), "GrandBudapest1"))+
  coord_cartesian(ylim=c(0,80000))+
  guides(col=F)+
  xlab("")+ylab("")

sfig4 <- ggarrange(sfig4a, sfig4b,
                   ncol=2,nrow=1, common.legend = TRUE)

pdf(file=paste0("figures/suppFigs/suppfig4A.pdf"), width = 6, height = 4)
sfig4
dev.off()


##############
### Fig 4B ###
##############


library(GSEABase)
library(fgsea)

## evaluate QC (how many genes do we lose from gene Universe to DEobject in Seurat)

pathTODE <- "analysis/outputTabs/DEsinglecell/"
cosmiclogFC <- paste0(pathTODE, list.files(path=pathTODE,
                                           pattern = "allDEtable2_")) %>%
  map(readRDS) %>% 
  bind_rows()


logFC_universe <- sapply(unique(cosmiclogFC$annot), function(x){
  tmp <- subset(cosmiclogFC, annot==x)
  tmp <- setNames(tmp$log2FC, tmp$geneId)
  tmp
}, simplify=F)


#Download MSigDb Hallmarks for cancer (version 7.4)
gsea_gset <- getGmt("analysis/inputTabs/h.all.v7.4.symbols.gmt")

pathways <- sapply(logFC_universe, function(y){
  gsea_list <- sapply(gsea_gset, function(x){
    unique(names(y)[names(y) %in% geneIds(x)])
  }, simplify=F)
  names(gsea_list) <- names(gsea_gset)
  gsea_list <- gsea_list[sapply(gsea_list, function(x) length(x))>10]
  gsea_list
  
}, simplify=F)


fgseaResList <- sapply(seq(pathways), function(x){
  set.seed(123)
  fgseaRes <- fgseaMultilevel(pathways = pathways[[x]], stats=logFC_universe[[x]],
                              maxSize=500, nPermSimple = 10000, eps=0)
  fgseaRes
  
}, simplify=F)

names(fgseaResList) <- unique(cosmiclogFC$annot)

#low_bound <- quantile(unlist(sapply(fgseaResList, function(x) x[x$padj<0.05,]$NES, simplify=F)),0.05)
#high_bound <- quantile(unlist(sapply(fgseaResList, function(x) x[x$padj<0.05,]$NES, simplify=F)),0.95)
fgseaResList_sub=sapply(fgseaResList, function(x) x[x$padj<0.05 & (x$NES>1.5 | x$NES<=-1.5),]$pathway, simplify=F)
pathwaysToShow=unique(unlist(fgseaResList_sub))
#names(pathwaysToShow) <- c("NA","NA","Hs","Hs","")

fgseaResList_sub <- sapply(fgseaResList, function(x){
  x <- x[x$pathway %in% pathwaysToShow,]
  x[x$padj>=0.05,]$NES <- NA
  x
}, simplify=F)

heatmapCosmic <- sapply(seq(fgseaResList_sub), function(y){
  
  fgseaResList_sub[[y]]$annot <- names(fgseaResList_sub)[y]
  fgseaResList_sub[[y]][,c("pathway","padj","NES","annot")]
  
}, simplify=F)

heatmapCosmic <- do.call("rbind", heatmapCosmic)
rownames(heatmapCosmic) <- NULL
heatmapCosmic <- as.data.frame(heatmapCosmic)


maxNES <- ceiling(max(heatmapCosmic$NES, na.rm=T))
minNES <- floor(min(heatmapCosmic$NES, na.rm=T))

heatmapCosmic$pathway <- gsub("_"," ", gsub("HALLMARK_","",heatmapCosmic$pathway))
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))

suppFig4b <- ggplot(heatmapCosmic,
                      aes(y=pathway, x=annot, fill=NES))+
  geom_tile(colour = "black")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=14),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(angle=90, size=15),
        plot.title=element_text(hjust=0.5, face="bold", size=13),
        legend.title=element_text(face="bold", size=14),
        legend.text=element_text(size=14))+
  xlab("")+
  ylab("")+
  #ggtitle("Oncogenic signature enrichment")+
  scale_fill_gradientn(name="Normalised\nenrichment\nscore (NES)",
                       colours=myPalette(100),
                       limits=c(minNES,maxNES),
                       labels=c(minNES,-1.5,1.5,maxNES),
                       breaks=c(minNES,-1.5,1.5,maxNES),
                       na.value = 'grey90')+
  guides(fill=guide_legend(order=1))

pdf(file="figures/suppFigs/suppfig4B.pdf", width=8, height = 6)
plot(suppFig4b)
dev.off()



