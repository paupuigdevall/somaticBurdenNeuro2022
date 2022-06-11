
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

source("functionsToImport.R")

summData <- readRDS("analysis/outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)

#############
### fig5A ###
#############

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


fig5a <- ggplot(summData, aes(annot, cfrac, col=outlier))+
  geom_jitter(width=0.15, size=0.3, alpha=0.5)+
  facet_wrap(~tp)+
  ylab("Cell-type fraction")+
  xlab("")+
  theme_bw()+
  scale_y_continuous(breaks=seq(0,1,0.2), labels=seq(0,1,0.2))+
  #ggtitle("Cell-type fraction outlier identification per donor ")+
  #xlab("Cell types")+ 
  #theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+
  theme(axis.text.x=element_text(angle=90, size=11, hjust=1, vjust=0.5),
        axis.title.y=element_text(size=12),
        legend.position="top",
        legend.text=element_text(size=12),
        strip.text = element_text(size=13))+
  scale_colour_manual(name="",labels=c("Non-Outlier","Outlier"),values=c("grey60","red"))+
  guides(alpha=F, size=F, color = guide_legend(override.aes = list(size=5)))
  

pdf(file="figures/mainFigs/figure5A.pdf", width=8, height = 4)
plot(fig5a)
dev.off()


#############
### fig5B ###
#############



#############
### fig5B ###
#############

summData$pool <- sapply(strsplit(summData$donor_extended, "/"), function(x) x[length(x)])
diffmetrics <- read.table("suppTabs/suppTable1.txt",
                          header=T)
outlier_analysis <- data.frame(donor_extended=unique(summData$donor_extended),
                               outlier_alltp_z2=NA)
avail <- sapply(unique(summData$donor_extended), function(x){ paste(unique(subset(summData, donor_extended==x)$tp), collapse="-")}, simplify=T)
outlier_analysis$avail <- unname(avail[match(outlier_analysis$donor_extended, names(avail))])

#all outliers (all tp)
outliers_z2 <- unique(summData[summData$outlier==TRUE,]$donor_extended)
outlier_analysis$outlier_alltp_z2 <- outlier_analysis$donor_extended %in% outliers_z2
outlier_analysis$donor_id <- sapply(strsplit(outlier_analysis$donor_extended,"/"), function(x) paste(x[-length(x)],collapse="/"))
outlier_analysis$allMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$all
outlier_analysis$delMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$deleterious

outlierLong <- outlier_analysis[,c("donor_extended","outlier_alltp_z2","allMutBurden","delMutBurden")]
outlierLong <- as.data.frame(outlierLong %>% pivot_longer(-c(1,2), names_to="MutClass", values_to="Burden"))
outlierLong <- subset(outlierLong, !is.na(Burden))

facet_labeller = labeller(
  MutClass = c("allMutBurden" = "All variants",
               "delMutBurden" = "Deleterious"))
outlierLong$MutClass <- factor(outlierLong$MutClass, levels=unique(outlierLong$MutClass))

vec4 <- c("Outlier","Non-outlier")
names(vec4) <- c("TRUE","FALSE")
outlierLong$outlier_alltp_z2 <- unname(vec4[as.character(outlierLong$outlier_alltp_z2)])

my_comparisons <- list(c("Outlier","Non-outlier"))

figure5b_row2 <- ggplot(outlierLong, aes(x=outlier_alltp_z2, y=Burden))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.25, alpha=0.5, aes(col=outlier_alltp_z2))+
  facet_wrap(~MutClass, scales="free_y", labeller = facet_labeller)+
  xlab("")+
  theme_bw()+
  ylab("Mutational burden")+
  theme(strip.text =element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13))+
  geom_signif(test="wilcox.test",
              comparisons=my_comparisons,
              map_signif_level = T,
              textsize=4.5,
              col="black")+
  scale_colour_manual(values=c("black","red"))+
  guides(col=FALSE)

facet_bounds <- read.table(header=TRUE,
                           text=                           
                             "MutClass ymin  ymax  breaks
                           allMutBurden  0 300 50
                           delMutBurden 0 100 20",
                           stringsAsFactors=FALSE)

ff <- with(facet_bounds,
           data.frame(Burden=c(ymin,ymax),
                      MutClass=c(MutClass,MutClass)))

figure5b_row2 <- figure5b_row2 + geom_point(data=ff, x=NA)


q <- ggplot_build(figure5b_row2)

mask_panel <- q$data[[3]]$PANEL==1

q$data[[3]]$y[mask_panel] <- q$data[[3]]$y[mask_panel]-40
q$data[[3]]$yend[mask_panel] <- q$data[[3]]$yend[mask_panel]-40

q$data[[3]]$y[!mask_panel] <- q$data[[3]]$y[!mask_panel]-20
q$data[[3]]$yend[!mask_panel] <- q$data[[3]]$yend[!mask_panel]-20

figure5b_row2 <- ggplot_gtable(q)



#Outlier analysis
outliersOnly <- subset(summData, outlier==TRUE)
outAnalysis <- sapply(unique(outliersOnly$donor_extended), function(x) {
  tmp <- subset(outliersOnly, donor_extended==x)
  data.frame(donor_extended=x,
             nTotalOutliers=dim(tmp)[1],
             availTP11=dim(subset(summData, tp=="D11" & donor_extended==x))[1]>0,
             availTP30=dim(subset(summData, tp=="D30" & donor_extended==x))[1]>0,
             availTP52=dim(subset(summData, tp=="D52" & donor_extended==x))[1]>0,
             nOutliersTP11=dim(subset(tmp, tp=="D11"))[1],
             nOutliersTP30=dim(subset(tmp, tp=="D30"))[1],
             nOutliersTP52=dim(subset(tmp, tp=="D52"))[1],
             ctypesTP11=paste(subset(tmp, tp=="D11")$annot,collapse=","),
             ctypesTP30=paste(subset(tmp, tp=="D30")$annot,collapse=","),
             ctypesTP52=paste(subset(tmp, tp=="D52")$annot,collapse=","))
  
}, simplify=F)

outAnalysis <- do.call("rbind", outAnalysis)
rownames(outAnalysis) <- NULL

##focusing on those lines in the same pool that we have info on the three timepoints
outAnalysis$availOverall <- paste0(outAnalysis$availTP11,"-", outAnalysis$availTP30,"-", outAnalysis$availTP52)
outAnalysis_onlyavail <- subset(outAnalysis, availOverall=="TRUE-TRUE-TRUE")

compute_mean_se <- function(vector, tp="D11"){
  error <- qt(0.95,df=length(vector)-1)*sd(vector)/sqrt(length(vector))
  res <- data.frame(tp=tp,
                    mean=mean(vector),
                    error=error)
  return(res)
  
}
df_nOut <- rbind(compute_mean_se(outAnalysis_onlyavail$nTotalOutliers, tp="All TP"),
                 compute_mean_se(outAnalysis_onlyavail$nOutliersTP11, tp="D11"),
                 compute_mean_se(outAnalysis_onlyavail$nOutliersTP30, tp="D30"),
                 compute_mean_se(outAnalysis_onlyavail$nOutliersTP52, tp="D52"))

df_nOut$tp <- factor(df_nOut$tp, levels=rev(c("D52","D30","D11","All TP")))
figure5b_row1 <- ggplot(df_nOut, aes(x=tp, y=mean)) + 
  geom_pointrange(aes(ymin=mean-error, ymax=mean+error, fatten = 10))+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,3,0.5), lim = c(0, 3), name = "Outlier events per outlier line")+
  xlab("")+
  theme(axis.title=element_text(size=13),
        axis.text=element_text(size=12))




outAnalysis_onlyavail$trajectory <- paste0(outAnalysis_onlyavail$nOutliersTP11>0,"-",outAnalysis_onlyavail$nOutliersTP30>0,"-", outAnalysis_onlyavail$nOutliersTP52>0)

vec <- c("D11, D30, D52", "D11, D30","D11, D52","D11", "D30, D52", "D30","D52")
names(vec) <- unique(outAnalysis_onlyavail$trajectory)
outAnalysis_onlyavail$trajectory <- as.list(unname(vec[outAnalysis_onlyavail$trajectory]))
outAnalysis_onlyavail$trajectory <- sapply(outAnalysis_onlyavail$trajectory, function(x) strsplit(x,", "))

figure5b_3 <- outAnalysis_onlyavail[,c("donor_extended","trajectory")] %>%
  #distinct(donor_extended, .keep_all=TRUE) %>%
  ggplot(aes(x=trajectory)) +
  theme_bw()+
  geom_bar(fill="grey70", col="black")+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, size=3.5)+
  scale_x_upset(n_intersections=length(sapply(unique(outAnalysis_onlyavail$trajectory), function(x) x)))+
  scale_y_continuous(breaks = seq(0,40,5), lim = c(0, 35), name = "Number of outlier lines")+
  xlab("TP with outlier events")+
  theme_combmatrix(combmatrix.label.text =element_text(size=11))+
  theme(axis.title=element_text(size=13),
        axis.text.y=element_text(size=12))

figure5b_row1 <- ggarrange(figure5b_row1, figure5b_3,
                           ncol=2,nrow=1, common.legend = FALSE)


figure5b <- ggarrange(figure5b_row1, figure5b_row2,
                      ncol=1, nrow=2)

pdf(file="figures/mainFigs/figure5B.pdf")
plot(figure5b)
dev.off()


#############
### fig5C ###
#############


readSpecifiedDir <- function(tp="D11"){

  dirToSave <- paste0("analysis/outputTabs/outAnalysis/")
  #setwd(dirToSave)

  df <- paste0(dirToSave, list.files(path=dirToSave, pattern = "D11_.+.RDS")) %>%
    map(readRDS) %>%
    bind_rows()

  return(df)

}

dfAnnot <- rbind(readSpecifiedDir(tp="D11"),
                 readSpecifiedDir(tp="D30"),
                 readSpecifiedDir(tp="D52"))
dfAnnot <- dfAnnot[!(is.na(dfAnnot$pval) & is.na(dfAnnot$corrPearson)),]
dfAnnot$pAdjusted <- p.adjust(dfAnnot$pval, "BH")

## Download Cosmic database, version 90 (GRCh37)
cat( "Annotating COSMIC Tier 1 gene-level information \n")
path_to_cosmic_db <- "inputTabs/CosmicMutantExportCensus.tsv"
cosmic_db <- read.table(path_to_cosmic_db, sep="\t", header=TRUE)
cosmic_db <- cosmic_db[cosmic_db$Tier==1,]
cosmic_db <- as.character(unique(cosmic_db$Gene.name))

## Download DDG2P genes (24.1.2020)
cat( "Annotating DDD gene-level curated list \n")
path_to_ddd_genes <- "inputTabs/DDD/DDG2P_24_1_2020.tsv"
ddd_curated <- read.csv(path_to_ddd_genes, sep="\t")
ddd_dominantMOI <- as.character(unique(ddd_curated[ddd_curated$allelic.requirement=="monoallelic",]$gene.symbol))
ddd_curated <- as.character(unique(ddd_curated$gene.symbol))

dfAnnot$CosmicT1 <- dfAnnot$gene %in% cosmic_db
dfAnnot$DDD <- dfAnnot$gene %in% ddd_curated
# saveRDS(dfAnnot, file="analysis/outputTabs/outAnalysis/allgenesPerTPAnnot.RDS")
# dfAnnot <- readRDS("analysis/outputTabs/outAnalysis/allgenesPerTPAnnot.RDS")

exampleCorr <- readRDS("analysis/outputTabs/outAnalysis/exampleCorr_KMT2D.RDS")
exampleCorr <- exampleCorr[!is.na(exampleCorr$meanExp),]
exampleAnnot <- subset(dfAnnot, tp=="D11" & gene=="KMT2D")
clusters=c("FPP-1","DA")
exampleAnnot <- exampleAnnot[match(clusters, exampleAnnot$annot),]
#exampleAnnot$annotNew <- paste0(exampleAnnot$tp,": ", exampleAnnot$annot, " (",exampleAnnot$gene," gene)")
exampleCorr <- sapply(unique(exampleCorr$annot), function(z){
  
  tmp <- subset(exampleCorr, annot==z)
  tmp$zscoMeanExp <- (tmp$meanExp-mean(tmp$meanExp))/sd(tmp$meanExp)
  tmp$zscoCfrac <- (tmp$cfrac-mean(tmp$cfrac))/sd(tmp$cfrac)
  
  return(tmp)
  
}, simplify=F)

exampleCorr <- do.call("rbind", exampleCorr)
rownames(exampleCorr) <- NULL

facet_names <- c(
  'DA'="D11: DA (gene KMT2D)",
  'FPP-1'="D11: FPP-1 (gene KMT2D)"
)

figure5c_1 <- ggplot(exampleCorr, aes(x=zscoCfrac, y=zscoMeanExp))+
  geom_point(size=4, alpha=0.5)+
  facet_grid(~annot, labeller=as_labeller(facet_names))+
  theme_bw()+
  xlab("Z-score Mean Expression")+
  ylab("Z-score Cell-type fraction")+
  theme(axis.title=element_text(size=13),
        axis.text=element_text(size=12),
        strip.text = element_text(size=12))

tt <- ggplot_build(figure5c_1)

figure5c_1 <- figure5c_1 + geom_text(data = exampleAnnot, aes(x = mean(tt$layout$panel_scales_x[[1]]$range$range),
                                     y = 9,
                                     label = paste0("italic(R) == ",signif(corrPearson,3))),
            parse = TRUE, col="black",size=5)+
  geom_text(data = exampleAnnot, aes(x = mean(tt$layout$panel_scales_x[[1]]$range$range),
                                     y = 7,
                                     label = paste0("italic(pAdj) == ",signif(pAdjusted,1))),
            parse = TRUE, col="black",size=5)



geneDistrPerTpPerAnnot <- function(dfAnnot, timePoint="D11", clusters=c("FPP-1","DA")){
  
  tmp_tp <- subset(dfAnnot, tp==timePoint)
  tmp_tp$annotNum <- paste0(tmp_tp$tp,": ",tmp_tp$annot, " (n=", tmp_tp$nDonors, " donors)")
  
  annotPval=sapply(clusters, function(x){
    
    tmp_tpAnnot <- subset(tmp_tp, annot==x)
    data.frame(tp=timePoint,
               annot=x,
               annotNum=unique(tmp_tpAnnot$annotNum),
               minCorr=min(tmp_tpAnnot[tmp_tpAnnot$pAdjusted>0.05,]$corrPearson),
               maxCorr=max(tmp_tpAnnot[tmp_tpAnnot$pAdjusted>0.05,]$corrPearson))
    
  }, simplify=F)
  
  annotPval <- do.call("rbind", annotPval)
  rownames(annotPval) <- NULL

  #tmp_tp
  
  figPerTP <- ggplot(tmp_tp[tmp_tp$annot %in% clusters,], aes(x=corrPearson))+
    geom_histogram(color="black", fill="white")+
    facet_wrap(~annotNum)+
    geom_vline(data=annotPval,
               mapping=aes(xintercept=minCorr), col="purple", linetype="dashed", size=0.75)+
    geom_vline(data=annotPval,
               mapping=aes(xintercept=maxCorr), col="purple", linetype="dashed", size=0.75)+
    theme_bw()+
    ylab("Number of genes")+
    xlab("Pearson coefficient")+
    #ggtitle("Z-score correlation: Gene expression ~ Cell-type fraction")+
    theme(plot.title=element_text(hjust=0.5, face="bold"))+
    geom_text(data=subset(annotPval, annot=="DA"),
              aes(x=minCorr, y=750), col="purple", angle=90, label="pAdj<0.05", vjust=-1, size=5)+
    geom_text(data=subset(annotPval, annot=="DA"),
              aes(x=maxCorr, y=750), col="purple", angle=90, label="pAdj<0.05", vjust=2, size=5)+
    theme(axis.title=element_text(size=13),
          axis.text=element_text(size=12),
          strip.text = element_text(size=12))
  
  return(figPerTP)

}


figure5c_2 <- geneDistrPerTpPerAnnot(dfAnnot, timePoint="D11", clusters=c("FPP-1","DA"))
  
figure5c <- ggarrange(figure5c_1, figure5c_2,
                      ncol=1, nrow=2)

pdf(file="figures/mainFigs/figure5C.pdf")
plot(figure5c)
dev.off()



#############
### fig5D ###
#############
dfAnnot <- readRDS("analysis/outputTabs/outAnalysis/allgenesPerTPAnnot.RDS")

enrichTab <- function(dfAnnot, geneSet="DDD"){
  
  allTpEnrich <- sapply(unique(dfAnnot$tp), function(x){
    
    tmp <- subset(dfAnnot, tp==x)
    
    annotPerTp <- sapply(unique(tmp$annot), function(y){
      
      tmp2 <- subset(tmp, annot==y)
      mtrixEnrich <- matrix(NA, nrow=2, ncol=2)
      rownames(mtrixEnrich) <- c("signif","nonsignif")
      colnames(mtrixEnrich) <- c("inGeneSet","notGeneSet")
      matchCol <- match(geneSet, colnames(tmp2))
      mtrixEnrich["signif","inGeneSet"] <- sum(tmp2$pAdjusted<0.05 & tmp2[,matchCol])
      mtrixEnrich["signif","notGeneSet"] <- sum(tmp2$pAdjusted<0.05 & !tmp2[,matchCol])
      mtrixEnrich["nonsignif","inGeneSet"] <- sum(!tmp2$pAdjusted<0.05 & tmp2[,matchCol])
      mtrixEnrich["nonsignif","notGeneSet"] <- sum(!tmp2$pAdjusted<0.05 & !tmp2[,matchCol])

      tmp <- data.frame(tp=x,
                        annot=y,
                        #pval=chisq.test(mtrixEnrich)$p.value,
                        pval=chisq.test(mtrixEnrich, simulate.p.value = TRUE)$p.value,
                        #pval=prop.test(mtrixEnrich)$p.value,
                        nDon=unique(tmp2$nDonors),
                        geneSetTested=geneSet)

      tmp
      
    }, simplify=F)
    
    annotPerTp <- do.call("rbind", annotPerTp)
    rownames(annotPerTp) <- NULL
    annotPerTp
    
  }, simplify=F)
  
  
  allTpEnrich <- do.call("rbind", allTpEnrich)
  rownames(allTpEnrich) <- NULL

  return(allTpEnrich)
  
}

enrichRes <- rbind(enrichTab(dfAnnot, geneSet="DDD"),
                   enrichTab(dfAnnot, geneSet="CosmicT1"))
  
enrichRes$pvalAdj <- p.adjust(enrichRes$pval,"BH")
enrichRes$signif <- pvalConverter(enrichRes$pvalAdj)
enrichRes$combination <- paste0(enrichRes$tp,"-", enrichRes$geneSetTested)
enrichRes$combination <- factor(enrichRes$combination, levels=rev(paste0(rep(c("D11","D30","D52"),2),"-", c(rep("DDD",3),rep("CosmicT1",3)))))
enrichRes$combination2 <- paste0(enrichRes$tp,"-", enrichRes$annot)
enrichRes <- addMissingNA(enrichRes)

myPalette <- colorRampPalette(brewer.pal(3, "YlOrRd"))

enrichRes$signif <- factor(enrichRes$signif, levels=c("**","*","ns"))
fig5d <- ggplot(enrichRes, aes(y=combination, x=annot, fill=signif))+
  geom_tile(colour = "black")+
  theme_classic()+
  geom_text(aes(x=annot,y=combination, label=nDon))+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=12, hjust=1),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=13, face="bold", hjust=0.5),
        legend.text=element_text(size=13),
        legend.position = "top")+
  xlab("")+ylab("")+
  scale_fill_manual(name="Gene set enrichment on\n(anti)-correlated genes",
                    breaks = levels(enrichRes$signif),
                    values=rev(myPalette(3)),
                    na.value = 'grey90')


pdf(file="figures/mainFigs/figure5D.pdf", width=10, height = 5)
plot(fig5d)
dev.off()


#############
### fig5E ###
#############

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

####
####

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

#outlierLongest <- as.data.frame(outlierLong %>% pivot_longer(-c(1), names_to="outlierClass", values_to="Group"))
colnames(outlierLong)[2:5] <- c("agg_allTP", "Day11", "Day30", "Day52")


vec4 <- c("Outlier","Non-outlier")
names(vec4) <- c("TRUE","FALSE")
outlierLong$agg_allTP <- unname(vec4[as.character(outlierLong$agg_allTP)])
outlierLong$Day11 <- unname(vec4[as.character(outlierLong$Day11)])
outlierLong$Day30 <- unname(vec4[as.character(outlierLong$Day30)])
outlierLong$Day52 <- unname(vec4[as.character(outlierLong$Day52)])


##########

outlierLong$donor_extended <- gsub("/pool","/Pool",outlierLong$donor_extended)

mutTab_logReg <- readRDS("analysis/outputTabs/mutTab_logReg2.RDS")

outlierLong$donor_id <- sapply(strsplit(outlierLong$donor_extended, "/"), function(x) x[1])
mask <- elementNROWS(strsplit(outlierLong$donor_extended, "/"))==3
outlierLong$donor_id[mask] <- gsub("/Pool.+","",outlierLong$donor_extended[elementNROWS(strsplit(outlierLong$donor_extended, "/"))==3])

######################
######################

allmerged <- merge(outlierLong, ratiosTab, by.y="cline_expanded", by.x="donor_extended")

## all TP
resAll <- glm(formula= as.factor(agg_allTP)~ratio_d52_d0, family = binomial, data=allmerged)

pval_wilcox <- function(outcomeSelected="agg_allTP", growthRatio="ratio_d52_d0"){

  df <- data.frame(outlierDef=outcomeSelected,
                   growthRatioUsed=growthRatio,
                   pvalWilcox=wilcox.test(subset(allmerged, get(outcomeSelected)=="Outlier")[,growthRatio],
                                          subset(allmerged, get(outcomeSelected)=="Non-outlier")[,growthRatio])$p.value)

  return(df)
}

pvalWilcox=rbind(pval_wilcox(outcomeSelected="agg_allTP", growthRatio="ratio_d52_d0"),
                 pval_wilcox(outcomeSelected="Day11", growthRatio="ratio_d11_d0"),
                 pval_wilcox(outcomeSelected="Day30", growthRatio="ratio_d30_d0"),
                 pval_wilcox(outcomeSelected="Day52", growthRatio="ratio_d52_d0"))



pval_logReg <- function(outcomeSelected="agg_allTP", growthRatio="ratio_d52_d0"){
  
  res <- glm(formula= as.factor(get(outcomeSelected))~get(growthRatio), family = binomial, data=allmerged)
  df <- data.frame(test.name=outcomeSelected,
                   growthRatioUsed=growthRatio,
                   pval_logReg=summary(res)$coef["get(growthRatio)","Pr(>|z|)"],
                   effectSize=summary(res)$coef["get(growthRatio)", "Estimate"])
  return(df)
}

pValLofReg <- rbind(pval_logReg(outcomeSelected="agg_allTP", growthRatio="ratio_d52_d0"),
                    pval_logReg(outcomeSelected="Day11", growthRatio="ratio_d11_d0"),
                    pval_logReg(outcomeSelected="Day30", growthRatio="ratio_d30_d0"),
                    pval_logReg(outcomeSelected="Day52", growthRatio="ratio_d52_d0"))


pvalResults <- merge(pValLofReg, pvalWilcox, by.x="test.name", by.y="outlierDef")
expSup <- function(w, digits=0) {
  sprintf(paste0("%.", digits, "f*x*10^%d"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}

logRegOutliers <- ggplot(allmerged, aes(x=Day52, y=ratio_d52_d0))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.3, alpha=0.3, size=4)+
  theme_bw()+
  ylab("Proliferation rate (Day52/D0)")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text=element_text(size=14))+
  geom_text(data=subset(pvalResults,test.name=="Day52"),
            aes(label=paste0("p==",expSup(pval_logReg, digits = 2))), parse=T, x=1.5,y=9.5, size=6)+
  # geom_text(data=subset(pvalResults,test.name=="Day52"),
  #           aes(x=1.5,y=9.5, label=paste0("Pval=", formatC(pval_logReg, format = "e", digits = 2))), size=7)+
  # geom_text(data=subset(pvalResults,test.name=="Day52"),
  #           aes(x=1.5,y=9.5, label=paste0("Pval=", formatC(pval_logReg, format = "e", digits = 2))), size=7)+
  geom_text(data=subset(pvalResults,test.name=="Day52"), aes(x=1.5,y=10.25), size=5,label="Logistic Regression")

pdf(file="figures/mainFigs/figure5E.pdf", width=10, height = 5)
plot(logRegOutliers)
dev.off()

