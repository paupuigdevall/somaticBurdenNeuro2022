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
library(RColorBrewer)

source("functionsToImport.R")

###############
### SFig 5A ###
###############

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
diffmetrics <- read.table("suppTabs/suppTable1.txt",
                          header=T)
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
outlier_analysis$allMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$all
outlier_analysis$delMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$deleterious


tpOutlierSet <- function(summData, timePoint="D11", colname="outlier_tp11_z2"){
  outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE)$donor_extended)
  colIndex <- match(colname, colnames(outlier_analysis))
  outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
  return(outlier_analysis)
}

outlier_analysis <- tpOutlierSet(summData, timePoint="D11", colname="outlier_tp11_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D30", colname="outlier_tp30_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D52", colname="outlier_tp52_z2")

tpOutlierSetPercentile <- function(summData, timePoint="D11", colname="outlier_tp11_up5"){
  
  sense <- gsub("[0-9]","",sapply(strsplit(colname,"_"), function(x) x[3]))
  
  if (sense=="up"){
    
    outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE & top5==TRUE)$donor_extended)
    colIndex <- match(colname, colnames(outlier_analysis))
    outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
    
  } else if (sense=="down"){
    
    outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE & bot5==TRUE)$donor_extended)
    colIndex <- match(colname, colnames(outlier_analysis))
    outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
    
  }
  
  return(outlier_analysis)
}


outlierLong <- outlier_analysis[,c("donor_extended","outlier_alltp_z2","outlier_tp11_z2","outlier_tp30_z2","outlier_tp52_z2","allMutBurden","delMutBurden")]
outlierLong <- as.data.frame(outlierLong %>% pivot_longer(-c(1:5), names_to="MutClass", values_to="Burden"))
outlierLong <- subset(outlierLong, !is.na(Burden))

facet_labeller = labeller(
  MutClass = c("allMutBurden" = "All variants",
               "delMutBurden" = "Deleterious"))
outlierLong$MutClass <- factor(outlierLong$MutClass, levels=unique(outlierLong$MutClass))
outlierLongest <- as.data.frame(outlierLong %>% pivot_longer(-c(1,6,7), names_to="outlierClass", values_to="Group"))

vec <- c("Agg All TP", "Day 11", "Day 30", "Day52")
names(vec) <-  unique(outlierLongest$outlierClass)
outlierLongest$outlierClass <- unname(vec[outlierLongest$outlierClass])

vec4 <- c("Outlier","Non-outlier")
names(vec4) <- c("TRUE","FALSE")
outlierLongest$Group <- unname(vec4[as.character(outlierLongest$Group)])
my_comparisons <- list(c("Outlier","Non-outlier"))

myplot <- ggplot(outlierLongest, aes(x=Group, y=Burden))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.25, alpha=0.25, aes(col=Group), size=3)+
  facet_grid(MutClass~outlierClass, scales="free_y", labeller = facet_labeller)+
  xlab("")+ylab("")+
  theme_bw()+
  theme(strip.text =element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20, angle=90,hjust=1, vjust=0.5),
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

myplot <- myplot + geom_point(data=ff, x=NA)

q <- ggplot_build(myplot)

mask_panel <- q$data[[3]]$PANEL==1

q$data[[3]]$y[mask_panel] <- q$data[[3]]$y[mask_panel]-40
q$data[[3]]$yend[mask_panel] <- q$data[[3]]$yend[mask_panel]-40

q$data[[3]]$y[!mask_panel] <- q$data[[3]]$y[!mask_panel]-20
q$data[[3]]$yend[!mask_panel] <- q$data[[3]]$yend[!mask_panel]-20


pv <- generics::tidy(with(outlierLongest, pairwise.wilcox.test(Burden, interaction(Group,outlierClass,MutClass), p.adjust.method = "BH")))

pv_final <- as.data.frame(pv %>% 
                            separate(group1, c("g1_1", "g1_2", "g1_3"), sep="\\.") %>% 
                            separate(group2, c("g2_1", "g2_2", "g2_3"), sep="\\.") %>% 
                            filter(g1_1 != g2_1) %>% filter(g1_2 == g2_2) %>% filter(g1_3 == g2_3) %>%
                            mutate(p=ifelse(p.value > 0.05, "", ifelse(p.value > 0.01,"*", "**"))))


# each pvalue is repeated three times in this dataset.  s
q$data[[3]]$annotation <- rep(pv_final$p, each=3)
# remove non significants
q$data[[3]] <- q$data[[3]][q$data[[3]]$annotation != "",]
# and the final plot
figure <- ggplot_gtable(q)

pdf(file=paste0("figures/suppFigs/suppfig5A.pdf"))
plot(figure)
dev.off()


###############
### SFig 5B ###
###############


## column 1
dfAnnot <- readRDS("analysis/outputTabs/outAnalysis/allgenesPerTPAnnot.RDS")
exampleCorr <- readRDS("analysis/outputTabs/outAnalysis/exampleCorr_BCOR_D11_FPP-3.RDS")
exampleCorr <- exampleCorr[!is.na(exampleCorr$meanExp),]
exampleAnnot <- subset(dfAnnot, tp=="D11" & gene=="BCOR")
clusters=c("FPP-3")
exampleAnnot <- exampleAnnot[match(clusters, exampleAnnot$annot),]

exampleCorr <- sapply(unique(exampleCorr$annot), function(z){
  
  tmp <- subset(exampleCorr, annot==z)
  tmp$zscoMeanExp <- (tmp$meanExp-mean(tmp$meanExp))/sd(tmp$meanExp)
  tmp$zscoCfrac <- (tmp$cfrac-mean(tmp$cfrac))/sd(tmp$cfrac)
  
  return(tmp)
  
}, simplify=F)

exampleCorr <- do.call("rbind", exampleCorr)
rownames(exampleCorr) <- NULL

## column 2
dfAnnot <- readRDS("analysis/outputTabs/outAnalysis/allgenesPerTPAnnot.RDS")
exampleCorr2 <- readRDS("analysis/outputTabs/outAnalysis/exampleCorr_BCOR_D52_DA.RDS")
exampleCorr2 <- exampleCorr2[!is.na(exampleCorr2$meanExp),]
exampleAnnot2 <- subset(dfAnnot, tp=="D52" & gene=="BCOR")
clusters=c("DA")
exampleAnnot2 <- exampleAnnot2[match(clusters, exampleAnnot2$annot),]
#exampleAnnot$annotNew <- paste0(exampleAnnot$tp,": ", exampleAnnot$annot, " (",exampleAnnot$gene," gene)")
exampleCorr2 <- sapply(unique(exampleCorr2$annot), function(z){
  
  tmp <- subset(exampleCorr2, annot==z)
  tmp$zscoMeanExp <- (tmp$meanExp-mean(tmp$meanExp))/sd(tmp$meanExp)
  tmp$zscoCfrac <- (tmp$cfrac-mean(tmp$cfrac))/sd(tmp$cfrac)
  
  return(tmp)
  
}, simplify=F)

exampleCorr2 <- do.call("rbind", exampleCorr2)
rownames(exampleCorr2) <- NULL

exampleCorr <- rbind(exampleCorr, exampleCorr2)
exampleAnnotNew <- rbind(exampleAnnot, exampleAnnot2)

vec <- c(
  'FPP-3'="D11: FPP-3 (gene BCOR)",
  'DA'="D52: DA (gene BCOR)"
)

exampleCorr$annot <- unname(vec[exampleCorr$annot])
exampleCorr$annot <- factor(exampleCorr$annot, levels=c("D11: FPP-3 (gene BCOR)",
                                                        "D52: DA (gene BCOR)"))
exampleAnnotNew$annot <- unname(vec[exampleAnnotNew$annot])
exampleAnnotNew$annot <- factor(exampleAnnotNew$annot, levels=c("D11: FPP-3 (gene BCOR)",
                                                                "D52: DA (gene BCOR)"))

figure5c_1 <- ggplot(exampleCorr, aes(x=zscoCfrac, y=zscoMeanExp))+
  geom_point(size=4, alpha=0.5)+
  facet_grid(~annot)+
  theme_bw()+
  xlab("Z-score Mean Expression")+
  ylab("Z-score Cell-type fraction")+
  theme(axis.title=element_text(size=13),
        axis.text=element_text(size=12),
        strip.text = element_text(size=12))

tt <- ggplot_build(figure5c_1)

figure5c_1 <- figure5c_1 + geom_text(data = exampleAnnotNew, aes(x = mean(tt$layout$panel_scales_x[[1]]$range$range),
                                                                 y = 9,
                                                                 label = paste0("italic(R) == ",signif(corrPearson,3))),
                                     parse = TRUE, col="black",size=5)+
  geom_text(data = exampleAnnotNew, aes(x = mean(tt$layout$panel_scales_x[[1]]$range$range),
                                        y = 7,
                                        label = paste0("italic(pAdj) == ",signif(pAdjusted,3))),
            parse = TRUE, col="black",size=5)


geneDistrPerTpPerAnnot <- function(dfAnnot, timePoint=c("D11","D52"),  combinations=list(c("D11","FPP-3"),c("D52","DA"))){
  
  tmp_tp <- subset(dfAnnot, tp %in% timePoint)
  tmp_tp$annotNum <- paste0(tmp_tp$tp,": ",tmp_tp$annot, " (n=", tmp_tp$nDonors, " donors)")
  
  
  annotPval=sapply(combinations, function(x){
    
    tmp_tpAnnot <- subset(tmp_tp, annot==x[2] & tp==x[1])
    data.frame(tp=x[1],
               annot=x[2],
               annotNum=unique(tmp_tpAnnot$annotNum),
               minCorr=min(tmp_tpAnnot[tmp_tpAnnot$pAdjusted>0.05,]$corrPearson),
               maxCorr=max(tmp_tpAnnot[tmp_tpAnnot$pAdjusted>0.05,]$corrPearson))
    
  }, simplify=F)
  
  annotPval <- do.call("rbind", annotPval)
  rownames(annotPval) <- NULL
  
  #tmp_tp
  tmp_tpNew <- sapply(combinations, function(x){
    
    tmp_tpAnnot <- subset(tmp_tp, annot==x[2] & tp==x[1])
    tmp_tpAnnot
  }, simplify=F)
  
  tmp_tpNew <- do.call("rbind", tmp_tpNew)
  rownames(tmp_tpNew) <- NULL
  
  clusters <- sapply(combinations, function(y) y[2])
  
  exampleAnnotNew2 <- exampleAnnotNew
  exampleAnnotNew2$annot <- sapply(strsplit(sapply(strsplit(as.character(exampleAnnotNew2$annot), ": "), function(x) x[2]), " "), function(x) x[1])
  exampleAnnotNew2 <-  merge(exampleAnnotNew2, annotPval)
  
  figPerTP <- ggplot(tmp_tpNew[tmp_tpNew$annot %in% clusters,], aes(x=corrPearson))+
    geom_histogram(color="black", fill="white")+
    facet_wrap(~annotNum)+
    geom_vline(data=annotPval,
               mapping=aes(xintercept=minCorr), col="purple", linetype="dashed", size=0.75)+
    geom_vline(data=annotPval,
               mapping=aes(xintercept=maxCorr), col="purple", linetype="dashed", size=0.75)+
    theme_bw()+
    ylab("Number of genes")+
    xlab("Pearson coefficient")+
    theme(plot.title=element_text(hjust=0.5, face="bold"))+
    geom_text(data=subset(annotPval, annot==annot),
              aes(x=minCorr, y=2000), col="purple", angle=90, label="pAdj<0.05", vjust=-1, size=5)+
    geom_text(data=subset(annotPval, annot==annot),
              aes(x=maxCorr, y=2000), col="purple", angle=90, label="pAdj<0.05", vjust=2, size=5)+
    geom_point(data=subset(exampleAnnotNew2, annotNum==annotNum), aes(x=corrPearson, y=1000), col="red", size =4)+
    geom_text(data=subset(exampleAnnotNew2, annotNum==annotNum), aes(x=corrPearson, y=1000, label=gene, vjust=-1, hjust=0.15))+
    theme(axis.title=element_text(size=13),
          axis.text=element_text(size=12),
          strip.text = element_text(size=12))
  
  return(figPerTP)
  
}


figure5c_2 <- geneDistrPerTpPerAnnot(dfAnnot, timePoint=c("D11","D52"),  combinations=list(c("D11","FPP-3"),c("D52","DA")))

figure5b <- ggarrange(figure5c_1, figure5c_2,
                      ncol=1, nrow=2)

pdf(file=paste0("figures/suppFigs/suppfig5B.pdf"))
plot(figure5b)
dev.off()

















