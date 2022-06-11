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

pdf(file=paste0("figures/suppFigs/suppfig1C.pdf"))
plot(bcorPlot)
dev.off()

