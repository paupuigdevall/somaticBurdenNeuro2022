
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

source("functionsToImport.R")

##############
### Fig 2A ###
##############

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

plot2aFunction <- function(outcomes, var="all", mainTitle="All variants (SNVs+dinuc+indels+CNVs)"){
  
  figure2aTab_all <- subset(outcomes, varCategory==var)
  newvec <- c("Sensory N","Macrophages","DopaN (Obs)","DopaN (Pred)")
  names(newvec) <- unique(outcomes$differentiation)
  figure2aTab_all$differentiation2 <- unname(newvec[figure2aTab_all$differentiation])
  
  my_comparisons <- list(c("Failed","Successful"))
  
  plot2a_1 <- ggplot(figure2aTab_all, aes(x=differentiation2, y=mutBurden, fill=outcome))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha=0.2, size=0.25)+
    theme_bw()+
    xlab("")+
    ylab("")+
    scale_fill_manual(name="Differentiation outcome",values=c('#E69F00', '#56B4E9'))+
    ggtitle(mainTitle)+
    theme(axis.text.x=element_text(angle=60, size=13, vjust=0.5, colour="black"),
          axis.text.y=element_text(size=15, colour="black"),
          axis.title = element_text(size=16),
          plot.title=element_text(hjust=0.5, face="bold", size=15),
          legend.title=element_text(face="bold", size=14),
          legend.text=element_text(size=13),
          legend.position="bottom")+
    stat_compare_means(label = "p.format", size=4.25, method = "wilcox.test")
  
  plot(plot2a_1)
  
  return(plot2a_1)
  
}

p2a_1 <- plot2aFunction(outcomes, var="all", mainTitle="All variants")
p2a_2 <- plot2aFunction(outcomes, var="deleterious", mainTitle="Deleterious variants")

figure2a <- ggarrange(p2a_1, p2a_2,
                      ncol=2, nrow=1, common.legend = TRUE, legend="top")

#figure2a <- annotate_figure(figure2a, left=text_grob("Somatic mutational burden per cell-line", rot=90, size=16, face="bold"))

pdf(file=paste0("figures/mainFigs/figure2A.pdf"))
plot(figure2a)
dev.off()

########################
########################

##############
### Fig 2B ###
##############

colnames_of_interest <- c("donor_iPSC","endoderm","all","synonymous","other","deleterious","missPatho","ptv")
newendoderm <- diffmetrics[,match(c(colnames_of_interest), colnames(diffmetrics))]
matchCol <- match("endoderm", colnames(newendoderm))
colnames(newendoderm)[matchCol] <- "outcome"
newendoderm <- as.data.frame(newendoderm %>% pivot_longer(-c("donor_iPSC","outcome"), names_to="varCategory", values_to="mutBurden"))
newendoderm <- newendoderm[!is.na(newendoderm$mutBurden) & !is.na(newendoderm$outcome),]

newvec <- c("All variants","Synonymous","Non-coding","Deleterious","Missense Pathogenic","Loss-of-function")
names(newvec) <- unique(newendoderm$varCategory)
newendoderm$varCategory <- unname(newvec[newendoderm$varCategory])
newendoderm$varCategory <- factor(newendoderm$varCategory, levels=unique(newendoderm$varCategory))

## pearson correlation between diff.efficiencies and somatic burden
corr <- sapply(unique(newendoderm$varCategory), function(x){
  varClass <- subset(newendoderm, varCategory==x)
  reg <- cor.test(varClass$outcome, varClass$mutBurden, 
                  method = "pearson")
  obj <- data.frame(varCategory=x,
                    coef=reg$estimate,
                    pval=reg$p.value)
  obj
  
}, simplify=F)
corr <- do.call("rbind", corr)
rownames(corr) <- NULL
corr$pvalAdj <- p.adjust(corr$pval, "BH")

plot2bFunction <- function(newendoderm, var="All variants", mainTitle="All variants (SNVs+dinuc+indels+CNVs)"){
  
  figure2bTab_all <- subset(newendoderm, varCategory==var)
  corr_all <- subset(corr, varCategory==var)
  endodermPlot <- ggplot(figure2bTab_all, aes(x=outcome, y=mutBurden))+
    geom_point(size=3, alpha=0.5, col="black")+
    theme_bw()+
    scale_x_continuous(breaks=seq(0,1,0.2),
                       labels=seq(0,1,0.2))+
    coord_cartesian(xlim=c(0,1))
  
  ycurrent=ggplot_build(endodermPlot)$layout$panel_scales_y[[1]]$range$range[2]
  y1=ycurrent-0.05*ycurrent
  y2=ycurrent-0.*ycurrent
  
  endodermPlot <- endodermPlot+
    xlab("")+
    ylab("")+
    geom_text(data=corr_all,
              aes(x=0.2, y=y1, label=paste0("R= ",signif(coef,2))), col="black",size=5, inherit.aes = T)+
    geom_text(data=corr_all,
              aes(x=0.2, y=y2, label=paste0("pAdj= ",signif(pvalAdj,2))), col="black",size=5, inherit.aes = T)+
    ggtitle(mainTitle)+
    theme(plot.title=element_text(hjust=0.5, face="bold", size=15),
          axis.text=element_text(size=15, colour="black"))
  
  return(endodermPlot)
  
}

p2b_1 <- plot2bFunction(newendoderm, var="All variants", mainTitle="All variants")
p2b_2 <- plot2bFunction(newendoderm, var="Deleterious", mainTitle="Deleterious variants")


figure2b <- ggarrange(p2b_1, p2b_2,
                      ncol=2, nrow=1, common.legend = TRUE, legend="top")

# figure2b <- annotate_figure(figure2b,
#                             left=text_grob("Somatic mutational burden per cell-line", rot=90, size=14, face="bold"),
#                             bottom=text_grob("Endoderm differentiation efficiency", size=14, face="bold"))
# 

pdf(file=paste0("figures/mainFigs/figure2B.pdf"))
plot(figure2b)
dev.off()


##############
### Fig 2C ###
##############

geneUniverse <- readRDS("analysis/outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
mutTab_logReg <- readRDS("analysis/outputTabs/mutTab_logReg2.RDS")
pathToFile <- "analysis/outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
synMutBurden <- readRDS(paste0(pathToFile, "mutBurden_synonymous.RDS"))
delMutBurden <- readRDS(paste0(pathToFile, "mutBurden_del.RDS"))



corrDataSet_logReg <- function(varSet="ptvMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_obs"){
  
  regSummary <- sapply(rownames(ptvMutBurden), function(x){
    #print(x)
    mutTable$gene <- unname(get(varSet)[x,][match(mutTable$cline, names(get(varSet)[x,]))])
    
    pvalWilcox <- wilcox.test(subset(mutTable, get(dataSet)=="Failed")$gene, subset(mutTable, get(dataSet)=="Successful")$gene)$p.value
    pvalWilcox[is.nan(pvalWilcox)] <- NA
    
    df <- data.frame(gene=x,
                     pval=pvalWilcox,
                     nFailed=length(subset(mutTable, get(dataSet)=="Failed")$gene),
                     nSucc=length(subset(mutTable, get(dataSet)=="Successful")$gene),
                     mutFailed=sum(subset(mutTable, get(dataSet)=="Failed")$gene),
                     mutSucc=sum(subset(mutTable, get(dataSet)=="Successful")$gene),
                     genomicLength=geneUniverse[match(x, geneUniverse$symbol),]$genomic_length,
                     totBurden=sum(subset(mutTable, get(dataSet)=="Failed" | get(dataSet)=="Successful")$gene),
                     mutState=NA,
                     norm_mutFailed=NA,
                     norm_mutSucc=NA,
                     log2FC=NA)
    
    df$norm_mutFailed <- df$mutFailed*1000/(df$genomicLength*df$nFailed)
    df$norm_mutSucc <- df$mutSucc*1000/(df$genomicLength*df$nSucc)
    
    df
    
  }, simplify=F)
  
  regSummary <- do.call("rbind", regSummary)
  rownames(regSummary) <- NULL
  regSummary$adjPval <- NA
  mask <- !is.na(regSummary$pval)
  regSummary$adjPval[mask] <- p.adjust(regSummary$pval[mask],"BH")
  regSummary$signif <- pvalConverter(regSummary$pval)
  regSummary$signifAdj <- pvalConverter(regSummary$adjPval)
  
  pseudocount <- c(regSummary$norm_mutFailed, regSummary$norm_mutSucc)
  pseudocount <- pseudocount[pseudocount!=0]
  pseudocount <- min(pseudocount)
  regSummary$log2FC <- log2((regSummary$norm_mutFailed+pseudocount)/(regSummary$norm_mutSucc+pseudocount))
  regSummary$dataset <- dataSet
  regSummary$varCategory <- varSet
  
  ##FC>2.5, FC<1/2.5
  up <- log(2.5)/log(2)
  down <- -log(2.5)/log(2)
  regSummary[regSummary$log2FC>up,]$mutState <- "up"
  regSummary[regSummary$log2FC<down,]$mutState <- "down"
  regSummary[ (regSummary$log2FC<=up & regSummary$log2FC>=down),]$mutState <- "eq"
  
  return(regSummary)
  
}


allRes <- rbind(corrDataSet_logReg(varSet="ptvMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_obs"),
                corrDataSet_logReg(varSet="ptvMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_pred"),
                corrDataSet_logReg(varSet="ptvMutBurden", mutTable=mutTab_logReg, dataSet="macrophages_outcome"),
                corrDataSet_logReg(varSet="synMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_obs"),
                corrDataSet_logReg(varSet="synMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_pred"),
                corrDataSet_logReg(varSet="synMutBurden", mutTable=mutTab_logReg, dataSet="macrophages_outcome"),
                corrDataSet_logReg(varSet="delMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_obs"),
                corrDataSet_logReg(varSet="delMutBurden", mutTable=mutTab_logReg, dataSet="NeuroSeq_pred"),
                corrDataSet_logReg(varSet="delMutBurden", mutTable=mutTab_logReg, dataSet="macrophages_outcome"))
saveRDS(allRes, file="analysis/outputTabs/mutTab/diffMutRatio.RDS")

#allRes <- readRDS("analysis/outputTabs/mutTab/diffMutRatio.RDS")
up <- log(2.5)/log(2)
down <- -log(2.5)/log(2)


vecMutState <- c("None","Successful","Failed")
names(vecMutState) <- unique(allRes$mutState)
allRes$mutState <- unname(vecMutState[allRes$mutState])
allRes$mutState <- factor(allRes$mutState, levels=c("Successful","None","Failed"))

vecDataset <- c("DopaminergicN (Obs)", "DopaminergicN (Pred)", "Macrophages")
names(vecDataset) <- unique(allRes$dataset)
allRes$dataset <- unname(vecDataset[allRes$dataset])
allRes$dataset <- factor(allRes$dataset, levels=c("DopaminergicN (Obs)", "DopaminergicN (Pred)", "Macrophages"))

allRes$varCategory <- as.character(allRes$varCategory)
vecVarCategory <- c("LoF variant","Synonymous","Deleterious")
names(vecVarCategory) <- unique(allRes$varCategory)
allRes$varCategory <- unname(vecVarCategory[allRes$varCategory])
allRes$varCategory <- factor(allRes$varCategory, levels=c("LoF variant", "Deleterious", "Synonymous"))

my_breaksy <- c(0.1,1,10,50)

allRes$diffExpressed <- "grey50"
allRes$diffExpressed[!is.na(allRes$signifAdj) & allRes$signifAdj<0.05 & allRes$log2FC>up] <- "red"
allRes$diffExpressed[!is.na(allRes$signifAdj) & allRes$signifAdj<0.05 & allRes$log2FC<down] <- "blue"


volcanoPlot <- ggplot(allRes, aes(x=log2FC, y=-log10(adjPval), col=diffExpressed, size=diffExpressed))+
  geom_point(alpha=0.5)+
  theme_bw()+
  geom_text_repel(data=subset(allRes, (signifAdj=="***" & (log2FC>up | log2FC<down))), 
                  aes(label=gene), 
                  col="black",
                  seed=123,
                  box.padding=1, max.overlaps = 999, size=4)+
  geom_vline(xintercept=c(down, up), col="grey60", linetype="dashed", alpha=0.7) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype="dashed", alpha=0.7)+
  facet_grid(varCategory~dataset)+
  theme(legend.position = "top")+
  scale_color_manual(name="Differentially mutated in:",
                     values=c("blue","grey60","red"), 
                     labels=c("Successful","None","Failed"))+
  scale_size_manual(name="", values=c(2,0.25,2),
                    labels=c("Successful","None","Failed"))+
  guides(size=F)+
  xlab("log2FC")+
  ylab("-log10(Adj.Pval)")+
  scale_x_continuous(breaks=seq(-15,15,5), labels=seq(-15,15,5))+
  scale_y_continuous(trans = log2_trans(), breaks = my_breaksy)+
  theme(strip.text = element_text(size = 11),
        axis.text.x=element_text(angle=90, size=11, colour="black", hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=11, colour="black"),
        legend.text = element_text(size=12, colour="black"),
        legend.key.size = unit(1, 'cm'))


pdf(file=paste0("figures/mainFigs/figure2C.pdf"), width=6, height = 6)
plot(volcanoPlot)
dev.off()


##############
### Fig 2D ###
##############

runGO <- list(c("observed_del","failedGenes", "Failed_Observed_Del"),
              c("predicted_del","failedGenes", "Failed_Predicted_Del"),
              c("observed_del","successfulGenes", "Succ_Observed_Del"),
              c("predicted_del","successfulGenes", "Succ_Predicted_Del")) 


GOresults <- read.table("suppTabs/suppTable3.txt",
                         sep="\t", header=T)

GOresults <- sapply(unique(GOresults$label), function(x) subset(GOresults, label==x), simplify=F)

highlightTopOrAndNeuroGO_adapted <- function(enrichment_failed_all){
  
  enrichment_failed_all$GeneSyms <- NULL
  enrichment_failed_all$position <- c(1:dim(enrichment_failed_all)[1])
  enrichment_failed_all$highlight <- FALSE
  #enrichment_failed_all$highlight[1:5] <- TRUE
  
  vecMatch <- c("axon","neuron","glial","brain",
                "hindbrain","forebrain","midbrain","synapse","chromatin",
                "cerebellum", "neural","cortex","neurogenesis",
                "axonogenesis","nervous","hippocampus","neurotransmitter",
                "dopaminergic", "axenome", "action potential","synaptic")
  
  mask <- rowSums(sapply(vecMatch, function(x){
    grepl(x,enrichment_failed_all$Term)
  }))>0
  if (sum(mask)){
    enrichment_failed_all[mask,]$highlight <- TRUE
    enrichment_failed_all <- subset(enrichment_failed_all, highlight==T)
    return(enrichment_failed_all)
  } else {
    return(NA)
  }
}

GOresultsHigh <- sapply(GOresults, function(y){
  highlightTopOrAndNeuroGO_adapted(y)
}, simplify=F)

GOresultsHigh <- GOresultsHigh[!sapply(GOresultsHigh, function(x) all(is.na(x)))]
GOresultsHigh <- do.call("rbind", GOresultsHigh)
rownames(GOresultsHigh) <- NULL
GOresultsHigh$label <- gsub("_"," ", GOresultsHigh$label)
GOresultsHigh$label <- gsub(" Del","", GOresultsHigh$label)

#GOresultsHigh$highlight[is.na(GOresultsHigh$highlight)] <- F
labels=gsub(" Del","",gsub("_"," ",sapply(runGO, function(x) x[3])))


fillNA_fig2 <- function(enrichmentGOexclusive, highlight=F, labels=labels){
  alltabs <- sapply(unique(enrichmentGOexclusive$Term), function(x){
    tabs <- sapply(labels, function(y){
      dimensions <- dim(subset(enrichmentGOexclusive, Term==x & label==y))[1]
      if (dimensions==0){
        df <- data.frame(GOBPID=NA,
                         Pvalue="ns",
                         OddsRatio=NA,
                         ExpCount=NA,
                         Count=NA,
                         Size=NA,
                         Term=x,
                         label=y,
                         position=NA,
                         highlight=F)
      } else {
        df <- NA
      }
      return(df)
    }, simplify=F)
    
    tabs <- tabs[!sapply(tabs, function(x) all(is.na(x)))]
    tabs <- do.call("rbind", tabs)
    rownames(tabs) <- NULL
    tabs
  }, simplify=F)
  
  alltabs <- do.call("rbind", alltabs)
  rownames(alltabs) <- NULL
  return(alltabs)
}


enrichmentGOexclusive <- rbind(GOresultsHigh,
                               fillNA_fig2(GOresultsHigh, highlight=F, labels=labels))

maxOdds <- ceiling(max(enrichmentGOexclusive$OddsRatio, na.rm=T))
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))

mask_long <- nchar(enrichmentGOexclusive$Term)>30
enrichmentGOexclusive$short <- enrichmentGOexclusive$Term
enrichmentGOexclusive$short[mask_long] <- sapply(enrichmentGOexclusive[mask_long,]$Term, function(x){
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
})

enrichmentGOexclusive$label <- factor(enrichmentGOexclusive$label,
                                      levels=unique(enrichmentGOexclusive$label))

ipscGOEnrichement <- ggplot(enrichmentGOexclusive,
                            aes(y=short, x=label, fill=OddsRatio))+
  geom_tile(colour = "black")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1, size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        plot.title=element_text(hjust=0.5, face="bold", size=14),
        axis.title=element_text(size=12, face="bold"),
        legend.position="top", 
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  scale_fill_gradientn(name="OddsRatio",
                       colours=myPalette(100),
                       limits=c(0,maxOdds),
                       labels=seq(0,maxOdds,0.5),
                       breaks=seq(0,maxOdds,0.5),
                       na.value = 'grey90',
                       guide = guide_colourbar(barwidth = 12))+
  geom_text(aes(label=position), size=5)+
  xlab("")+ylab("")

pdf(paste0("figures/mainFigs/figure2D.pdf"))
plot(ipscGOEnrichement)
dev.off()









