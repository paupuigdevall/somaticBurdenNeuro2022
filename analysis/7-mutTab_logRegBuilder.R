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

##############################
####  Logistic regression ####
##############################

pathToFile <- "outputTabs/iPSC/"
ptvMutBurden <- readRDS(paste0(pathToFile, "mutBurden_ptv.RDS"))
bcorPositive <- names(which(ptvMutBurden["BCOR",]>0))

mutTab_logReg <- data.frame(cline=colnames(ptvMutBurden),
                            NeuroSeq_obs=NA,
                            NeuroSeq_pred=NA,
                            sensoryNeurons_outcome=NA,
                            macrophages_outcome=NA)

##################
## NeuroSeq_obs ##
##################

summData <- readRDS("outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)
### define successful and failed lines based on our data (to include KO, basically!)
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
summDataDiff <- subset(summDataDiff, outcome2!="Discordant")

mutTab_logReg$NeuroSeq_obs <- summDataDiff[match(mutTab_logReg$cline, summDataDiff$donor_id),]$outcome



###################
## NeuroSeq_pred ##
###################

dopaminergic_neurons <- "inputTabs/dopNeurons/predDopaminergic.csv"
dopaminergic_neurons <- read.csv(dopaminergic_neurons)
mutTab_logReg$NeuroSeq_pred <- dopaminergic_neurons[match(mutTab$cline, dopaminergic_neurons$donor_id),]$model_score>0.02231
vec <- c("Failed","Successful")
names(vec) <- c("FALSE","TRUE")
mutTab_logReg$NeuroSeq_pred <- unname(vec[as.character(mutTab_logReg$NeuroSeq_pred)])


############################
## sensoryNeurons_outcome ##
############################

sensory_neurons <- "inputTabs/sensory_neurons/41588_2017_5_MOESM3_ESM.xlsx"
sensory_neurons <- as.data.frame(read_excel(sensory_neurons, sheet = 2))
sensory_neurons <- sensory_neurons[!is.na(sensory_neurons$NeuronQuality),]
#move NonNeuronal label to Poor (later will be converted to Failed)
sensory_neurons[sensory_neurons$NeuronQuality %in% "NonNeuronal",]$NeuronQuality <- rep("Poor", sum(sensory_neurons$NeuronQuality %in% "NonNeuronal"))

#sample_id that are duplicated in the data frame
duplicated.list <- sapply(names(table(sensory_neurons$HipsciID)[table(sensory_neurons$HipsciID)>1]),
                          function(x) unique(subset(sensory_neurons, sensory_neurons$HipsciID==x)$NeuronQuality), simplify=FALSE)
discordant_list <- names(duplicated.list)[elementNROWS(duplicated.list)>1]
mask_dup <- elementNROWS(sapply(split(sensory_neurons$NeuronQuality, sensory_neurons$HipsciID), function(x) unique(x)))==1
sensory_neurons_sub <- data.frame(sample_id=names(sapply(split(sensory_neurons$NeuronQuality, sensory_neurons$HipsciID)[mask_dup], function(x) unique(x))),
                                  sensory_neurons=unname(sapply(split(sensory_neurons$NeuronQuality, sensory_neurons$HipsciID)[mask_dup], function(x) unique(x))))
vector <- c("Failed","Successful")
names(vector) <- c("Poor","Good")
sensory_neurons_sub$sensory_neurons <- unname(vector[sensory_neurons_sub$sensory_neurons])

mutTab_logReg$sensoryNeurons_outcome <- sensory_neurons_sub[match(mutTab_logReg$cline, sensory_neurons$HipsciID),]$sensory_neurons

#########################
## macrophages_outcome ##
#########################

macrophages <- "inputTabs/macrophages/41588_2018_46_MOESM3_ESM.xlsx"
macrophages <- as.data.frame(read_excel(macrophages, sheet="S1"))
macrophages <- subset(macrophages, status!="RNA_QC_fail")

duplicated.list <- sapply(names(table(macrophages$genotype_id)[table(macrophages$status)>1]),
                          function(x) unique(subset(macrophages, macrophages$genotype_id==x)$status), simplify=FALSE)
discordant_list <- duplicated.list[elementNROWS(duplicated.list)>1]
to_discard <- names(which(!sapply(discordant_list, function(x) all(grepl("fail",x)))))
macrophages <- macrophages[!macrophages$genotype_id %in% to_discard,]

duplicated.list <- sapply(names(table(macrophages$genotype_id)[table(macrophages$status)>1]),
                          function(x) unique(subset(macrophages, macrophages$genotype_id==x)$status), simplify=FALSE)
duplicated.list[elementNROWS(duplicated.list)>1][[1]] <- duplicated.list[elementNROWS(duplicated.list)>1][[1]][2]
stopifnot(all(elementNROWS(duplicated.list)==1))

macrophages_sub <- data.frame(sample_id=names(unlist(duplicated.list)),
                              macrophages=unname(unlist(duplicated.list)))

vector <- c("Failed","Failed","Successful")
names(vector) <- c("Diff_fail","FC_QC_fail","Success")
macrophages_sub$macrophages <- unname(vector[macrophages_sub$macrophages])
mutTab_logReg$macrophages_outcome <- macrophages_sub[match(mutTab_logReg$cline, macrophages_sub$sample_id),]$macrophages

mutTab_logReg$BCORpos <- 
  unname(ptvMutBurden[match("BCOR",
                            rownames(ptvMutBurden)),][match(mutTab_logReg$cline, names(ptvMutBurden[match("BCOR", rownames(ptvMutBurden)),]))]>0)


library(readxl)
metadata_table <- read.table("inputTabs/hipsci.qc1_sample_info.20170927.tsv",
                             header=TRUE, sep="\t")
mutTab_logReg$phenotype <- metadata_table[match(mutTab_logReg$cline, metadata_table$name),]$disease
saveRDS(mutTab_logReg, file="outputTabs/mutTab_logReg2.RDS")





