

#########################################################
### Differentiation efficiencies (different datasets) ###
#########################################################

library(readxl)

exomes_allinfo <- readRDS("outputTabs/exomes_allinfo_nonhyper_v2.RDS")

sensory_neurons <- "inputTabs/sensory_neurons/41588_2017_5_MOESM3_ESM.xlsx"
sensory_neurons <- as.data.frame(read_excel(sensory_neurons, sheet = 2))

macrophages <- "inputTabs/macrophages/41588_2018_46_MOESM3_ESM.xlsx"
macrophages <- as.data.frame(read_excel(macrophages, sheet="S1"))

dopaminergic_neurons <- "inputTabs/dopNeurons/predDopaminergic.csv"
dopaminergic_neurons <- read.csv(dopaminergic_neurons)

endoderm <- "inputTabs/endoderm/diff_capacity/endoderm_diff.csv"
endoderm <- read.csv(endoderm)
endoderm <- endoderm[order(endoderm$diff_efficiency),]


################
### NeuroSeq ###
################

sensory_neurons_sub <- data.frame(sample_id=sensory_neurons$HipsciID,
                                  line_id=tolower(sensory_neurons$Sample),
                                  donor=sensory_neurons$donor,
                                  replicate=sensory_neurons$Differentiaton.Replicate,
                                  protocol=sensory_neurons$`Protocol for differentiatio`,
                                  status=sensory_neurons$NeuronQuality,
                                  study="sensory_neurons")

#remove NA
sensory_neurons_sub <- sensory_neurons_sub[!is.na(sensory_neurons_sub$status),]
#move NonNeuronal label to Poor (later will be converted to Failed)
sensory_neurons_sub[sensory_neurons_sub$status %in% "NonNeuronal",]$status <- rep("Poor", sum(sensory_neurons_sub$status %in% "NonNeuronal"))

#sample_id that are duplicated in the data frame
duplicated.list <- sapply(names(table(sensory_neurons_sub$sample_id)[table(sensory_neurons_sub$sample_id)>1]),
                          function(x) unique(subset(sensory_neurons_sub, sensory_neurons_sub$sample_id==x)$status), simplify=FALSE)

discordant_list <- names(duplicated.list)[elementNROWS(duplicated.list)>1]
if (length(discordant_list)){
  warning(paste0("Cell-line/s ",paste(discordant_list, collapse=",") ," with discordant differentation outputs"))
}

mask_dup <- elementNROWS(sapply(split(sensory_neurons_sub$status, sensory_neurons_sub$sample_id), function(x) unique(x)))==1

sensory_neurons <- data.frame(sample_id=names(sapply(split(sensory_neurons_sub$status, sensory_neurons_sub$sample_id)[mask_dup], function(x) unique(x))),
                              sensory_neurons=unname(sapply(split(sensory_neurons_sub$status, sensory_neurons_sub$sample_id)[mask_dup], function(x) unique(x))))

vector <- c(paste0("Failed"),
            paste0("Successful"))
names(vector) <- c("Poor","Good")
sensory_neurons$sensory_neurons <- vector[sensory_neurons$sensory_neurons]

table(sensory_neurons$sensory_neurons)
# Failed Successful 
# 5        100

exomes_allinfo$sensory_neurons <- NA
to_add <- sensory_neurons[!sensory_neurons$sample_id %in% exomes_allinfo$donor_iPSC,]$sample_id
to_add_df <- as.data.frame(matrix(NA, length(to_add), dim(exomes_allinfo)[2]))
colnames(to_add_df) <- colnames(exomes_allinfo)
to_add_df$donor_iPSC <- to_add
exomes_allinfo <- rbind(exomes_allinfo, to_add_df)
exomes_allinfo[match(sensory_neurons$sample_id, exomes_allinfo$donor_iPSC),]$sensory_neurons <- sensory_neurons$sensory_neurons

###################
### Macrophages ###
###################

macrophages_sub <- data.frame(sample_id=macrophages$genotype_id,
                              line_id=macrophages$line_id,
                              donor=macrophages$donor,
                              replicate=macrophages$replicate,
                              protocol=NA,
                              status=macrophages$status,
                              study="macrophages")

macrophages_sub <- subset(macrophages_sub, !(macrophages_sub$status=="FC_QC_fail" | macrophages_sub$status=="RNA_QC_fail" ))

stopifnot(all(elementNROWS(sapply(unique(macrophages_sub[macrophages_sub$replicate==2,]$sample_id),
                                  function(x) unique(subset(macrophages_sub, macrophages_sub$sample_id==x)$status), simplify=FALSE))==1))

#sample_id that are duplicated in the data frame
duplicated.list <- sapply(names(table(macrophages_sub$sample_id)[table(macrophages_sub$sample_id)>1]), function(x) unique(subset(macrophages_sub, macrophages_sub$sample_id==x)$status), simplify=FALSE)
discordant_list <- names(duplicated.list)[elementNROWS(duplicated.list)>1]
if (length(discordant_list)){
  warning(paste0("Cell-line ",discordant_list ," with discordant differentation outputs"))
}

#remove RNA_QC_fail and FC_QC_fail
#FC_QC_fail mean low macrophage purity
#RNA_QC_fail means degraded RNA
macrophages <- data.frame(sample_id=names(sapply(split(macrophages_sub$status, macrophages_sub$sample_id), function(x) unique(x))),
                          macrophages=unname(sapply(split(macrophages_sub$status, macrophages_sub$sample_id), function(x) unique(x))))

vector <- c(paste0("Failed"),
            paste0("Successful"))
names(vector) <- c("Diff_fail","Success")
macrophages$macrophages <- vector[macrophages$macrophages]

table(macrophages$macrophages)
# Failed Successful 
# 23         90 

exomes_allinfo$macrophages <- NA
to_add <- macrophages[!macrophages$sample_id %in% exomes_allinfo$donor_iPSC,]$sample_id
to_add_df <- as.data.frame(matrix(NA, length(to_add), dim(exomes_allinfo)[2]))
colnames(to_add_df) <- colnames(exomes_allinfo)
to_add_df$donor_iPSC <- to_add
exomes_allinfo <- rbind(exomes_allinfo, to_add_df)
exomes_allinfo[match(macrophages$sample_id, exomes_allinfo$donor_iPSC),]$macrophages <- macrophages$macrophages


###############################
### Neuroseq (Experimental) ###
###############################

summData <- readRDS("outputTabs/summData_donorCtype_Cells.RDS")
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
summDataDiff <- subset(summDataDiff, outcome2!="Discordant")

table(sapply(split(summDataDiff$outcome, summDataDiff$donor_id), function(x) unique(x)))
# Failed Successful 
# 56        150 


summDataDiff <- summDataDiff[summDataDiff$donor_id %in% exomes_allinfo$donor_iPSC,]
summDataDiff <- summDataDiff[!duplicated(summDataDiff$donor_id),]
exomes_allinfo$dopaminergic_neurons_exp <- NA
exomes_allinfo[match(summDataDiff$donor_id, exomes_allinfo$donor_iPSC),]$dopaminergic_neurons_exp <- summDataDiff$outcome



############################
### Neuroseq (Predicted) ###
############################

dopaminergic_neurons <- "inputTabs/dopNeurons/predDopaminergic.csv"
dopaminergic_neurons <- read.csv(dopaminergic_neurons)


# predicted differentiation efficiency
dopaminergic_neurons$pred_diff_efficiency <- NA
dopaminergic_neurons[dopaminergic_neurons$model_score>0.02231,]$pred_diff_efficiency <- "Successful"
dopaminergic_neurons[dopaminergic_neurons$model_score<0.02231,]$pred_diff_efficiency <- "Failed"

table(dopaminergic_neurons$pred_diff_efficiency)
# Failed Successful 
# 188        624 

exomes_allinfo$dopaminergic_neurons_pred <- NA
to_add <- dopaminergic_neurons[!dopaminergic_neurons$donor_id %in% exomes_allinfo$donor_iPSC,]$donor_id
to_add_df <- as.data.frame(matrix(NA, length(to_add), dim(exomes_allinfo)[2]))
colnames(to_add_df) <- colnames(exomes_allinfo)
to_add_df$donor_iPSC <- to_add
exomes_allinfo <- rbind(exomes_allinfo, to_add_df)
exomes_allinfo[match(dopaminergic_neurons$donor_id, exomes_allinfo$donor_iPSC),]$dopaminergic_neurons_pred <- dopaminergic_neurons$pred_diff_efficiency


####################################
### Endoderm (Cutoff - Low 20% ) ###
####################################

metadata_table <- read.table("inputTabs/hipsci.qc1_sample_info.20170927.tsv",
                             header=TRUE, sep="\t")

stopifnot(all(!is.na(sapply(endoderm$donor, function(x) match(x, gsub(".+-","", metadata_table$name))))))
endoderm$sample_id <- metadata_table[sapply(endoderm$donor, function(x) match(x, gsub(".+-","", metadata_table$name))),]$name

exomes_allinfo$endoderm <- NA
to_add <- endoderm[!endoderm$sample_id %in% exomes_allinfo$donor_iPSC,]$sample_id
to_add_df <- as.data.frame(matrix(NA, length(to_add), dim(exomes_allinfo)[2]))
colnames(to_add_df) <- colnames(exomes_allinfo)
to_add_df$donor_iPSC <- to_add
exomes_allinfo <- rbind(exomes_allinfo, to_add_df)
exomes_allinfo[match(endoderm$sample_id, exomes_allinfo$donor_iPSC),]$endoderm <- endoderm$diff_efficiency

write.table(exomes_allinfo, file="../suppTabs/suppTable1.txt",
            quote=F,
            row.names=F,
            col.names = T,
            sep="\t")


### Diffentiation cell-lines per dataset (Overlap burden + outcome) ###


sum(!is.na(exomes_allinfo$all) & !is.na(exomes_allinfo$sensory_neurons))
# [1] 85
# 
sum(!is.na(exomes_allinfo$all) & !is.na(exomes_allinfo$endoderm))
# [1] 86
# 
sum(!is.na(exomes_allinfo$all) & !is.na(exomes_allinfo$macrophages))
# [1] 102
# 
sum(!is.na(exomes_allinfo$all) & !is.na(exomes_allinfo$dopaminergic_neurons_exp))
# [1] 126
# 
sum(!is.na(exomes_allinfo$all) & !is.na(exomes_allinfo$dopaminergic_neurons_pred))
# [1] 349







