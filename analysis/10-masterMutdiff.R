library(VariantAnnotation)
library(dplyr)
library(tidyverse)

## each replicate pair is processed
args = commandArgs(trailingOnly=TRUE)
failedVcf <- args[1]
succVcf <- args[2]
ipsFailed <- gsub("\\..+","",sapply(strsplit(failedVcf,"\\/"), function(x) x[length(x)]))
ipsSucc <- gsub("\\..+","",sapply(strsplit(succVcf,"\\/"), function(x) x[length(x)]))

geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
vcfFileFailed <- readVcf(failedVcf, genome="hg19")
vcfFileSucc <- readVcf(succVcf, genome="hg19")

## annotation

## We look for mutations exclusively mutated in certain genes in either failed or successful cell-lines
vcfFileFailed_exclusive <- vcfFileFailed[!names(rowRanges(vcfFileFailed)) %in% names(rowRanges(vcfFileSucc)),]
vcfFileSucc_exclusive <- vcfFileSucc[!names(rowRanges(vcfFileSucc)) %in% names(rowRanges(vcfFileFailed)),]

annotateExclusiveMutations <- function(vcfFile, ipsFailed){
  csqList <- info(vcfFile)$CSQ
  names(csqList) <- rownames(info(vcfFile))
  
  vecMask <- sapply(csqList, function(x){
    sapply(strsplit(x,"\\|"), function(y) y[4] %in%  geneUniverse$symbol)
  },simplify=T)
  
  vcfFile <- vcfFile[vecMask]
  
  csqList <- info(vcfFile)$CSQ
  names(csqList) <- rownames(info(vcfFile))
  
  vep_severity <- read.csv("inputTabs/vep_severity.csv", header=TRUE)
  consequence <- as.character(1:length(vep_severity$SO_term))
  names(consequence) <- gsub("_variant","",vep_severity$SO_term)
  undef <- as.character(length(consequence)+1) 
  names(undef) <- "."
  consequence_bcftools <- c(consequence, undef)

  annot <- sapply(csqList, function(x){
    sapply(strsplit(x,"\\|"), function(y) y[2])
  },simplify=T)
  annot <- sapply(annot, function(x) strsplit(x,"\\&"))
  annot <- sapply(annot, function(x) gsub("_variant","",x))
  maskTest <- sapply(annot, function(x) x %in% names(consequence_bcftools))
  #filter by geneUniverse
  
  stopifnot(all(unlist(maskTest)))
  annot <- sapply(annot, function(x) names(which.min(consequence_bcftools[x])), simplify=T)
  
  
  tmpConsequence <- data.frame(var=names(annot),
                               varCategoryRaw=unname(annot))
  tmpConsequence$REF <- sapply(strsplit(names(annot), "_"), function(x) gsub("/.+","",x[2]))
  tmpConsequence$ALT <- sapply(strsplit(names(annot), "_"), function(x) gsub(".+/","",x[2]))
  
  tmpConsequence$class <- NA
  mask_snvs <- nchar(tmpConsequence$REF)==1 & nchar(tmpConsequence$ALT)==1
  tmpConsequence$class[mask_snvs] <- "snvs_dinuc"
  mask_dinuc <- nchar(tmpConsequence$REF)==2 & nchar(tmpConsequence$ALT)==2
  tmpConsequence$class[mask_dinuc] <- "snvs_dinuc"
  mask_insertions <- (nchar(tmpConsequence$REF)-nchar(tmpConsequence$ALT))>0
  tmpConsequence$class[mask_insertions] <- "indels"
  mask_deletions <- (nchar(tmpConsequence$REF)-nchar(tmpConsequence$ALT))<0
  tmpConsequence$class[mask_deletions] <- "indels"
  
  stopifnot(sum(table(tmpConsequence$class))==dim(tmpConsequence)[1])
  
  tmpConsequence$varCategoryFiltered <- NA
  tmpConsequence$geneSymbol <- unname(sapply(csqList, function(x){
    sapply(strsplit(x,"\\|"), function(y) y[4])
  },simplify=T))
  tmpConsequence$geneEnsemblId <- unname(sapply(csqList, function(x){
    sapply(strsplit(x,"\\|"), function(y) y[5])
  },simplify=T))
  
  tmpConsequence$ips <- ipsFailed
  
  ptv <- c("frameshift", "stop_gained", "splice_acceptor", "splice_donor")
  missense <- c("missense", "start_lost","protein_altering")
  tmpConsequence[!is.na(match(tmpConsequence$varCategoryRaw, ptv)),]$varCategoryFiltered <- "ptv"
  tmpConsequence[!is.na(match(tmpConsequence$varCategoryRaw, missense)),]$varCategoryFiltered <- "missense"
  tmpConsequence[!is.na(match(tmpConsequence$varCategoryRaw, "synonymous")),]$varCategoryFiltered <- "synonymous"
  tmpConsequence[is.na(tmpConsequence$varCategoryFiltered),]$varCategoryFiltered <- "other"
  
  table(tmpConsequence$varCategoryFiltered, tmpConsequence$varCategoryRaw)
  
  ## To determine damaging missense variants, we need to use the threshold of CADD Phred > 15. These frozen CADD
  ## score values are only available for single-point variants.
  ## Consequently, the PTV VS damaging missense variants comparison will only include single-point variants.
  
  mask_missense_snv <- tmpConsequence$varCategoryFiltered=="missense" & nchar(tmpConsequence$REF)==1 & nchar(tmpConsequence$ALT)==1
  loc <- tmpConsequence[mask_missense_snv,]$var
  loc2 <- GRanges(seqnames=sapply(strsplit(loc,"\\:"), function(x) x[1]),
                  IRanges(start=as.numeric(gsub("_.+","",gsub(".+:","",loc))),
                          width=width(tmpConsequence[mask_missense_snv,]$REF)))
  names(loc2) <- loc
  seqlevelsStyle(loc2) <- "UCSC"
  ##Need to download the CADD socres for GRCh37 in inputTabs folder
  bigwig_path <- "inputTabs/CADD_GRCh37-v1.6.bw"
  import_bigwig_region <- import(con = bigwig_path, which = loc2)
  hits <- findOverlaps(loc2, import_bigwig_region)
  loc2$cadd_bigwig <- NA
  loc2$cadd_bigwig <- import_bigwig_region[subjectHits(hits)]$score
  loc2 <- loc2[loc2$cadd_bigwig>15,]
  tmpConsequence[match(names(loc2), tmpConsequence$var),]$varCategoryFiltered <- "missPatho"
  
  table(tmpConsequence$varCategoryFiltered)
  stopifnot(all(tmpConsequence$geneSymbol %in% geneUniverse$symbol))
  
  return(tmpConsequence)
  
}


annotFail <- annotateExclusiveMutations(vcfFileFailed_exclusive, ipsFailed)
annotSucc <- annotateExclusiveMutations(vcfFileSucc_exclusive, ipsSucc)


##* all includes those variants that change significantly the frequency in culture 
##donor ptv_F del_F syn_F ptv_S del_S syn_S

mutCounter <- function(annotFail, varCategory="ptv", outcome="F"){
  
  if (varCategory=="del"){
    
    annotFail$varCategoryFiltered2 <- annotFail$varCategoryFiltered
    annotFail[annotFail$varCategoryFiltered2=="missPatho" | annotFail$varCategoryFiltered2=="ptv",]$varCategoryFiltered2 <- "del"
    
    varTmp <- subset(annotFail, varCategoryFiltered2==varCategory)
    mtrix <- matrix(0, nrow=dim(geneUniverse)[1], ncol=1)
    rownames(mtrix) <- geneUniverse$symbol
    mtrix[match(names(table(varTmp$geneSymbol)), rownames(mtrix)),] <- unname(table(varTmp$geneSymbol))
    colnames(mtrix) <- paste0(varCategory,"_",outcome)
    
  } else {
    
    varTmp <- subset(annotFail, varCategoryFiltered==varCategory)
    mtrix <- matrix(0, nrow=dim(geneUniverse)[1], ncol=1)
    rownames(mtrix) <- geneUniverse$symbol
    mtrix[match(names(table(varTmp$geneSymbol)), rownames(mtrix)),] <- unname(table(varTmp$geneSymbol))
    colnames(mtrix) <- paste0(varCategory,"_",outcome)
    
  }
  
  return(mtrix)
  
}


replicateInfo <- cbind(mutCounter(annotFail, varCategory="ptv", outcome="F"),
                       mutCounter(annotFail, varCategory="del", outcome="F"),
                       mutCounter(annotFail, varCategory="synonymous", outcome="F"),
                       mutCounter(annotSucc, varCategory="ptv", outcome="S"),
                       mutCounter(annotSucc, varCategory="del", outcome="S"),
                       mutCounter(annotSucc, varCategory="synonymous", outcome="S"))

replicateInfo <- as.data.frame(replicateInfo)
      

summTabs <- readRDS("outputTabs/iPSC/discordant/summTabs.RDS")
summTabs <- summTabs[sapply(summTabs, function(x) all(names(x) %in% c(ipsFailed, ipsSucc)))]
parentalDonor <- names(summTabs)

pathToTab <- "outputTabs/iPSC/discordant/outRDS/"
saveRDS(replicateInfo, file=paste0(pathToTab, parentalDonor,"_exclusiveMut.RDS"))




