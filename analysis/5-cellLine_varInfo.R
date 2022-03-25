library(VariantAnnotation)
library(dplyr)
library(tidyverse)

## We download the VCF files with WES data from 832 HipSci lines
## They are processed one at a time
## The files are annotated with the most severe VEP consequence

args = commandArgs(trailingOnly=TRUE)
fileVcf <- args[1]
ipsLine <- gsub("\\..+","",sapply(strsplit(fileVcf,"\\/"), function(x) x[length(x)]))

geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
vcfFile <- readVcf(fileVcf, genome="hg19")


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
#names(consequence_bcftools) <- gsub("UTR","utr",names(consequence_bcftools))

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

tmpConsequence$ips <- ipsLine

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
bigwig_path <- "inputTabs/CADD_GRCh37-v1.6.bw"
import_bigwig_region <- import(con = bigwig_path, which = loc2)
hits <- findOverlaps(loc2, import_bigwig_region)
loc2$cadd_bigwig <- NA
loc2$cadd_bigwig <- import_bigwig_region[subjectHits(hits)]$score
loc2 <- loc2[loc2$cadd_bigwig>15,]
tmpConsequence[match(names(loc2), tmpConsequence$var),]$varCategoryFiltered <- "missPatho"

table(tmpConsequence$varCategoryFiltered)
stopifnot(all(tmpConsequence$geneSymbol %in% geneUniverse$symbol))


mtrix <- matrix(NA, nrow=nrow(geneUniverse), ncol = 7)
rownames(mtrix) <- geneUniverse$symbol
colnames(mtrix) <- c("All","missPatho","ptv","missense","synonymous","other","del") 


myTabGenerator <- function(tmpConsequence, label="All"){
  
  if (label!="All"){
    print(label)
    tmpConsequence2 <- subset(tmpConsequence, varCategoryFiltered==label)
    
  } else {
    print("All")
    tmpConsequence2 <- tmpConsequence
    
  }
  
  burden <- table(tmpConsequence2$geneSymbol)
  
  index_mtch <- match(label, colnames(mtrix))
  mtrix[,index_mtch] <- unname(burden[match(rownames(mtrix), names(burden))])
  mtrix[,index_mtch][is.na(mtrix[,index_mtch])] <- 0
  return(mtrix[,index_mtch])
}

mtrix <- cbind(myTabGenerator(tmpConsequence, label="All"),
               myTabGenerator(tmpConsequence, label="missPatho"),
               myTabGenerator(tmpConsequence, label="ptv"),
               myTabGenerator(tmpConsequence, label="missense"),
               myTabGenerator(tmpConsequence, label="synonymous"),
               myTabGenerator(tmpConsequence, label="other"))
colnames(mtrix) <- c("All","missPatho","ptv","missense","synonymous","other")   
mtrix_df <- as.data.frame(mtrix)
mtrix_df$del <- mtrix_df$missPatho+mtrix_df$ptv
mtrix_df$geneSymbol <- rownames(mtrix_df)
rownames(mtrix_df) <- NULL

correct_order <- c("geneSymbol", "All", "missPatho","ptv","del","missense","synonymous","other")
mtrix_df <- mtrix_df[,match(correct_order, colnames(mtrix_df))]
fileToPath="outputTabs/iPSC/outRDS/"
saveRDS(mtrix_df, 
        file=paste0(fileToPath,"mutTable.",ipsLine,".RDS"))



