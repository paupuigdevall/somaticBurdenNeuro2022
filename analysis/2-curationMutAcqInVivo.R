library(ggplot2)
library(ggsignif)
library(GenomicScores)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(cowplot)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(nnet)

options(stringsAsFactors = FALSE)

################
#### EXOMES ####
################

########################################################
#### Exomes - SNV point mutations and dinucleotides #### 
########################################################

path_to_directory <- "inputTabs/Rouhani_2021/"
wes.SNVs.637 <- read.table(paste0(path_to_directory,"wes.SNVs.637bams.2019-09-02.txt.gz"),
                           comment.char = "#",
                           stringsAsFactors = FALSE)
colnames(wes.SNVs.637) <- c("chrom", "pos", "fibro", "ips", "Fisher_test_pval",
                            "REF", "ALT", "nREF_fibro", "nALT_fibro", "nREF_ips",
                            "nALT_ips", "AF_1000GP", "AF_ExAC", "clinsing", "SIFT_prob_deleterious",
                            "PolyPhen_prob_deleterious", "gene_name", "ENSID", "bcftools_csq", "VEP",
                            "FILTER")

#chromosome X was run separately
wes.SNVs.637.chrX <- read.table(file=paste0(path_to_directory,"wes.SNVs.637bams.chrX.2021-01-06.reAnnot.txt"), stringsAsFactors = FALSE, sep="\t", header=T)
wes.SNVs.637.chrX$AF_1000GP <- sapply(strsplit(wes.SNVs.637.chrX$AF_1000GP,","), function(x) max(x))
wes.SNVs.637.chrX$AF_ExAC <- sapply(strsplit(wes.SNVs.637.chrX$AF_ExAC,","), function(x) max(x))
wes.SNVs.637.chrX[wes.SNVs.637.chrX$bcftools_csq=="",]$bcftools_csq <- "."

wes.SNVs.637 <- rbind(wes.SNVs.637, wes.SNVs.637.chrX)

## Check if all variants fall within the covered baits
path_to_bait <- paste0(path_to_directory,"Agilent_human_exome_v5_S04380110/S04380110_Covered.baits.nochr.nr.bed")
baitFile <- import(path_to_bait)
tmpGR <- GRanges(wes.SNVs.637$chrom, IRanges(start=wes.SNVs.637$pos, width=1))
width(tmpGR[nchar(wes.SNVs.637$REF)==2]) <- 2
stopifnot(length(findOverlaps(tmpGR, baitFile))==length(tmpGR))

## Create dinucleotides column
wes.SNVs.637$dinucleotides <- NA
wes.SNVs.637[(nchar(wes.SNVs.637$REF)==2 & nchar(wes.SNVs.637$ALT)==2),]$dinucleotides <- "yes"
wes.SNVs.637[!(nchar(wes.SNVs.637$REF)==2 & nchar(wes.SNVs.637$ALT)==2),]$dinucleotides <- "no"

#############################
### Filtering of variants ###
#############################

## Germline mutation (>0.1% MAF) are removed 
wes.SNVs.637$AF_1000GP <- as.numeric(wes.SNVs.637$AF_1000GP)
wes.SNVs.637$AF_ExAC <- as.numeric(wes.SNVs.637$AF_ExAC)
wes.SNVs.637 <- wes.SNVs.637[-sort(unique(c(which(wes.SNVs.637$AF_ExAC>0.001), which(wes.SNVs.637$AF_1000GP>0.001)))),]

## Mutation calls are only considered as a significant change in allele frequency between the fibroblast and IPS cell sample
wes.SNVs.637 <- wes.SNVs.637[wes.SNVs.637$Fisher_test_pval<1.6e-4,]

## Filtration of probable spurious mutation calls (allele frequency cutoff >0.6)
mask_06_ips <- wes.SNVs.637$nALT_ips/(wes.SNVs.637$nREF_ips+wes.SNVs.637$nALT_ips)>0.6
wes.SNVs.637 <- wes.SNVs.637[!mask_06_ips,]
mask_06_fibro <- wes.SNVs.637$nALT_fibro/(wes.SNVs.637$nREF_fibro+wes.SNVs.637$nALT_fibro)>0.6
wes.SNVs.637 <- wes.SNVs.637[!mask_06_fibro,]

## Select variants that are acquired or positively selected in iPSC respect to fibroblasts (ratio>1)
mask_ratio <- ((wes.SNVs.637$nALT_ips/(wes.SNVs.637$nREF_ips+wes.SNVs.637$nALT_ips))/(wes.SNVs.637$nALT_fibro/(wes.SNVs.637$nREF_fibro+wes.SNVs.637$nALT_fibro)))>1
wes.SNVs.637 <- wes.SNVs.637[mask_ratio,]

## Only consider high-quality variants ("PASS" filter)
wes.SNVs.637 <- wes.SNVs.637[wes.SNVs.637$FILTER=="PASS",]

## Remove variants found in more than one donor (Petr's dataset already filtered for that)
wes.SNVs.637$varcolumn <- paste0(wes.SNVs.637$chrom,
                                 "-",
                                 wes.SNVs.637$pos,
                                 "-",
                                 wes.SNVs.637$REF,
                                 "-",
                                 wes.SNVs.637$ALT)

stopifnot(all(sapply(unique(wes.SNVs.637[duplicated(wes.SNVs.637$varcolumn),]$varcolumn),
                     function(x) elementNROWS(unique(subset(wes.SNVs.637, wes.SNVs.637$varcolumn==x)$fibro)))==1))

wes.SNVs.637$varcolumn <- NULL
summary(as.numeric(table(wes.SNVs.637$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   24.00   37.50   66.92   74.25  806.00

summary(as.numeric(table(wes.SNVs.637[wes.SNVs.637$dinucleotides=="yes",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    6.00   15.65   17.00  206.00 

summary(as.numeric(table(wes.SNVs.637[wes.SNVs.637$dinucleotides=="no",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   23.00   35.00   56.08   64.25  600.00 


## Annotating hyper-mutated lines (Z-score>2)
lines_zscores <- (as.numeric(table(wes.SNVs.637$ips))-mean(as.numeric(table(wes.SNVs.637$ips))))/(sd(as.numeric(table(wes.SNVs.637$ips))))
names(lines_zscores) <- names(table(wes.SNVs.637$ips))

hyper <- names(lines_zscores)[abs(lines_zscores)>2]
# [1] "HPSI0114i-joxm_1" "HPSI0214i-eiwy_1" "HPSI0413i-ougl_1" "HPSI0413i-ougl_3" "HPSI0413i-uahf_3" "HPSI0414i-ceik_1" "HPSI0414i-uawq_2" "HPSI0613i-funp_2" "HPSI0614i-ciwj_1"
# [10] "HPSI0614i-ciwj_2" "HPSI0813i-wots_3" "HPSI0814i-bokz_5" "HPSI0814i-bokz_6" "HPSI0913i-bulb_1" "HPSI0914i-zerv_8" "HPSI1013i-garx_2" "HPSI1014i-tixi_4"

## Remove hyper-mutated lines
wes.SNVs.637 <- wes.SNVs.637[!wes.SNVs.637$ips %in% hyper,]

summary(as.numeric(table(wes.SNVs.637$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.0    23.5    37.0    52.9    67.0   235.0 

summary(as.numeric(table(wes.SNVs.637[wes.SNVs.637$dinucleotides=="yes",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    5.00   10.43   15.00   66.00 

summary(as.numeric(table(wes.SNVs.637[wes.SNVs.637$dinucleotides=="no",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   23.00   34.00   45.83   58.00  195.00

## Annotating coding and non-coding variants based on location
wes.SNVs.637$coding <- NA


##NOTE! bcftools-csq does not annotate the consequence properly for all coding dinucleotides 
mask_dinuc <- wes.SNVs.637$dinucleotides=="yes"
mask_dinuc_coding <- sapply(strsplit(wes.SNVs.637$bcftools_csq[mask_dinuc],","), function(x) max(elementNROWS(sapply(strsplit(x,"\\|"), function(y) y, simplify=FALSE))))==7
wes.SNVs.637$coding[mask_dinuc][mask_dinuc_coding] <- "yes"
wes.SNVs.637$coding[mask_dinuc][!mask_dinuc_coding] <- "no"
## Frome 1,874 coding dinucleotides, 836 are annotated as single-point variants
table(mask_dinuc_coding)
# mask_dinuc_coding
# FALSE  TRUE 
# 723  1874 
mask_dinuc_coding_failed <- !sapply(strsplit(wes.SNVs.637$bcftools_csq[mask_dinuc][mask_dinuc_coding],","), function(x) any(sapply(strsplit(x,"\\|"), function(y) grepl("\\+",y[7]))))
table(mask_dinuc_coding_failed)
# mask_dinuc_coding_failed
# FALSE  TRUE 
# 1038   836 

## 1,874 out of 2,597 dinucleotides are coding
## 12,907 out of 16,818 single-point variants are coding
mask_snv_coding <- sapply(strsplit(wes.SNVs.637$bcftools_csq[!mask_dinuc],","), function(x) max(elementNROWS(sapply(strsplit(x,"\\|"), function(y) y, simplify=FALSE))))==7
wes.SNVs.637$coding[!mask_dinuc][mask_snv_coding] <- "yes"
wes.SNVs.637$coding[!mask_dinuc][!mask_snv_coding] <- "no"
table(wes.SNVs.637$coding)
# no   yes 
# 6508 12907 

vep_severity <- read.csv("inputTabs/vep_severity.csv", header=TRUE)
consequence <- as.character(1:length(vep_severity$SO_term))
names(consequence) <- gsub("_variant","",vep_severity$SO_term)
undef <- as.character(length(consequence)+1) 
names(undef) <- "."
consequence_bcftools <- c(consequence, undef)
names(consequence_bcftools) <- gsub("UTR","utr",names(consequence_bcftools))


max_consequence <- sapply(strsplit(wes.SNVs.637$bcftools_csq, ","), function(x) unique(sapply(strsplit(x,"\\|"), function(y) sapply(strsplit(gsub("\\*","", gsub("@.+","",y[1])),"&"), function(z) z[1]))), simplify=T)
max_consequence <- sapply(max_consequence, function(x) names(which.min(consequence_bcftools[x[!is.na(x)]])))
max_consequence[elementNROWS(max_consequence)==0] <- NA
wes.SNVs.637$varCategoryRaw <- unlist(max_consequence)
wes.SNVs.637[is.na(wes.SNVs.637$varCategoryRaw),]$varCategoryRaw <- "other"
wes.SNVs.637[wes.SNVs.637$varCategoryRaw==".",]$varCategoryRaw <- "other"


## Annotation of protein-truncating variants (PTV, loss of function) and missense pathogenic variants
wes.SNVs.637$varCategoryFiltered <- NA

## PTV (loss-of-function)
#1-frameshift
#2-stop-gain
#3-splice-acceptor
#4-splice-donor

## missense_pathogenic
#1-missense variants annotation
#2-initiator codon removed (start_lost) annotation
##1&2 must have CADD Phred>15 to be considered damaging missense

ptv <- c("frameshift", "stop_gained", "splice_acceptor", "splice_donor")
missense <- c("missense", "start_lost")
wes.SNVs.637[!is.na(match(wes.SNVs.637$varCategoryRaw, ptv)),]$varCategoryFiltered <- "ptv"
wes.SNVs.637[!is.na(match(wes.SNVs.637$varCategoryRaw, missense)),]$varCategoryFiltered <- "missense"
wes.SNVs.637[!is.na(match(wes.SNVs.637$varCategoryRaw, "synonymous")),]$varCategoryFiltered <- "synonymous"
wes.SNVs.637[is.na(wes.SNVs.637$varCategoryFiltered),]$varCategoryFiltered <- "other"

## To determine damaging missense variants, we need to use the threshold of CADD Phred > 15. These frozen CADD
## score values are only available for single-point variants.
## Consequently, the PTV VS damaging missense variants comparison will only include single-point variants.

mask_missense_snv <- wes.SNVs.637$varCategoryFiltered=="missense" & nchar(wes.SNVs.637$REF)==1 & nchar(wes.SNVs.637$ALT)==1
wes.SNVs.637$varcolumn <- paste0(wes.SNVs.637$ips,
                                 "-",
                                 wes.SNVs.637$chrom,
                                 ":",
                                 wes.SNVs.637$pos,
                                 "-",
                                 wes.SNVs.637$REF,
                                 "-",
                                 wes.SNVs.637$ALT)

loc <- wes.SNVs.637[mask_missense_snv,]$varcolumn

loc2 <- GRanges(sapply(strsplit(loc, "-"), function(x) x[3]))
names(loc2) <- loc
seqlevelsStyle(loc2) <- "UCSC"

##Need to download the CADD socres for GRCh37 in inputTabs folder
bigwig_path <- "inputTabs/CADD_GRCh37-v1.6.bw"
import_bigwig_region <- import(con = bigwig_path, which = loc2)
hits <- findOverlaps(loc2, import_bigwig_region)
loc2$cadd_bigwig <- NA
loc2$cadd_bigwig <- import_bigwig_region[subjectHits(hits)]$score
loc2 <- loc2[loc2$cadd_bigwig>15,]
wes.SNVs.637[match(names(loc2), wes.SNVs.637$varcolumn),]$varCategoryFiltered <- "missPatho"

table(wes.SNVs.637$varCategoryFiltered)
# missense  missPatho      other        ptv synonymous 
# 2779       6069       6473        745       3349 

#########################
#### Exomes - Indels ####
#########################

wes.SNVs.637.indels <- read.table(paste0(path_to_directory,"wes.indels.mpileup.637bams.2019-09-02.txt"),
                                  comment.char = "#",
                                  stringsAsFactors = FALSE)
colnames(wes.SNVs.637.indels) <- c("chrom", "pos", "fibro", "ips", "Fisher_test_pval",
                                   "REF", "ALT", "nREF_fibro", "nALT_fibro", "nREF_ips",
                                   "nALT_ips", "AF_1000GP", "AF_ExAC", "clinsing", "SIFT_prob_deleterious",
                                   "PolyPhen_prob_deleterious", "gene_name", "ENSID", "bcftools_csq", "VEP",
                                   "FILTER")

## Germline mutation are removed >1% MAF
wes.SNVs.637.indels$AF_1000GP <- as.numeric(wes.SNVs.637.indels$AF_1000GP)
wes.SNVs.637.indels$AF_ExAC <- as.numeric(wes.SNVs.637.indels$AF_ExAC)
wes.SNVs.637.indels <- wes.SNVs.637.indels[-sort(unique(c(which(wes.SNVs.637.indels$AF_ExAC>0.001), which(wes.SNVs.637.indels$AF_1000GP>0.001)))),]

## Mutation calls are only considered as a significant change in allele frequency between the fibroblast and IPS cell sample
wes.SNVs.637.indels <- wes.SNVs.637.indels[wes.SNVs.637.indels$Fisher_test_pval<1.6e-4,]

## Filtration of probable spurious mutation calls (allele frequency cutoff >0.6)
mask_06_ips <- wes.SNVs.637.indels$nALT_ips/(wes.SNVs.637.indels$nREF_ips+wes.SNVs.637.indels$nALT_ips)>0.6
wes.SNVs.637.indels <- wes.SNVs.637.indels[!mask_06_ips,]
mask_06_fibro <- wes.SNVs.637.indels$nALT_fibro/(wes.SNVs.637.indels$nREF_fibro+wes.SNVs.637.indels$nALT_fibro)>0.6
wes.SNVs.637.indels <- wes.SNVs.637.indels[!mask_06_fibro,]

## Select variants that are acquired or positively selected in iPSC respect to fibroblasts (ratio>1)
mask_ratio <- ((wes.SNVs.637.indels$nALT_ips/(wes.SNVs.637.indels$nREF_ips+wes.SNVs.637.indels$nALT_ips))/(wes.SNVs.637.indels$nALT_fibro/(wes.SNVs.637.indels$nREF_fibro+wes.SNVs.637.indels$nALT_fibro)))>1
wes.SNVs.637.indels <- wes.SNVs.637.indels[mask_ratio,]

## Only consider high-quality variants ("PASS" filter)
wes.SNVs.637.indels <- wes.SNVs.637.indels[wes.SNVs.637.indels$FILTER=="PASS",]

## Remove variants found in more than one donor (Petr's dataset already filtered for that)
wes.SNVs.637.indels$varcolumn <- paste0(wes.SNVs.637.indels$chrom,
                                        "-",
                                        wes.SNVs.637.indels$pos,
                                        "-",
                                        wes.SNVs.637.indels$REF,
                                        "-",
                                        wes.SNVs.637.indels$ALT)

stopifnot(all(sapply(unique(wes.SNVs.637.indels[duplicated(wes.SNVs.637.indels$varcolumn),]$varcolumn),
                     function(x) elementNROWS(unique(subset(wes.SNVs.637.indels, wes.SNVs.637.indels$varcolumn==x)$fibro)))==1))

wes.SNVs.637.indels$varcolumn <- NULL

summary(as.numeric(table(wes.SNVs.637.indels$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   2.304   3.000   8.000 

## Remove hyper-mutated lines
wes.SNVs.637.indels <- wes.SNVs.637.indels[!wes.SNVs.637.indels$ips %in% hyper,]
summary(as.numeric(table(wes.SNVs.637.indels$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   2.304   3.000   8.000 

## Annotating coding and non-coding variants based on location
wes.SNVs.637.indels$coding <- NA
mask_coding <- sapply(strsplit(wes.SNVs.637.indels$VEP,","), function(x) sum(unique(sapply(strsplit(x, "\\|"), function(y) !("" %in% y[16])))))>0
wes.SNVs.637.indels$coding[mask_coding] <- "yes"
wes.SNVs.637.indels$coding[!mask_coding] <- "no"


max_consequence <- sapply(strsplit(wes.SNVs.637.indels$VEP, ","), function(x) unique(sapply(strsplit(x,"\\|"), function(y) sapply(strsplit(gsub("_variant","", y[2]),"&"), function(z) z[1]))))
max_consequence <- sapply(max_consequence, function(x) x[which.min(as.numeric(consequence[x[!is.na(x)]]))])
wes.SNVs.637.indels$varCategoryRaw <- unlist(max_consequence)
wes.SNVs.637.indels$varCategoryFiltered <- NA

ptv <- c("frameshift", "stop_gained", "splice_acceptor", "splice_donor")
missense <- c("missense", "start_lost")
wes.SNVs.637.indels[!is.na(match(wes.SNVs.637.indels$varCategoryRaw, ptv)),]$varCategoryFiltered <- "ptv"
#wes.SNVs.637.indels[!is.na(match(wes.SNVs.637.indels$varCategoryRaw, missense)),]$varCategoryFiltered <- "missense"
#wes.SNVs.637.indels[!is.na(match(wes.SNVs.637.indels$varCategoryRaw, "synonymous")),]$varCategoryFiltered <- "synonymous"
wes.SNVs.637.indels[is.na(wes.SNVs.637.indels$varCategoryFiltered),]$varCategoryFiltered <- "other"

table(wes.SNVs.637.indels$varCategoryFiltered)
# other   ptv 
# 387   304

table(wes.SNVs.637.indels$coding)
# no yes 
# 315 367

summary(as.numeric(table(wes.SNVs.637.indels[wes.SNVs.637.indels$coding=="yes",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.631   2.000   6.000

summary(as.numeric(table(wes.SNVs.637.indels[wes.SNVs.637.indels$varCategoryFiltered=="ptv",]$ips)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.483   2.000   5.000


#######################
#### Exomes - CNVs #### 
#######################
library(readxl)
metadata_table <- read.table("inputTabs/hipsci.qc1_sample_info.20170927.tsv",
                             header=TRUE, sep="\t")
stopifnot(all(!is.na(match(unique(wes.SNVs.637$ips), metadata_table$name))))
metadata_table <- metadata_table[match(unique(wes.SNVs.637$ips), metadata_table$name),]
select_columns <- c(grep("name", colnames(metadata_table)),
                    grep("derived_from$", colnames(metadata_table)),
                    grep("cnv", colnames(metadata_table)))
cnv_tab <- metadata_table[,select_columns]

cnv_hipsci <- "inputTabs/nature22403-s3.xlsx"
cnv_hipsci <- as.data.frame(read_excel(cnv_hipsci, sheet = 1))
cnv_hipsci <- cnv_hipsci[cnv_hipsci$line %in% gsub(".+-","",cnv_tab$name),]

num_cnvRegions <- data.frame(donor=metadata_table$name,
                             cnv_hipsci=metadata_table$cnv_num_different_regions,
                             cnv_nature_K200=NA,
                             cnv_nature_K500=NA,
                             cnv_nature_K1000=NA)

cnv_hipsci$full <- cnv_tab[match(cnv_hipsci$line, gsub(".+-","",cnv_tab$name)),]$name
cnv_hipsci <- cnv_hipsci %>% filter( CNA_threshold == "K200")


#############################
## QC using the annotation ##
#############################

allvariants <- wes.SNVs.637
allvariants$class <- "snvs_dinuc"

tmp <- wes.SNVs.637.indels
tmp$class <- "indels"

allvariants <- allvariants[,match(intersect(colnames(allvariants), colnames(tmp)), colnames(allvariants))]
tmp <- tmp[,match(intersect(colnames(allvariants), colnames(tmp)), colnames(tmp))]

stopifnot(colnames(allvariants)==colnames(tmp))
allvariants <- rbind(allvariants, tmp)


###############################################################
### Max.CDS per gene (obtained from longest transcript CDS) ###
###############################################################

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)

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
saveRDS(maxCdsxGene, file="outputTabs/maxCdsxGene_UCSC_v2.RDS")

dim(maxCdsxGene)
# [1] 19809     3

## 8,107 out of the 19,809 genes from the gene universe are mutated (geneUniverse using the UCSC annotation)
table(unique(allvariants$gene_name) %in% maxCdsxGene$symbol)
# FALSE  TRUE 
# 1072  8107 


### We determine the max CDS length per gene with the corresponding tx using version 87 of Ensembl Gene Annotation
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
## Max CDS per gene (n=20,180)
maxCdsxGene$tx <- unname(txMaxCds)
saveRDS(maxCdsxGene, file="outputTabs/maxCdsxGene_Ensembl_v82_gff3.RDS")

## Cds by transcript
cdsbytx <- cdsbytx[match(maxCdsxGene$tx, names(cdsbytx))]
saveRDS(cdsbytx, file="outputTabs/cdsbytx_Ensembl_v82_gff3.RDS")


## Gene universe using 87 gff3 version of the GRCh37 human genome (protein-coding genes)
genes <- genes[genes$biotype=="protein_coding"]
saveRDS(genes, file="outputTabs/geneGenomicRanges_Ensembl_v82_gff3.RDS")
length(genes)
# [1] 20356

table(unique(allvariants$gene_name) %in% genes$Name)
# FALSE  TRUE 
# 571  8608

### The variants from the genes that are not accounted by the used annotation in downstream analysis are filtered-out
##Remove variants in genes not present in the annotation
allvariants <- allvariants[allvariants$gene_name %in% genes$Name,]

##Remove variants that do not fall within the genes of the annotation
tmpwes <- GRanges(allvariants$chrom, IRanges(allvariants$pos, width=nchar(allvariants$ALT)), strand="*")
hits2 <- findOverlaps(tmpwes, genes)
stopifnot(length(unique(queryHits(hits2)))==length(tmpwes))
saveRDS(allvariants, file="outputTabs/allvariants_v2.RDS")
dim(allvariants)
# [1] 18539    25


#############################
###### Summary table ########
#############################

### All burden
exomes_allinfo <- as.data.frame(table(allvariants$ips))
colnames(exomes_allinfo) <- c("donor_iPSC", "all")

### Synonymous
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="synonymous",]$ips))
exomes_allinfo$synonymous <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$synonymous))){
  exomes_allinfo[is.na(exomes_allinfo$synonymous),]$synonymous <- 0
}

### Missense (pathogenic and non-pathogenic)
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="missense" | allvariants$varCategoryFiltered=="missPatho",]$ips))
exomes_allinfo$missense <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$missense))){
  exomes_allinfo[is.na(exomes_allinfo$missense),]$missense <- 0
}

### Deleterious (Missense pathogenic + ptv)
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="ptv" | allvariants$varCategoryFiltered=="missPatho",]$ips))
exomes_allinfo$deleterious <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$deleterious))){
  exomes_allinfo[is.na(exomes_allinfo$deleterious),]$deleterious <- 0
}

### Missense pathogenic (all)
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="missPatho",]$ips))
exomes_allinfo$missPatho <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$missPatho))){
  exomes_allinfo[is.na(exomes_allinfo$missPatho),]$missPatho <- 0
}

### Protein-truncating variants (ptv)
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="ptv",]$ips))
exomes_allinfo$ptv <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$ptv))){
  exomes_allinfo[is.na(exomes_allinfo$ptv),]$ptv <- 0
}

### Other
tmp <- as.data.frame(table(allvariants[allvariants$varCategoryFiltered=="other",]$ips))
exomes_allinfo$other <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$other))){
  exomes_allinfo[is.na(exomes_allinfo$other),]$other <- 0
}

### SNV (including dinuc)
tmp <- as.data.frame(table(allvariants[nchar(allvariants$REF)==2 & nchar(allvariants$ALT)==2 & allvariants$class=="snvs_dinuc",]$ips))
exomes_allinfo$dinuc <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$dinuc))){
  exomes_allinfo[is.na(exomes_allinfo$dinuc),]$dinuc <- 0
}

### SNV (without dinuc)
tmp <- as.data.frame(table(allvariants[nchar(allvariants$REF)==1 & nchar(allvariants$ALT)==1 & allvariants$class=="snvs_dinuc",]$ips))
exomes_allinfo$onlySnvs <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$onlySnvs))){
  exomes_allinfo[is.na(exomes_allinfo$onlySnvs),]$onlySnvs <- 0
}

### indels
tmp <- as.data.frame(table(allvariants[allvariants$class=="indels",]$ips))
exomes_allinfo$indels <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$indels))){
  exomes_allinfo[is.na(exomes_allinfo$indels),]$indels <- 0
}

### indels (ptv)
tmp <- as.data.frame(table(allvariants[allvariants$class=="indels" & allvariants$varCategoryFiltered=="ptv",]$ips))
exomes_allinfo$indels_ptv <- tmp[match(exomes_allinfo$donor_iPSC, tmp$Var1),]$Freq
if (sum(is.na(exomes_allinfo$indels_ptv))){
  exomes_allinfo[is.na(exomes_allinfo$indels_ptv),]$indels_ptv <- 0
}

## CNV
tmp <- data.frame(name=cnv_tab$name,
                  cnv_num_different_regions=cnv_tab$cnv_num_different_regions,
                  cnv_length_different_regions_Mbp=cnv_tab$cnv_length_different_regions_Mbp,
                  cnv_length_shared_differences_Mbp=cnv_tab$cnv_length_shared_differences_Mbp)
tmp[is.na(tmp)] <- 0
exomes_allinfo <- cbind(exomes_allinfo, tmp[match(exomes_allinfo$donor_iPSC, tmp$name),-1])           
saveRDS(exomes_allinfo, file="outputTabs/exomes_allinfo_nonhyper_v2.RDS")


##########################################################
### Table of mutational counts per gene (TOTAL BURDEN) ###
##########################################################

samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
mutationsGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(mutationsGeneSample) <- genes
colnames(mutationsGeneSample) <- samples
tmp <- sapply(samples, function(x) table(subset(allvariants, allvariants$ips==x)$gene_name))
for (i in 1:length(tmp)){
  mutationsGeneSample[match(names(tmp[[i]]), rownames(mutationsGeneSample)),
                      match(names(tmp)[i], colnames(mutationsGeneSample))] <- tmp[[i]]
}

stopifnot(sum(colSums(mutationsGeneSample))==sum(exomes_allinfo$all))
write.table(mutationsGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/geneMutBurden_v2.txt")



###############################################
### Table of mutational counts per SNVs-ptv ###
###############################################

ptvVariants <- subset(allvariants, varCategoryFiltered=="ptv")
samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
ptvGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(ptvGeneSample) <- genes
colnames(ptvGeneSample) <- unique(wes.SNVs.637$ips)

tmp <- sapply(samples, function(x) table(subset(ptvVariants, ptvVariants$ips==x)$gene_name))

for (i in 1:length(tmp)){
  if (length(tmp[[i]])>0){
    ptvGeneSample[match(names(tmp[[i]]), rownames(ptvGeneSample)),
                  match(names(tmp)[i], colnames(ptvGeneSample))] <- tmp[[i]]
  }
}

stopifnot(sum(colSums(ptvGeneSample))==sum(exomes_allinfo$ptv))
write.table(ptvGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/ptvMutBurden_v2.txt")

##########################################################
### Table of mutational counts per pathogenic missense ###
##########################################################

missPathoVariants <- subset(allvariants, varCategoryFiltered=="missPatho")
samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
missGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(missGeneSample) <- genes
colnames(missGeneSample) <- unique(wes.SNVs.637$ips)

tmp <- sapply(samples, function(x) table(subset(missPathoVariants, missPathoVariants$ips==x)$gene_name))

for (i in 1:length(tmp)){
  if (length(tmp[[i]])>0){
    missGeneSample[match(names(tmp[[i]]), rownames(missGeneSample)),
                  match(names(tmp)[i], colnames(missGeneSample))] <- tmp[[i]]
  }
}
stopifnot(sum(colSums(missGeneSample))==sum(exomes_allinfo$missPatho))
write.table(missGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/missMutBurden_v2.txt")




#############################################################
### Table of mutational counts per "deleterious" variants ###
#############################################################

delVariants <- subset(allvariants, varCategoryFiltered=="missPatho" | varCategoryFiltered=="ptv")
samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
delGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(delGeneSample) <- genes
colnames(delGeneSample) <- unique(allvariants$ips)

tmp <- sapply(samples, function(x) table(subset(delVariants, delVariants$ips==x)$gene_name))

for (i in 1:length(tmp)){

  if (length(tmp[[i]])>0){
    delGeneSample[match(names(tmp[[i]]), rownames(delGeneSample)),
                  match(names(tmp)[i], colnames(delGeneSample))] <- tmp[[i]]
  }
  
}

stopifnot(sum(colSums(delGeneSample))==sum(exomes_allinfo$deleterious))


write.table(delGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/delMutBurden_v2.txt")



##########################################################
### Table of mutational counts per synonymous-variants ###
##########################################################

synVariants <- subset(allvariants, varCategoryFiltered=="synonymous")
samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
synGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(synGeneSample) <- genes
colnames(synGeneSample) <- unique(allvariants$ips)

tmp <- sapply(samples, function(x) table(subset(synVariants, synVariants$ips==x)$gene_name))

for (i in 1:length(tmp)){

  if (length(tmp[[i]])>0){
    synGeneSample[match(names(tmp[[i]]), rownames(synGeneSample)),
                  match(names(tmp)[i], colnames(synGeneSample))] <- tmp[[i]]
  }
  
  
}

stopifnot(sum(colSums(synGeneSample))==sum(exomes_allinfo$synonymous))

write.table(synGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/synMutBurden_v2.txt")



###################################################################################################
### Table of mutational counts per Other variants (mainly intronic, 5' & 3'utr + splice region) ###
###################################################################################################


otherVariants <- subset(allvariants, varCategoryFiltered=="other")
samples <- unique(allvariants$ips)
genes <- unique(allvariants$gene_name)
othersGeneSample <- matrix(0, nrow = length(genes), ncol = length(samples))
rownames(othersGeneSample) <- genes
colnames(othersGeneSample) <- unique(allvariants$ips)

tmp <- sapply(samples, function(x) table(subset(otherVariants, otherVariants$ips==x)$gene_name))

for (i in 1:length(tmp)){

  if (length(tmp[[i]])>0){
    othersGeneSample[match(names(tmp[[i]]), rownames(othersGeneSample)),
                     match(names(tmp)[i], colnames(othersGeneSample))] <- tmp[[i]]
  }


}

stopifnot(sum(colSums(othersGeneSample))==sum(exomes_allinfo$other))

write.table(othersGeneSample,
            sep="\t",
            quote=FALSE,
            file="outputTabs/othersMutBurden_v2.txt")




######################################
#### Summary figures in the paper ####
######################################


### Mutational burden summary ###
allburden <- exomes_allinfo[,c("all","cnv_num_different_regions")]
allburden <- rowSums(allburden)
allburden <- allburden[!is.na(allburden)]
sum(allburden)
# [1] 18999

summary(allburden)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   24.00   37.00   51.77   65.00  214.00 


allburden_exceptCNV <- exomes_allinfo[,c("all")]
allburden_exceptCNV <- allburden_exceptCNV[!is.na(allburden_exceptCNV)]
summary(allburden_exceptCNV)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   23.00   35.00   50.51   63.50  213.00 


allburden_del <- exomes_allinfo[,c("deleterious")]
allburden_del <- allburden_del[!is.na(allburden_del)]
summary(allburden_del)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    9.00   14.00   18.32   24.00   70.00 

(sum(exomes_allinfo$missense)+sum(exomes_allinfo$ptv))*100/sum(exomes_allinfo$all)
# [1] 50.46119
