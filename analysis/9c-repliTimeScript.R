library(ggplot2)
library(ggsignif)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

repliFile <- import("inputTabs/RepliSeq_out.bed")
test <- import("inputTabs/Homo_sapiens.GRCh37.87.gff3.gz")
genes <- test[test$type=="gene"]
genes <- genes[genes$biotype=="protein_coding"]

allvariants <- readRDS("outputTabs/allvariants_v2.RDS")
allvariants <- allvariants[allvariants$gene_name %in% genes$Name,]
##Remove variants that do not fall within the genes of the annotation
tmpwes <- GRanges(allvariants$chrom, IRanges(allvariants$pos, width=nchar(allvariants$ALT)), strand="*")
seqlevelsStyle(tmpwes) <- "NCBI"
hitsGenes <- findOverlaps(tmpwes, genes)
stopifnot(length(unique(queryHits(hitsGenes)))==length(tmpwes))

genesMut <- genes[unique(subjectHits(hitsGenes))]
genesUnMut <- genes[-unique(subjectHits(hitsGenes))]

genesMut <- genesMut[as.logical(seqnames(genesMut) %in% levels(seqnames(tmpwes)))]
genesUnMut <- genesUnMut[as.logical(seqnames(genesUnMut) %in% levels(seqnames(tmpwes)))]


## compute the gene average (Repliscore) per timePoint

seqlevelsStyle(repliFile) <- "NCBI"

## mutated genes
hitsMut <- findOverlaps(genesMut, repliFile)
genesMut$avgRepliSeq <- NA

stopifnot(all(1:length(genesMut) %in% unique(queryHits(hitsMut))))
new_genesMut <- sapply(unique(queryHits(hitsMut)), function(x){
  
  genesMut[x,]$avgRepliSeq <- mean(repliFile[subjectHits(hitsMut[queryHits(hitsMut)==x,]),]$score)
  genesMut[x,]
  
}, simplify=F)

new_genesMut <- do.call("c", new_genesMut)


## unmutated genes

hitsUnMut <- findOverlaps(genesUnMut, repliFile)
genesUnMut$avgRepliSeq <- NA

stopifnot(all(1:length(genesUnMut) %in% unique(queryHits(hitsUnMut))))
new_genesUnMut <- sapply(unique(queryHits(hitsUnMut)), function(x){
  
  genesUnMut[x,]$avgRepliSeq <- mean(repliFile[subjectHits(hitsUnMut[queryHits(hitsUnMut)==x,]),]$score)
  genesUnMut[x,]
  
}, simplify=F)

new_genesUnMut <- do.call("c", new_genesUnMut)

length(genesUnMut)
length(new_genesUnMut)

saveRDS(new_genesMut, file="outputTabs/genesMutAnnotRepli.RDS")
saveRDS(new_genesUnMut, file="outputTabs/genesUnMutAnnotRepli.RDS")


