library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library(qvalue)
library(ggpubr)


## Download Cosmic database, version 90 (GRCh37)
cat( "Annotating COSMIC Tier 1 gene-level information \n")
path_to_cosmic_db <- "inputTabs/CosmicMutantExportCensus.tsv"
cosmic_db <- read.table(path_to_cosmic_db, sep="\t", header=TRUE)
cosmic_db <- cosmic_db[cosmic_db$Tier==1,]
cosmic_db <- as.character(unique(cosmic_db$Gene.name))

## Download DDG2P genes (24.1.2020)
cat( "Annotating DDD gene-level curated list \n")
path_to_ddd_genes <- "inputTabs/DDG2P_24_1_2020.tsv"
ddd_curated <- read.csv(path_to_ddd_genes, sep="\t")
ddd_dominantMOI <- as.character(unique(ddd_curated[ddd_curated$allelic.requirement=="monoallelic",]$gene.symbol))
ddd_curated <- as.character(unique(ddd_curated$gene.symbol))


## Genomic gene length ##
genes <- readRDS("analysis/outputTabs/geneGenomicRanges_Ensembl_v82_gff3.RDS")
path_to_bait <-"inputTabs/Rouhani_2021/Agilent_human_exome_v5_S04380110/S04380110_Covered.baits.nochr.nr.bed"
baitFile <- import(path_to_bait)
hits <- findOverlaps(genes, baitFile)

##Remove genes that do not fall in the exome baits
id_hits <-  1:length(genes)
id_hits_torm <- id_hits[is.na(match(id_hits, unique(queryHits(hits))))]
genes <- genes[-id_hits_torm]
maxCdsxGene <- data.frame(gene=genes$gene_id,
                          genomic_length=width(genes),
                          symbol=genes$Name)

## Keep only longest genes when symbol ids are duplicated (only 75)
duplicated <- names(split(maxCdsxGene$gene, maxCdsxGene$symbol)[elementNROWS(split(maxCdsxGene$gene, maxCdsxGene$symbol))>1])
duplicated_to_rm <- sapply(duplicated, function(x) maxCdsxGene[maxCdsxGene$symbol==x,][which.min(maxCdsxGene[maxCdsxGene$symbol==x,]$genomic_length),]$gene, simplify=TRUE)
maxCdsxGene <- maxCdsxGene[!maxCdsxGene$gene %in% duplicated_to_rm,]

dim(maxCdsxGene)
# [1] 19653     3

## Download Shet scores from Weghorn
shet <- read.table("inputTabs/Shet_Weghorn.txt", header=TRUE)
## Download GnomAD lof metrics per gene (version 2.1.1)
pli <- read.table("inputTabs/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t", header=TRUE)
table(shet$s_het_drift>0.15)
mask_shet <- shet$s_het_drift>0.15
table(pli$pLI>0.9)
mask_pli <- pli$pLI>0.9

maxCdsxGene$ddd <-  maxCdsxGene$symbol %in% ddd_curated
maxCdsxGene$cosmic <-  maxCdsxGene$symbol %in% cosmic_db
maxCdsxGene$ddd_dominantMOI <- maxCdsxGene$symbol %in% ddd_dominantMOI
maxCdsxGene$shet_drift <-  maxCdsxGene$symbol %in% shet[mask_shet,]$Gene
maxCdsxGene$pli <-  maxCdsxGene$symbol %in% pli[mask_pli,]$gene

sapply(maxCdsxGene[,4:8], function(x) table(x))
# ddd cosmic ddd_dominantMOI shet_drift   pli
# FALSE 17715  19095           18961      17670 16593
# TRUE   1938    558             692       1983  3060

saveRDS(maxCdsxGene, file="outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
