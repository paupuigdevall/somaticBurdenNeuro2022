
library(ggplot2)
library(gridExtra)
library(scales)
library(ggpubr)
library(viridis)

#####################################
### QUALITY CONTROL for all pools ###
#####################################

################
### flagstat ###
################

pathTOfile <- "../demuxlet/deconvolution/"

bamQC <- read.table(file=paste0(pathTOfile,"bam_flagstatQC.txt"), header=F, sep="\t")
bamQC <- data.frame(pathToBamFile=bamQC$V1[grepl(".bam", bamQC$V1)],
                 perc_mapped=gsub(" ", "",bamQC$V1[!grepl(".bam", bamQC$V1)]))
stopifnot(all(grepl("Pool", sapply(strsplit(bamQC$pathToBamFile, "/"), function(x) x[11]))))
poolnames <- sapply(strsplit(bamQC$pathToBamFile, "/"), function(x) x[11])

bamQC$protocol <- NA
bamQC$protocol[grepl("D119", poolnames)] <- "organoids"
bamQC$protocol[!grepl("D119", poolnames)] <- "da_neurons"

bamQC$pool <- gsub("[a-b]-.+", "", poolnames)

mask_biorep1 <- !grepl("BIO2", poolnames)
mask_biorep2 <- grepl("BIO2", poolnames)
bamQC$bioRep <- NA
bamQC$bioRep[mask_biorep1] <- 1
bamQC$bioRep[mask_biorep2] <- 2

mask_techrepA <- grepl("a-", poolnames)
mask_techrepB <- grepl("b-", poolnames)
bamQC$techRep <- NA
bamQC$techRep[mask_techrepA] <- 1
bamQC$techRep[mask_techrepB] <- 2

mask_tenXRep1 <- !grepl("REP2", poolnames)
mask_tenXRep2 <- grepl("REP2", poolnames)
bamQC$tenXRep <- NA
bamQC$tenXRep[mask_tenXRep1] <- 1
bamQC$tenXRep[mask_tenXRep2] <- 2


bamQC$timepoint <- gsub("ROT", "", sapply(strsplit(poolnames, "-") , function(x) x[grepl("^D.+", x)]))

bamQC$rotenone_treatment <- NA
bamQC$rotenone_treatment[grepl("ROT", poolnames)] <- "yes"
bamQC$rotenone_treatment[!grepl("ROT", poolnames)] <- "no"


bamQC$sampleId <- sapply(strsplit(bamQC$pathToBamFile, "/"), function(x) x[12])
bamQC$runId <- sapply(strsplit(bamQC$pathToBamFile, "/"), function(x) x[10])


library(ggplot2)
bamQCmod <- subset(bamQC, bamQC$timepoint!="D119")
bamQCmod$pool <- factor(bamQCmod$pool, levels=paste0("Pool", c(1:17,20,21)))
bamQCmod$timepoint <- factor(bamQCmod$timepoint, levels=paste0("D", c(11,30,52,119)))
bamQCmod$perc_mapped <- as.numeric(gsub("%", "", as.character(bamQCmod$perc)))

########################
### QC summary table ###
########################

#totalnum_reads
pathTofile2 <- "../demuxlet/deconvolution/log/"
system(paste0("grep 'Total number input reads :' ", pathTofile2, "demuxlet.*.err | cut -d':' -f1,5 > ", pathTofile2,"total_number_reads.txt"))
totalnum_reads <- read.table(paste0(pathTofile2, "total_number_reads.txt"))

#droplet
system(paste0("grep 'droplet/cell barcodes to consider' ", pathTofile2, "demuxlet.*.err | awk '{print $1,$7}' | sed 's/NOTICE//g' > ", pathTofile2,"droplet_num.txt"))
droplet_num <- read.table(paste0(pathTofile2, "droplet_num.txt"))

#totalnum_passfilt_reads
system(paste0("grep 'Total number of pass-filtered reads :' ", pathTofile2, "demuxlet.*.err | cut -d':' -f1,5 > ", pathTofile2,"total_passfilt_reads.txt"))
passfilt_reads <- read.table(paste0(pathTofile2, "total_passfilt_reads.txt"))

#totalnum_markers
system(paste0("grep 'Total number valid SNPs observed' ", pathTofile2, "demuxlet.*.err | cut -d':' -f1,5 > ", pathTofile2,"valid_SNPs.txt"))
valid_snps <- read.table(paste0(pathTofile2, "valid_SNPs.txt"))

totalnum_reads$V1 <- gsub(":", "", gsub(".+log\\/", "", totalnum_reads$V1))
droplet_num$V1 <- gsub(":", "", gsub(".+log\\/", "", droplet_num$V1))
passfilt_reads$V1 <- gsub(":", "", gsub(".+log\\/", "", passfilt_reads$V1))
valid_snps$V1 <- gsub(":", "", gsub(".+log\\/", "", valid_snps$V1))

allfiles <- read.table(paste0(pathTOfile,"pathToPossortedBam_And_PoolID_per10xSample.txt"))
allfiles$fileID <- gsub("/projects/NeuroSeq/scRNAseqData/fastq/", "", allfiles$V1)
allfiles$fileID <- gsub("/cellranger-hg19/outs/possorted_genome_bam.bam", "", allfiles$fileID)
allfiles$index <- as.character(1:length(allfiles$V1))

#totalnum_reads
allfiles$totalnum_reads <- NA_integer_
allfiles[match(sapply(strsplit(totalnum_reads$V1, "\\."), function(x) x[4]), allfiles$index),]$totalnum_reads <- totalnum_reads$V2

#droplet_num
allfiles$droplet_num <- NA_integer_
allfiles[match(sapply(strsplit(droplet_num$V1, "\\."), function(x) x[4]), allfiles$index),]$droplet_num <- droplet_num$V2

#passfilt_reads
allfiles$totalnum_passfilt_reads <- NA_integer_
allfiles[match(sapply(strsplit(passfilt_reads$V1, "\\."), function(x) x[4]), allfiles$index),]$totalnum_passfilt_reads <- passfilt_reads$V2

#valid_snps
allfiles$totalnum_markers <- NA_integer_
allfiles[match(sapply(strsplit(valid_snps$V1, "\\."), function(x) x[4]), allfiles$index),]$totalnum_markers <- valid_snps$V2


cnames_wanted <- c("index", "totalnum_reads", "droplet_num", "totalnum_passfilt_reads", "totalnum_markers")
bamQC <- cbind(bamQC[match(allfiles$V1, bamQC$pathToBamFile),], allfiles[,match(cnames_wanted, colnames(allfiles))])
bamQC$refgenome <- "hg19"

bamQC$AMB <- NA
bamQC$DBL <- NA
bamQC$SNG <- NA
bamQC$droplet_sng <- NA

irodsPath <- "/projects/NeuroSeq/scRNAseqData/fastq/"
system(paste0("ls ",irodsPath,"run*/*/*/cellranger-hg19/outs/output.demuxlet.doublet0.5.best > ",pathTOfile,"QC/demuxlet/files_with_demuxlet_output.txt"))

##iterate over file in files_with_demuxlet_output.txt
files_with_demuxlet_output <- read.table(paste0(pathTOfile, "QC/demuxlet/files_with_demuxlet_output.txt"))
bamQCshort <- gsub("possorted_genome_bam.bam","",bamQC$pathToBamFile)
demuxletshort <- gsub("output.demuxlet.doublet0.5.best","",files_with_demuxlet_output$V1)
files_with_demuxlet_output <- files_with_demuxlet_output[match(bamQCshort, demuxletshort),]


i=1
for (file in as.character(files_with_demuxlet_output)){
  print(i)
  tmptable <- read.table(file, header=TRUE)
  cells <- gsub("-.+", "", as.character(tmptable$BEST))
  amb <- round(table(cells)["AMB"]*100/sum(table(cells)),2)
  dbl <- round(table(cells)["DBL"]*100/sum(table(cells)),2)
  sng <- round(table(cells)["SNG"]*100/sum(table(cells)),2)

  bamQC[i,]$AMB <- amb
  bamQC[i,]$DBL <- dbl
  bamQC[i,]$SNG <- sng
  bamQC[i,]$droplet_sng <- unname(table(cells)["SNG"])
  
  i=i+1
}


bamQC$altId <- paste0(bamQC$pool, "-", bamQC$bioRep, "-", bamQC$techRep, "-", bamQC$tenXRep, "-", bamQC$timepoint, "-",
                      bamQC$rotenone_treatment)

## perfect replicates (for timepoint, pool, rot treatment, and other replicates)
bamQC$techReplabel <- NA
bamQC$bioReplabel <- NA
bamQC$tenReplabel <- NA

list_techRep <- bamQC[bamQC$techRep==2,]$altId
i=1
for (element in list_techRep){
  alt_techRep <- sapply(strsplit(element, "-"), function(y) {y[3] <- 1; paste(y, collapse="-")}, simplify=TRUE)
  bamQC[bamQC$altId %in% c(element, alt_techRep),]$techReplabel <- as.character(i)
  i=i+1
}

list_bioRep <- bamQC[bamQC$bioRep==2,]$altId
i=1
for (element in list_bioRep){
  alt_bioRep <- sapply(strsplit(element, "-"), function(y) {y[2] <- 1; paste(y, collapse="-")}, simplify=TRUE)
  if (alt_bioRep %in% bamQC$altId){
    bamQC[bamQC$altId %in% c(element, alt_bioRep),]$bioReplabel <- as.character(i)
    i=i+1
  }
}


list_tenXRep <- bamQC[bamQC$tenXRep==2,]$altId
i=1
for (element in list_tenXRep){
  alt_tenXRep <- sapply(strsplit(element, "-"), function(y) {y[4] <- 1; paste(y, collapse="-")}, simplify=TRUE)
  if (alt_tenXRep %in% bamQC$altId){
    bamQC[bamQC$altId %in% c(element, alt_tenXRep),]$tenReplabel <- as.character(i)
    i=i+1
  }
}

#table summary
write.table(bamQC,
            file ="outputTabs/QC_tab.txt", 
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)




####################################
### Cell-line abundance per pool ###
####################################

bamQC <- read.table("outputTabs/QC_tab.txt", header=TRUE)

pathToDemuxlet <- pathTOfile
files <- list.files(pattern=paste0(pathToDemuxlet,"pool*"))
list_pools <- dir(pathToDemuxlet, pattern="*_sample_list.txt")
list_pools <- list_pools[order(as.numeric(sapply(strsplit(list_pools, "_"), function(x) gsub("pool","", x[1]))))]

poolSampleshg19 <- vector(mode="list", length=length(list_pools))
names(poolSampleshg19) <- gsub("_sample.+", "", list_pools)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

for (i in seq_len(length(poolSampleshg19))){
  
  tmppool <- read.table(paste0(pathToDemuxlet, list_pools[i]))
  pool_number <- firstup(gsub("_sample.+", "", list_pools[i]))
  allrunsPool <- gsub("possorted_genome_bam.bam","output.demuxlet.doublet0.5.best",bamQC[as.character(bamQC$pool) %in% pool_number,]$pathToBamFile)
  poolshg19 <- allrunsPool[grep("hg19", allrunsPool)]

  #hg19
  poolSampleshg19[[i]] <- matrix(NA, nrow=length(tmppool$V1), ncol=length(poolshg19))
  rownames(poolSampleshg19[[i]]) <- gsub(".+-", "", as.character(tmppool$V1))
  poolshg19 <- sapply(strsplit(poolshg19, "/"), function(x) paste(c(x[10],x[11],x[12]), collapse="/"))
  poolshg19 <- poolshg19[order(poolshg19)]
  colnames(poolSampleshg19[[i]]) <- poolshg19
  
}


files_with_demuxlet_output <- read.table(paste0(pathToDemuxlet, "files_with_demuxlet_output.txt"))
files_hg19 <- as.character(files_with_demuxlet_output$V1)[grep("hg19", as.character(files_with_demuxlet_output$V1))]

for (element in 1:length(poolSampleshg19)){
  
  tmp <- colnames(poolSampleshg19[[element]])

  for (el in tmp){
    
    tmpfile <- read.table(files_hg19[grep(paste0("/",el,"/"), files_hg19)], header=TRUE)
    tmpfile <- tmpfile[grep("SNG",tmpfile$BEST),]
    tmpcounter <- table(gsub(".+-", "", gsub("SNG-", "", tmpfile$BEST)))
    
    cnames <- match(el, colnames(poolSampleshg19[[element]]))
    tmpcounter <- tmpcounter[match(names(poolSampleshg19[[element]][,cnames]), names(tmpcounter))]
    tmpcounter[is.na(tmpcounter)] <- 0
    poolSampleshg19[[element]][,cnames] <- unname(tmpcounter)
      
  }

}

saveRDS(poolSampleshg19, "outputTabs/poolSampleshg19.RDS")


