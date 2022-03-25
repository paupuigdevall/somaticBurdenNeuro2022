
library(ggplot2)
library(gridExtra)
library(scales)
library(ggpubr)
library(viridis)
library(wesanderson)

source("functionsToImport.R")

##############
### Fig 2A ###
##############


poolSampleshg19_perc <- readRDS("analysis/outputTabs/poolSampleshg19.RDS")

percentageCompute <- function(poolSampleshg19_perc){
  
  for (i in 1:length(poolSampleshg19_perc)){
    
    tmp_cnum <- dim(poolSampleshg19_perc[[i]])[2]
    tmp_cnum <- 1:tmp_cnum
    
    for (el in tmp_cnum){
      poolSampleshg19_perc[[i]][,el] <- round(poolSampleshg19_perc[[i]][,el]*100/sum(poolSampleshg19_perc[[i]][,el]),2)
    }
    
  }
  
  return(poolSampleshg19_perc)
  
}

poolSampleshg19_perc <- percentageCompute(poolSampleshg19_perc)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


test <- poolSampleshg19_perc[["pool2"]]
test <- as.data.frame(test)
test$clines <- rownames(test)
rownames(test) <- NULL
colnames(test)[-length(colnames(test))] <- gsub(".*[/]([^.]+)[/].*", "\\1",  colnames(test)[-length(colnames(test))])

test <- test %>% gather(key=pool, value=percentage, -clines )
test$pool <- factor(test$pool, levels=levels(as.factor(test$pool))[order(gsub(".+-", "", levels(as.factor(test$pool))))])

tmp_plot <- ggplot(test, aes(x=pool, y=clines, fill=percentage))+ 
  geom_tile(col="black") +
  scale_fill_viridis(name="Percentage",
                     guide = guide_colourbar(barwidth = 15),
                     discrete = FALSE, limits=c(0,100), direction=-1) +
  #ggtitle("Pool2")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold"),
        axis.title=element_blank(),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=20),
        legend.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=18),
        legend.position = "top")


pdf(file=paste0("figures/suppFigs/suppfig2A.pdf"))
plot(tmp_plot)
dev.off()


##############
### Fig 2B ###
##############

QC_tab <- read.table("analysis/outputTabs/QC_tab.txt", header=TRUE)
QC_tab$pool_mod <- gsub("Pool", "", QC_tab$pool)
QC_tab$pool_mod <- factor(QC_tab$pool_mod, levels=as.character(c(1:17,20,21)))
QC_tab$techReplabel <- as.character(QC_tab$techReplabel)

pos <- position_jitter(width = 0.15, seed = 1)

consideredTp <- subset(QC_tab, timepoint!="D119")

tmp_plot <- ggplot(data=subset(consideredTp, timepoint=="D52"), aes(x=pool_mod, y=SNG, label=altId))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  #facet_wrap(~timepoint, ncol =3)+
  geom_point(data=subset(consideredTp, SNG<=51), position=pos, size=4, alpha=0.5, col="red")+
  geom_point(data=subset(consideredTp, SNG>51), position=pos, size=4, alpha=0.5, col="black")+
  ylab("Singletons (%)")+
  #ggtitle("Percentage of singletons demultiplexed per pool run")+
  theme_bw()+
  theme(axis.text.x= element_text(size=20, angle=90, vjust=0.5, hjust=0.5),
        axis.text.y= element_text(size=20))+
  #geom_label_repel(data=subset(QC_tab, QC_tab$SNG<70), aes(label=altId), position=pos, force=10)+
  geom_hline(yintercept=51, linetype="dashed", alpha=0.3)+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(breaks=seq(0,100,10))+
  xlab("Pool")+
  xlab("")+
  ylab("")

pdf(file=paste0("figures/suppFigs/suppfig2B.pdf"))
plot(tmp_plot)
dev.off()

##############
### Fig 2C ###
##############

#input number of droplets

tmp_plot <- ggplot(data=consideredTp, aes(x=pool_mod, y=droplet_num))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col=timepoint),alpha=0.5, width=0.2, size=4)+
  ylab("")+xlab("")+
  theme_bw()+
  theme(axis.text.x= element_text(size=20, angle=90, vjust=0.5, hjust=0.5),
        axis.text.y= element_text(size=20),
        legend.position="top",
        legend.title=element_text(face="bold", size=22),
        legend.text=element_text(size=18))+
  scale_y_continuous(labels=comma)+
  scale_color_manual(name="Timepoint",
                     values=wes_palette(n = length(unique(consideredTp$timepoint)), "GrandBudapest1"))+
  xlab("Pool")+
  ylab("Number of cells")


pdf(file=paste0("figures/suppFigs/suppfig2C.pdf"))
plot(tmp_plot)
dev.off()


##############
### Fig 2J ###
##############

metadata <- readRDS("analysis/outputTabs/suppData1.RDS")
heatmapQC <- sapply(unique(metadata$batch), function(x){
  
  tmp <- subset(metadata, batch==x)
  cfractmp <- table(tmp$annot)/sum(table(tmp$annot))
  newdf <- data.frame(annot=unique(metadata$annot),
                      cfrac=NA)
  newdf$cfrac <- as.numeric(cfractmp)[match(newdf$annot, names(cfractmp))]
  newdf$batch <- x
  newdf$timepoint <- unique(tmp$tp)
  newdf$comb <- paste0(newdf$timepoint,"-",newdf$batch)
  newdf
}, simplify=F)

heatmapQC <- do.call("rbind", heatmapQC)
rownames(heatmapQC) <- NULL
heatmapQC[is.na(heatmapQC$cfrac),]$cfrac <- 0
heatmapQC <- heatmapQC[order(heatmapQC$timepoint),]
heatmapQC$comb <- factor(heatmapQC$comb, levels=rev(unique(heatmapQC$comb)))
annotLevels <- c("FPP-1","FPP-2","FPP-3","proFPP-1","proFPP-2",
                 "Astro","DA","Epend-1","Sert-like","Pro.Sert-like",
                 "Unk-1","Unk-2")
heatmapQC$annot <- factor(heatmapQC$annot, levels=annotLevels)


tmp_plot <- ggplot(heatmapQC, aes(x=annot, y=comb, fill=cfrac))+ 
  geom_tile(col="black") +
  scale_fill_viridis(name="",
                     guide = guide_colourbar(barwidth = 20),
                     discrete = FALSE, limits=c(0,1), direction=-1) +
  ggtitle("")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=11),
        axis.title=element_blank(),
        axis.text.x= element_text(angle = 90, hjust = 1, size=20),
        axis.text.y= element_blank(),
        legend.title=element_text(size=22, face="bold"),
        axis.ticks.y=element_blank(),
        legend.position="top",
        legend.text=element_text(size=20))

pdf(file=paste0("figures/suppFigs/suppfig2J.pdf"))
plot(tmp_plot)
dev.off()


