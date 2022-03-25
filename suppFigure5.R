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
library(grid)
library(tidyverse)
library(ggsignif)
library(ggtext)
library(glue)
library(ggupset)
library(wesanderson)
library(tidyverse)

source("functionsToImport.R")

## outlier analysis
summData <- readRDS("analysis/outputTabs/summData_donorCtype_Cells.RDS")
summData <- subset(summData, quantile!=1)

summData <- sapply(unique(summData$tp), function(x){
  
  tmp <- subset(summData, tp==x)
  annotTmp <- sapply(unique(tmp$annot), function(y){
    
    tmp2 <- subset(tmp, annot==y)
    tmp2$zsco <- round((tmp2$cfrac-mean(tmp2$cfrac))/sd(tmp2$cfrac),2)
    tmp2$outlier <- abs(tmp2$zsco)>2
    tmp2
    
  }, simplify=F)
  
  annotTmp <- do.call("rbind", annotTmp)
  rownames(annotTmp) <- NULL
  annotTmp
  
}, simplify=F)

summData <- do.call("rbind", summData)
rownames(summData) <- NULL


summData$pool <- sapply(strsplit(summData$donor_extended, "/"), function(x) x[length(x)])
diffmetrics <- read.table("suppTabs/suppTable1.txt",
                          header=T)
outlier_analysis <- data.frame(donor_extended=unique(summData$donor_extended),
                               outlier_tp11_z2=NA,
                               outlier_tp30_z2=NA,
                               outlier_tp52_z2=NA,
                               outlier_alltp_z2=NA)
avail <- sapply(unique(summData$donor_extended), function(x){ paste(unique(subset(summData, donor_extended==x)$tp), collapse="-")}, simplify=T)
outlier_analysis$avail <- unname(avail[match(outlier_analysis$donor_extended, names(avail))])



#all outliers (all tp)
outliers_z2 <- unique(summData[summData$outlier==TRUE,]$donor_extended)
outlier_analysis$outlier_alltp_z2 <- outlier_analysis$donor_extended %in% outliers_z2
outlier_analysis$donor_id <- sapply(strsplit(outlier_analysis$donor_extended,"/"), function(x) paste(x[-length(x)],collapse="/"))
outlier_analysis$allMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$all
outlier_analysis$delMutBurden <- diffmetrics[match(outlier_analysis$donor_id, diffmetrics$donor_iPSC),]$deleterious


tpOutlierSet <- function(summData, timePoint="D11", colname="outlier_tp11_z2"){
  outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE)$donor_extended)
  colIndex <- match(colname, colnames(outlier_analysis))
  outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
  return(outlier_analysis)
}

outlier_analysis <- tpOutlierSet(summData, timePoint="D11", colname="outlier_tp11_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D30", colname="outlier_tp30_z2")
outlier_analysis <- tpOutlierSet(summData, timePoint="D52", colname="outlier_tp52_z2")

tpOutlierSetPercentile <- function(summData, timePoint="D11", colname="outlier_tp11_up5"){
  
  sense <- gsub("[0-9]","",sapply(strsplit(colname,"_"), function(x) x[3]))
  
  if (sense=="up"){
    
    outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE & top5==TRUE)$donor_extended)
    colIndex <- match(colname, colnames(outlier_analysis))
    outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
    
  } else if (sense=="down"){
    
    outliers_tp <- unique(subset(summData, tp==timePoint & outlier==TRUE & bot5==TRUE)$donor_extended)
    colIndex <- match(colname, colnames(outlier_analysis))
    outlier_analysis[,colIndex] <- outlier_analysis$donor_extended %in% outliers_tp
    
  }
  
  return(outlier_analysis)
}


outlierLong <- outlier_analysis[,c("donor_extended","outlier_alltp_z2","outlier_tp11_z2","outlier_tp30_z2","outlier_tp52_z2","allMutBurden","delMutBurden")]
outlierLong <- as.data.frame(outlierLong %>% pivot_longer(-c(1:5), names_to="MutClass", values_to="Burden"))
outlierLong <- subset(outlierLong, !is.na(Burden))

facet_labeller = labeller(
  MutClass = c("allMutBurden" = "All variants",
               "delMutBurden" = "Deleterious"))
outlierLong$MutClass <- factor(outlierLong$MutClass, levels=unique(outlierLong$MutClass))
outlierLongest <- as.data.frame(outlierLong %>% pivot_longer(-c(1,6,7), names_to="outlierClass", values_to="Group"))

vec <- c("Agg All TP", "Day 11", "Day 30", "Day52")
names(vec) <-  unique(outlierLongest$outlierClass)
outlierLongest$outlierClass <- unname(vec[outlierLongest$outlierClass])

vec4 <- c("Outlier","Non-outlier")
names(vec4) <- c("TRUE","FALSE")
outlierLongest$Group <- unname(vec4[as.character(outlierLongest$Group)])
my_comparisons <- list(c("Outlier","Non-outlier"))

myplot <- ggplot(outlierLongest, aes(x=Group, y=Burden))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.25, alpha=0.25, aes(col=Group), size=3)+
  facet_grid(MutClass~outlierClass, scales="free_y", labeller = facet_labeller)+
  xlab("")+ylab("")+
  theme_bw()+
  #ylab("Mutational burden")+
  theme(strip.text =element_text(size=20),
        axis.text.y=element_text(size=20,),
        axis.text.x=element_text(size=20, angle=90,hjust=1, vjust=0.5),
        axis.title=element_text(size=13))+
  geom_signif(test="wilcox.test",
              comparisons=my_comparisons,
              map_signif_level = T,
              textsize=4.5,
              col="black")+
  scale_colour_manual(values=c("black","red"))+
  guides(col=FALSE)

facet_bounds <- read.table(header=TRUE,
                           text=                           
                             "MutClass ymin  ymax  breaks
allMutBurden  0 300 50
delMutBurden 0 100 20",
                           stringsAsFactors=FALSE)

ff <- with(facet_bounds,
           data.frame(Burden=c(ymin,ymax),
                      MutClass=c(MutClass,MutClass)))

myplot <- myplot + geom_point(data=ff, x=NA)

q <- ggplot_build(myplot)

mask_panel <- q$data[[3]]$PANEL==1

q$data[[3]]$y[mask_panel] <- q$data[[3]]$y[mask_panel]-40
q$data[[3]]$yend[mask_panel] <- q$data[[3]]$yend[mask_panel]-40

q$data[[3]]$y[!mask_panel] <- q$data[[3]]$y[!mask_panel]-20
q$data[[3]]$yend[!mask_panel] <- q$data[[3]]$yend[!mask_panel]-20


pv <- generics::tidy(with(outlierLongest, pairwise.wilcox.test(Burden, interaction(Group,outlierClass,MutClass), p.adjust.method = "BH")))

pv_final <- as.data.frame(pv %>% 
                            separate(group1, c("g1_1", "g1_2", "g1_3"), sep="\\.") %>% 
                            separate(group2, c("g2_1", "g2_2", "g2_3"), sep="\\.") %>% 
                            filter(g1_1 != g2_1) %>% filter(g1_2 == g2_2) %>% filter(g1_3 == g2_3) %>%
                            mutate(p=ifelse(p.value > 0.05, "", ifelse(p.value > 0.01,"*", "**"))))


# each pvalue is repeated three times in this dataset.  s
q$data[[3]]$annotation <- rep(pv_final$p, each=3)
# remove non significants
q$data[[3]] <- q$data[[3]][q$data[[3]]$annotation != "",]
# and the final plot
figure <- ggplot_gtable(q)

pdf(file=paste0("figures/suppFigs/suppfig5.pdf"))
plot(figure)
dev.off()

