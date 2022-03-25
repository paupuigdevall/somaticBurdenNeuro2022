replicateLines <- function(mutationsGeneSample, dopaminergic_neurons, geneUniverse){
  
  colnames(mutationsGeneSample) <- gsub("\\.","-", colnames(mutationsGeneSample))
  
  ## Filter-out those without gene mutational burden
  metadata_table <- read.table("tabs/hipsci_metadata/hipsci.qc1_sample_info.20170927.tsv",
                               header=TRUE, sep="\t")
  replicates <- split(metadata_table$name, metadata_table$donor)
  replicates <- sapply(replicates, function(x) as.character(x[x %in% colnames(mutationsGeneSample)]))
  replicates <- replicates[elementNROWS(replicates)>1]
  #stopifnot(all(elementNROWS(replicates)==2))
  
  ## Filter-out those cell-lines without differentiation outcome and which are no replicates
  replicates <- replicates[sapply(replicates, function(x) length(x[x %in% dopaminergic_neurons$donor_id])>1)]
  #stopifnot(all(elementNROWS(replicates)==2))
  
  dop_vector <- dopaminergic_neurons$diff_efficiency_modelScore
  names(dop_vector) <- dopaminergic_neurons$donor_id
  outcomeRep <- sapply(replicates, function(x){
    tmp <- dop_vector[x]
    tmp <- tmp[!is.na(tmp)]
  }, simplify=FALSE) 
  
  outcomeRep_discordant <- outcomeRep[sapply(outcomeRep, function(x) length(unique(x))==2)]
  outcomeRep_discordant <- sapply(outcomeRep_discordant, function(x) sort(x), simplify=FALSE)
  
  outcomeRep_concordant <- outcomeRep[sapply(outcomeRep, function(x) length(unique(x))==1)]
  outcomeRep_succConcordant <- outcomeRep_concordant[sapply(outcomeRep_concordant, function(x) all(x %in% "succeeded"))]
  outcomeRep_failedConcordant <- outcomeRep_concordant[sapply(outcomeRep_concordant, function(x) all(x %in% "failed"))]
  
  to_add_list <- geneUniverse[is.na(match(geneUniverse$symbol, rownames(mutationsGeneSample))),]$symbol
  to_add_mtrix <- matrix(0, nrow=length(to_add_list), ncol = length(colnames(mutationsGeneSample)))
  colnames(to_add_mtrix) <- colnames(mutationsGeneSample)
  rownames(to_add_mtrix) <- to_add_list
  mutationsGeneSample <- rbind(mutationsGeneSample, to_add_mtrix)
  
  ## Discordant replicates: Genes exclusively mutated in failed
  failed_mut <- sapply(outcomeRep_discordant, function(x) {
    
    colIndexs <- match(names(x), colnames(mutationsGeneSample))
    tmp <- mutationsGeneSample[,colIndexs]
    rownames(tmp[tmp[,1]>0 & tmp[,2]==0,])
    
  }, simplify=FALSE)
  failed_mut <- unique(unname(unlist(failed_mut)))
  
  ## Discordant replicates: Genes exclusively mutated in succeeded
  succ_mut <- sapply(outcomeRep_discordant, function(x) {
    
    colIndexs <- match(names(x), colnames(mutationsGeneSample))
    tmp <- mutationsGeneSample[,colIndexs]
    rownames(tmp[tmp[,1]==0 & tmp[,2]>0,])
    
  }, simplify=FALSE)
  succ_mut <- unique(unname(unlist(succ_mut)))
  
  ## Failed concordant replicates: Both genes mutated
  failed_mut_conc <- sapply(outcomeRep_failedConcordant, function(x) {
    
    colIndexs <- match(names(x), colnames(mutationsGeneSample))
    tmp <- mutationsGeneSample[,colIndexs]
    rownames(tmp[tmp[,1]>0 & tmp[,2]>0,])
    
  }, simplify=FALSE)
  failed_mut_conc <- unique(unname(unlist(failed_mut_conc)))
  
  
  ## Succeeded concordant replicates: Both genes mutated
  succ_mut_conc <- sapply(outcomeRep_succConcordant, function(x) {
    
    colIndexs <- match(names(x), colnames(mutationsGeneSample))
    tmp <- mutationsGeneSample[,colIndexs]
    rownames(tmp[tmp[,1]>0 & tmp[,2]>0,])
    
  }, simplify=FALSE)
  succ_mut_conc <- unique(unname(unlist(succ_mut_conc)))
  
  
  only_failed <- unique(c(failed_mut, failed_mut_conc))
  only_succ <- unique(c(succ_mut, succ_mut_conc))
  
  only_failed <- setdiff(only_failed, only_succ)
  only_succ <- setdiff(only_succ, only_failed)
  
  only_failed <- geneUniverse[match(only_failed, geneUniverse$symbol),]$entrezid
  only_failed <- unique(only_failed[!is.na(only_failed)])
  
  only_succ <- geneUniverse[match(only_succ, geneUniverse$symbol),]$entrezid
  only_succ <- unique(only_succ[!is.na(only_succ)])
  
  outcomeRepAnalysis <- list("only_failed"= only_failed,
                             "only_successful"=only_succ)
  
  return(outcomeRepAnalysis)
  
}







GOenrichmentAndReport <- function(geneIds, universeGeneIds,
                                  minSize=3, maxSize=300, minCount=3,
                                  minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                  label="allDE"){
  
  
  GOparams <- new("GOHyperGParams", geneIds=geneIds,
                  universeGeneIds=universeGeneIds,
                  annotation="org.Hs.eg.db",
                  ontology="BP", pvalueCutoff=0.05, conditional=TRUE,
                  minSizeCutoff=3, maxSizeCutoff=300, orCutoff=1.5,
                  testDirection="over")
  
  hgOverGOBP <- hyperGTest(GOparams)
  
  report <- data.frame(GOBPID=as.vector(names(geneIdUniverse(hgOverGOBP))),
                       Pvalue=pvalues(hgOverGOBP),
                       OddsRatio=oddsRatios(hgOverGOBP),
                       ExpCount=expectedCounts(hgOverGOBP),
                       Count=geneCounts(hgOverGOBP),
                       Size=universeCounts(hgOverGOBP),
                       stringsAsFactors=FALSE)
  
  ## discard gene sets that do not meet a minimum and maximum number of genes
  report <- report[report$Size >= minSize & report$Size <= maxSize, , drop=FALSE]
  
  ## discard gene sets that show a p.value>0.05
  report <- report[report$Pvalue < p.value, , drop=FALSE]
  
  ## discard gene sets that do not satisfy the OR cutoff
  report <- report[report$OddsRatio >= minOddsRatio & report$Count >= minCount, , drop=FALSE]
  
  ## apply the maximum-number-of-GO-terms-reported cutoff and sort by odds ratio
  maxReported <- min(nrow(report))
  report <- report[sort(report$OddsRatio, decreasing=TRUE, index.return=TRUE)$ix[1:maxReported], ]
  
  if (dim(report[complete.cases(report),])[1]==0){
    message <- "No GO terms enriched"
    print(message)
    return(message)
  }
  
  ## add the symbol and GO term description, information
  reportGenes <- geneIdsByCategory(hgOverGOBP, report$GOBPID)
  reportGeneSyms <- lapply(reportGenes, annotate::getSYMBOL, "org.Hs.eg.db")
  highlightGeneMask <- lapply(reportGenes, function(x, hgenes, fmt) x %in% hgenes, highlightGenes)
  reportGeneSyms <- mapply(function(genes, mask, fmt) ifelse(mask, sprintf(fmt, genes), genes),
                           reportGeneSyms, highlightGeneMask, MoreArgs=list(fmt=highlightStr), SIMPLIFY=FALSE)
  reportGeneSyms <- sapply(reportGeneSyms, paste, collapse=", ")
  reportGenes <- sapply(reportGenes, paste, collapse=", ")
  report <- cbind(report,
                  Term=sapply(mget(report$GOBPID, GOstats:::GOenv("TERM")), Term),
                  GeneSyms=reportGeneSyms,
                  label=label)
  rownames(report) <- NULL
  
  return(report)
  
}







highlightTopOrAndNeuroGO <- function(enrichment_failed_all){
  
  enrichment_failed_all$GeneSyms <- NULL
  enrichment_failed_all$position <- c(1:dim(enrichment_failed_all)[1])
  enrichment_failed_all$highlight <- FALSE
  enrichment_failed_all$highlight[1:5] <- TRUE
  
  vecMatch <- c("axon","neuron","glial","brain",
                "hindbrain","forebrain","midbrain","synapse","chromatin",
                "cerebellum", "neural","cortex","neurogenesis",
                "axonogenesis","nervous","hippocampus","neurotransmitter",
                "dopaminergic", "axenome", "action potential","synaptic")
  
  mask <- rowSums(sapply(vecMatch, function(x){
    grepl(x,enrichment_failed_all$Term)
  }))>0
  
  enrichment_failed_all[mask,]$highlight <- TRUE
  enrichment_failed_all <- subset(enrichment_failed_all, highlight==T)
  enrichment_failed_all <- subset(enrichment_failed_all, position<=25)
  
  return(enrichment_failed_all)
}




fillNA_fig2 <- function(enrichmentGOexclusive, highlight=F){
  alltabs <- sapply(unique(enrichmentGOexclusive$Term), function(x){
    tabs <- sapply(unique(enrichmentGOexclusive$label), function(y){
      dimensions <- dim(subset(enrichmentGOexclusive, Term==x & label==y))[1]
      if (dimensions==0){
        df <- data.frame(GOBPID=NA,
                         Pvalue="ns",
                         OddsRatio=NA,
                         ExpCount=NA,
                         Count=NA,
                         Size=NA,
                         Term=x,
                         label=y,
                         position=NA,
                         highlight=F)
      } else {
        df <- NA
      }
      return(df)
    }, simplify=F)
    
    tabs <- tabs[!sapply(tabs, function(x) all(is.na(x)))]
    tabs <- do.call("rbind", tabs)
    rownames(tabs) <- NULL
    tabs
  }, simplify=F)
  
  alltabs <- do.call("rbind", alltabs)
  rownames(alltabs) <- NULL
  return(alltabs)
} 




buildOneSidedAllDelTab <- function(oneSidedTests, subClass="DDD_allMutBurden_greater"){
  
  tests <- oneSidedTests[subClass][[1]]
  subTests <- listRDS[tests]
  len <- length(subTests)
  
  df <- sapply(1:len, function(y){
    geneSet <- strsplit(names(subTests)[y],"_")[[1]][1]
    mutClass <- strsplit(names(subTests)[y],"_")[[1]][2]
    alternative <- strsplit(names(subTests)[y],"_")[[1]][3]
    if (is.null(dim(subTests[[y]]))){
      df <- subTests[[y]]
      df <- as.data.frame(df)
      colnames(df) <- names(subTests)[y]
    } else {
      df <- subTests[[y]]
      df <- as.data.frame(df)
      colnames(df) <- paste0(names(subTests)[y],"_",1:1000)
    }
    df
  },simplify=F)
  
  df <- do.call("cbind", df)
  df$perm <- paste0("perm",1:dim(df)[1])
  index <- grep("perm",colnames(df))
  df <- as.data.frame(df %>% pivot_longer(-c(index), names_to="test", values_to="pval"))
  df$perm <- NULL
  df$test <- gsub("_[0-9]+$","",df$test)
  df$geneSet <- sapply(strsplit(df$test,"_"), function(x) x[1])
  df$mutClass <- sapply(strsplit(df$test,"_"), function(x) x[2])
  df$alternative <- sapply(strsplit(df$test,"_"), function(x) x[3])
  
  return(df)
  
}





summgeneSetGenerator <- function(allSnvsDinuc, geneSet="DDD", out="Failed"){
  
  allSnvsDinuc_outcome <- subset(allSnvsDinuc, outcome==out)
  
  if (geneSet=="DDD"){
    
    summgeneSetDD <- as.data.frame(table(allSnvsDinuc_outcome$SampleconsensusSimplified, allSnvsDinuc_outcome$varCategoryFiltered, allSnvsDinuc_outcome$DDgenes))
    colnames(summgeneSetDD) <- c("chromState","varClass","geneSetDD","Freq")
    setChanger <- c("DD","non-DD")
    names(setChanger) <- c("TRUE","FALSE")
    summgeneSetDD$geneSetDD <- unname(setChanger[as.character(summgeneSetDD$geneSetDD)])
    summgeneSetDD$fill <- paste0(summgeneSetDD$chromState,"-", summgeneSetDD$geneSetDD)
    summgeneSetDD$fill <- factor(summgeneSetDD$fill, c("active-DD","active-non-DD",
                                                       "inactive-DD","inactive-non-DD",
                                                       "bivalent-DD","bivalent-non-DD"))
    summgeneSetDD$freqVarClass <- NA
    summgeneSetDD$varClass <- factor(summgeneSetDD$varClass, levels=c("synonymous","other","missense", "missPatho", "ptv"))
    summgeneSetDD$chromState <- factor(summgeneSetDD$chromState, levels=c("active","inactive","bivalent"))
    
    results_summgeneSetDD <- sapply(unique(summgeneSetDD$varClass), function(x){
      summgeneSetDD[summgeneSetDD$varClass==x,]$freqVarClass <- round(summgeneSetDD[summgeneSetDD$varClass==x,]$Freq/sum(summgeneSetDD[summgeneSetDD$varClass==x,]$Freq),4)
      summgeneSetDD[summgeneSetDD$varClass==x,]
    }, simplify=FALSE )
    results_summgeneSetDD <- do.call("rbind", results_summgeneSetDD)
    #results_summgeneSetDD <- subset(results_summgeneSetDD, varClass!="Other")
    
    return(results_summgeneSetDD)
    
  } else if (geneSet=="DDDdominant"){
    
    summgeneSetDDdom <- as.data.frame(table(allSnvsDinuc_outcome$SampleconsensusSimplified, allSnvsDinuc_outcome$varCategoryFiltered, allSnvsDinuc_outcome$DDdominantgenes))
    colnames(summgeneSetDDdom) <- c("chromState","varClass","geneSetDDdom","Freq")
    setChanger <- c("DDdominant","non-DDdominant")
    names(setChanger) <- c("TRUE","FALSE")
    summgeneSetDDdom$geneSetDDdom <- unname(setChanger[as.character(summgeneSetDDdom$geneSetDD)])
    summgeneSetDDdom$fill <- paste0(summgeneSetDDdom$chromState,"-", summgeneSetDDdom$geneSetDD)
    summgeneSetDDdom$fill <- factor(summgeneSetDDdom$fill, c("active-DDdominant","active-non-DDdominant",
                                                             "inactive-DDdominant","inactive-non-DDdominant",
                                                             "bivalent-DDdominant","bivalent-non-DDdominant"))
    summgeneSetDDdom$freqVarClass <- NA
    summgeneSetDDdom$varClass <- factor(summgeneSetDDdom$varClass, levels=c("synonymous","other","missense", "missPatho", "ptv"))
    summgeneSetDDdom$chromState <- factor(summgeneSetDDdom$chromState, levels=c("active","inactive","bivalent"))
    
    results_summgeneSetDDdom <- sapply(unique(summgeneSetDDdom$varClass), function(x){
      summgeneSetDDdom[summgeneSetDDdom$varClass==x,]$freqVarClass <- round(summgeneSetDDdom[summgeneSetDDdom$varClass==x,]$Freq/sum(summgeneSetDDdom[summgeneSetDDdom$varClass==x,]$Freq),4)
      summgeneSetDDdom[summgeneSetDDdom$varClass==x,]
    }, simplify=FALSE )
    results_summgeneSetDDdom <- do.call("rbind", results_summgeneSetDDdom)
    #results_summgeneSetDD <- subset(results_summgeneSetDD, varClass!="Other")
    
    return(results_summgeneSetDDdom)
    
    
  } else if (geneSet=="Cosmic"){
    
    summgeneSetCosmic <- as.data.frame(table(allSnvsDinuc_outcome$SampleconsensusSimplified, allSnvsDinuc_outcome$varCategoryFiltered, allSnvsDinuc_outcome$Cosmicgenes))
    colnames(summgeneSetCosmic) <- c("chromState","varClass","geneSetCosmic","Freq")
    setChanger <- c("Cosmic","non-Cosmic")
    names(setChanger) <- c("TRUE","FALSE")
    summgeneSetCosmic$geneSetCosmic <- unname(setChanger[as.character(summgeneSetCosmic$geneSetCosmic)])
    summgeneSetCosmic$fill <- paste0(summgeneSetCosmic$chromState,"-", summgeneSetCosmic$geneSetCosmic)
    summgeneSetCosmic$fill <- factor(summgeneSetCosmic$fill, c("active-Cosmic","active-non-Cosmic",
                                                               "inactive-Cosmic","inactive-non-Cosmic",
                                                               "bivalent-Cosmic","bivalent-non-Cosmic"))
    summgeneSetCosmic$freqVarClass <- NA
    summgeneSetCosmic$varClass <- factor(summgeneSetCosmic$varClass, levels=c("synonymous","other","missense", "missPatho", "ptv"))
    summgeneSetCosmic$chromState <- factor(summgeneSetCosmic$chromState, levels=c("active","inactive","bivalent"))
    
    results_summgeneSetCosmic <- sapply(unique(summgeneSetCosmic$varClass), function(x){
      summgeneSetCosmic[summgeneSetCosmic$varClass==x,]$freqVarClass <- round(summgeneSetCosmic[summgeneSetCosmic$varClass==x,]$Freq/sum(summgeneSetCosmic[summgeneSetCosmic$varClass==x,]$Freq),4)
      summgeneSetCosmic[summgeneSetCosmic$varClass==x,]
    }, simplify=FALSE )
    results_summgeneSetCosmic <- do.call("rbind", results_summgeneSetCosmic)
    #results_summgeneSetCosmic <- subset(results_summgeneSetCosmic, varClass!="Other")
    
    
    return(results_summgeneSetCosmic)
    
  }
  
}


varNum <- function(results_summgeneSetDD_f){
  varClassesDD <- sapply(unique(results_summgeneSetDD_f$varClass), function(x) sum(subset(results_summgeneSetDD_f, varClass==x)$Freq), simplify=TRUE)
  names(varClassesDD) <-  as.character(unique(results_summgeneSetDD_f$varClass))
  return(varClassesDD)
}

nullGenerator <- function(allSnvsDinuc, numSynVar=400, select_col="DDgenes", out="Failed"){
  
  if (out=="Failed"){
    allsynVar <- subset(allSnvsDinuc, varCategoryFiltered=="synonymous" & outcome=="Failed")
    
  } else if (out=="Successful"){
    allsynVar <- subset(allSnvsDinuc, varCategoryFiltered=="synonymous" & outcome=="Successful")
  }
  
  set.seed(123)
  perm <- sapply(1:10000, function(x) {
    sample(1:dim(allsynVar)[1], numSynVar)
  }, simplify=FALSE)
  
  nullPerm <- sapply(perm, function(x){
    match_col <- match(select_col, colnames(allsynVar))
    setgenes <- gsub("genes","", select_col)
    tt <- as.data.frame(table(allsynVar[x,]$SampleconsensusSimplified, allsynVar[x,match_col]))
    colnames(tt) <- c("chromState", "geneSet","Freq")
    vec <- c(setgenes,paste0("non-", setgenes))
    names(vec) <- c("TRUE","FALSE")
    tt$geneSet <- unname(vec[as.character(tt$geneSet)])
    tt$freqVarClass <- round(tt$Freq/sum(tt$Freq),4)
    tt$fill <- paste0(tt$chromState,"-",tt$geneSet)
    tt
    
  } , simplify=FALSE)
  
  return(nullPerm)
  
}

pval_tests <- function(summgeneSetDD, lof_null, variantClass="missPatho"){
  
  nperm=length(lof_null)
  matched_order <- sapply(lof_null, function(x) x$fill)[,1]
  tmp_summgeneSetDD <- subset(summgeneSetDD, varClass==variantClass)
  tmp_summgeneSetDD <- tmp_summgeneSetDD[match(tmp_summgeneSetDD$fill, matched_order),]
  
  ## comparative using frequencies, not counts due to: 
  ## n_synVar<n_missPatho
  enrichmentN <- sapply(lof_null, function(x) {
    stopifnot(sum(x$Freq)!=sum(tmp_summgeneSetDD$Freq))
    x$ratio <- x$freqVarClass/tmp_summgeneSetDD$freqVarClass
    x$ratio>1
  }, simplify=FALSE)
  enrichmentN <- do.call("rbind", enrichmentN)
  
  enrichmentN <- colSums(enrichmentN)/nperm
  names(enrichmentN) <- matched_order
  
  return(enrichmentN)
  
}




matrix_pval_generator <- function(ddd_vec, geneSet="DDD"){
  
  pval_ddd <- as.data.frame(sapply(ddd_vec, function(x) get(x)))
  colnames(pval_ddd) <- sapply(strsplit(colnames(pval_ddd),"_"), function(x) x[3])
  pval_ddd$chromState <- rownames(pval_ddd)
  rownames(pval_ddd) <- NULL
  mtch_col <- match("chromState", colnames(pval_ddd))
  pval_ddd <- pval_ddd[,c(mtch_col,which(!colnames(pval_ddd) %in% "chromState"))]
  
  pval_ddd <-gather(pval_ddd, key=outcome, value=pval, -c(1))
  vec <- c("Failed","Successful")
  names(vec) <- unique(pval_ddd$outcome)
  
  pval_ddd$outcome <- unname(vec[pval_ddd$outcome])
  pval_ddd <- pval_ddd[!grepl("bivalent",pval_ddd$chromState),]
  pval_ddd$geneSetgroup <- sapply(strsplit(pval_ddd$chromState,"-"), function(x) paste0(x[2:length(x)], collapse=""))
  pval_ddd$chromState2 <- sapply(strsplit(pval_ddd$chromState,"-"), function(x) x[1])
  return(pval_ddd)
}



pvalConverter <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.05)){
    symbolVec[vectorPval<0.05] <- "*"
  } 
  
  if (any(vectorPval<0.01)){
    symbolVec[vectorPval<0.01] <- "**"
  }
  
  if (any(vectorPval<0.001)){
    symbolVec[vectorPval<0.001] <- "***"
  }
  return(symbolVec)
}



firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


enrichmentCalc <- function(summResults, geneUniverse, geneSet="ddd"){
  
  annotTp <- sapply(unique(summResults$timepoint), function(x){
    
    tmp1 <- subset(summResults, timepoint==x)
    
    annot <- sapply(unique(tmp1$annot), function(y){
      
      tmp <- subset(tmp1, annot==y)
      stopifnot(all(tmp$geneDE %in% geneUniverse$symbol))
      mtrix <- matrix(NA, 2, 2)
      colnames(mtrix) <- c("DE","non-DE")
      rownames(mtrix) <- c("In-GeneSet","Out-GeneSet")
      mtrix[1,1] <- sum(geneUniverse[match(tmp$geneDE, geneUniverse$symbol),match(geneSet,colnames(geneUniverse))])
      mtrix[1,2] <- sum(geneUniverse[!(geneUniverse$symbol %in% tmp$geneDE), match(geneSet,colnames(geneUniverse))])
      mtrix[2,1] <- sum(!geneUniverse[match(tmp$geneDE, geneUniverse$symbol),match(geneSet,colnames(geneUniverse))])
      mtrix[2,2] <- sum(!geneUniverse[!(geneUniverse$symbol %in% tmp$geneDE), match(geneSet,colnames(geneUniverse))])
      
      set.seed(123)
      pval <- chisq.test(mtrix, simulate.p.value = T, B=100000)$p.value
      #pval <- prop.test(mtrix)$p.value
      
      if (geneSet=="ddd"){
        converter <- c(0.75,1.75,2.75)
        names(converter) <- c("D11","D30","D52")
      } else if (geneSet=="cosmic"){
        converter <- c(1.25,2.25,3.25)
        names(converter) <- c("D11","D30","D52")
      } else if (geneSet=="ddd_dominantMOI"){
        converter <- c(1,2,3)
        names(converter) <- c("D11","D30","D52")
      }
      
      dfSumm <- data.frame(annot=y,
                           timepoint=x,
                           pval=pval,
                           geneSet=geneSet,
                           xpos=unname(converter[x]))
      
    }, simplify=F)
    
    annot <- do.call("rbind", annot)
    rownames(annot) <- NULL
    annot
    
  }, simplify=F)
  
  annotTp <- do.call("rbind", annotTp)
  rownames(annotTp) <- NULL
  
  return(annotTp)
  
}

fillNA <- function(neuroTab, list_clustr_ddd, highlight="neuroRelated"){
  colMatch <- match(highlight, colnames(neuroTab))
  alltabs <- sapply(unique(neuroTab$Term), function(x){
    tabs <- sapply(c("allDE","allsignif",names(list_clustr_ddd)), function(y){
      dimensions <- dim(subset(neuroTab, Term==x & label==y))[1]
      if (dimensions==0){
        df <- data.frame(GOBPID=NA,
                         Pvalue="ns",
                         OddsRatio=NA,
                         ExpCount=NA,
                         Count=NA,
                         Size=NA,
                         Term=x,
                         label=y,
                         position=NA,
                         neuroRelated=FALSE,
                         top1=F,
                         top2=F,
                         mostShared=F)
        df[,colMatch] <- TRUE
      } else {
        df <- NA
      }
      return(df)
    }, simplify=F)
    
    tabs <- tabs[!sapply(tabs, function(x) all(is.na(x)))]
    tabs <- do.call("rbind", tabs)
    rownames(tabs) <- NULL
    tabs
  }, simplify=F)
  
  alltabs <- do.call("rbind", alltabs)
  rownames(alltabs) <- NULL
  return(alltabs)
} 




addMissingNA <- function(enrichRes){
  
  days <- unique(enrichRes$tp)
  ctypes <- unique(enrichRes$annot)
  allvec <-  paste0(as.character(sapply(days, function(x) rep(x,length(ctypes)))),"-", rep(ctypes,3))
  
  enrichRes <- sapply(unique(enrichRes$geneSetTested), function(x){
    tmp <- subset(enrichRes, geneSetTested==x)
    missing <- allvec[!allvec %in% tmp$combination2]
    newAdd <- data.frame(matrix(ncol=length(colnames(tmp)), nrow=length(missing)))
    colnames(newAdd) <- colnames(tmp)
    newAdd$combination2 <- missing
    newAdd$tp <- sapply(strsplit(missing, "-"), function(x) x[1])
    newAdd$annot <- sapply(strsplit(missing, "-"), function(x) paste(x[2:length(x)], collapse="-"))
    newAdd$combination <- paste0(sapply(strsplit(missing, "-"), function(x) x[1]), "-", x)
    newAdd$geneSetTested <- x
    rbind(tmp,newAdd)
  }, simplify=F)
  
  enrichRes <- do.call("rbind", enrichRes)
  rownames(enrichRes) <- NULL
  enrichRes
  
}




fillNotEnoughCells <- function(koTab, ctypes){
  
  #ctypes <- unique(koTab$annot)
  
  koTab <- sapply(unique(koTab$gene), function(x){
    #print(x)
    tmp <- subset(koTab, gene==x)
    tpoints <- unique(tmp$tp)
    combinations <- paste0(tpoints,"-",as.character(sapply(ctypes, function(x) rep(x,3))))
    missing <- combinations[!combinations %in% tmp$combinations]
    newAdd <- data.frame(matrix(ncol=length(colnames(tmp)), nrow=length(missing)))
    colnames(newAdd) <- colnames(tmp)
    newAdd$combinations <- missing
    newAdd$gene <- x
    newAdd$donor_extended <- unique(tmp$donor_extended)
    newAdd$annot <- unname(sapply(missing, function(x) sapply(strsplit(x,"-"), function(y) paste(y[2:length(y)],collapse="-"))))
    newAdd$tp <- unname(sapply(missing, function(x) sapply(strsplit(x,"-"), function(y) y[1])))
    newAdd$zscoValue2 <- "nCells<10"
    newAdd2 <- rbind(newAdd, newAdd)
    newAdd2$zscoType <- as.character(sapply(c("Expression","Cell-type fraction"), function(x) rep(x, dim(newAdd)[1])))
    tmp <- rbind(tmp, newAdd2)
    
  }, simplify=F)
  
  koTab <- do.call("rbind", koTab)
  rownames(koTab) <- NULL
  return(koTab)
  
  
}


enrichmentCalc_KO <- function(summResults, geneUniverse, geneSet="ddd"){
  
  annotTp <- sapply(unique(summResults$timepoint), function(x){
    print(x)
    tmp1 <- subset(summResults, timepoint==x)
    
    kodf <- sapply(unique(tmp1$KO), function(z){
      print(z)
      tmp2 <- subset(tmp1, KO==z)
      
      annot <- sapply(unique(tmp2$annot), function(y){
        print(y)
        tmp <- subset(tmp2, annot==y)
        
        if (all(is.na(tmp$numDE))){
          
          dfSumm <- data.frame(annot=y,
                               timepoint=x,
                               pval=NA,
                               geneSet=geneSet,
                               xpos=NA,
                               KO=z,
                               sense=NA)
          
        } else if (unique(tmp$numDE)==0){
          
          dfSumm <- data.frame(annot=y,
                               timepoint=x,
                               pval=NA,
                               geneSet=geneSet,
                               xpos=NA,
                               KO=z,
                               sense=NA)
          
          
        } else {
          
          stopifnot(all(tmp$geneDE %in% geneUniverse$symbol))
          mtrix <- matrix(NA, 2, 2)
          colnames(mtrix) <- c("DE","non-DE")
          rownames(mtrix) <- c("In-GeneSet","Out-GeneSet")
          mtrix[1,1] <- sum(geneUniverse[match(tmp$geneDE, geneUniverse$symbol),match(geneSet,colnames(geneUniverse))])
          mtrix[1,2] <- sum(geneUniverse[!(geneUniverse$symbol %in% tmp$geneDE), match(geneSet,colnames(geneUniverse))])
          mtrix[2,1] <- sum(!geneUniverse[match(tmp$geneDE, geneUniverse$symbol),match(geneSet,colnames(geneUniverse))])
          mtrix[2,2] <- sum(!geneUniverse[!(geneUniverse$symbol %in% tmp$geneDE), match(geneSet,colnames(geneUniverse))])
          set.seed(123)
          pval <- chisq.test(mtrix, simulate.p.value = TRUE, B=50000)$p.value
          
          #print(mtrix)
          #print(pval)
          
          if (pval<0.05){
            sense <- as.character(which.max(sapply(1:ncol(mtrix), function(x) {mtrix[1,x]/mtrix[2,x]})))
          } else {
            sense <- NA
          }
          
          vecSense <- c("more","less")
          names(vecSense) <- c("1","2")
          
          if (geneSet=="ddd"){
            converter <- (1:length(sort(unique(tmp1$KO))))-0.25
            names(converter) <- sort(unique(tmp1$KO))
          } else if (geneSet=="cosmic"){
            converter <- (1:length(sort(unique(tmp1$KO))))+0.25
            names(converter) <- sort(unique(tmp1$KO))
          } else if (geneSet=="ddd_dominantMOI"){
            converter <- (1:length(sort(unique(tmp1$KO))))
            names(converter) <- sort(unique(tmp1$KO))
          }
          
          dfSumm <- data.frame(annot=y,
                               timepoint=x,
                               pval=pval,
                               geneSet=geneSet,
                               xpos=unname(converter[z]),
                               KO=z,
                               sense=unname(vecSense[sense]))
          
        }
        
        
        
      }, simplify=F)
      
      annot <- do.call("rbind", annot)
      rownames(annot) <- NULL
      annot
      
    }, simplify=F)
    
    kodf <- do.call("rbind", kodf)
    rownames(kodf) <- NULL
    kodf
    
  }, simplify=F)
  
  annotTp <- do.call("rbind", annotTp)
  rownames(annotTp) <- NULL
  annotTp
  return(annotTp)
  
}

addPosition <- function(report_allclusters){
  report_allclusters$position <- 1:dim(report_allclusters)[1]
  return(report_allclusters)
}


fillNA_ko <- function(neuroTab, allTestsKO, highlight="neuroRelated"){
  colMatch <- match(highlight, colnames(neuroTab))
  alltabs <- sapply(unique(neuroTab$Term), function(x){
    tabs <- sapply(allTestsKO, function(y){
      dimensions <- dim(subset(neuroTab, Term==x & label==y))[1]
      if (dimensions==0){
        df <- data.frame(GOBPID=NA,
                         Pvalue="ns",
                         OddsRatio=NA,
                         ExpCount=NA,
                         Count=NA,
                         Size=NA,
                         Term=x,
                         label=y,
                         position=NA,
                         neuroRelated=FALSE)
        df[,colMatch] <- TRUE
      } else {
        df <- NA
      }
      return(df)
    }, simplify=F)
    
    tabs <- tabs[!sapply(tabs, function(x) all(is.na(x)))]
    tabs <- do.call("rbind", tabs)
    rownames(tabs) <- NULL
    tabs
  }, simplify=F)
  
  alltabs <- do.call("rbind", alltabs)
  rownames(alltabs) <- NULL
  return(alltabs)
} 



functionReproTab <- function(summTab, repro="bioRep"){
  
  colMatch <- colnames(summTab)[grep(repro, colnames(summTab))]
  tmp <- subset(summTab, get(colMatch)!=0)
  
  tmpLine <- sapply(unique(tmp$cell_lines), function(x) {
    
    tmp2 <- subset(tmp, cell_lines==x)
    
    if (repro=="bioRep"){
      
      tmp2[,c(colMatch)] <- gsub(".+-","", tmp2[,c(colMatch)])
      stopifnot(length(unique(tmp2[,c(colMatch)]))==2)
      data.frame(cell_lines=x,
                 tp=unique(tmp2$tp),
                 rep.x=sapply(split(tmp2$norm_numCells, tmp2[,c(colMatch)]), mean)[["1"]],
                 rep.y=sapply(split(tmp2$norm_numCells, tmp2[,c(colMatch)]), mean)[["2"]],
                 rep.x.num=table(tmp2[,c(colMatch)])[["1"]],
                 rep.y.num=table(tmp2[,c(colMatch)])[["2"]],
                 type=repro)
      
    } else if (repro=="techRep" | repro=="tenXRep") {
      
      tmp2$newCol <- gsub("-.+","",tmp2[,c(colMatch)])
      tmpRep <- sapply(unique(tmp2$newCol), function(y){
        tmp3 <- subset(tmp2, newCol==y)
        tmp3[,c(colMatch)] <- gsub(".+-","", tmp3[,c(colMatch)])
        stopifnot(length(unique(tmp3[,c(colMatch)]))==2)
        data.frame(cell_lines=x,
                   tp=unique(tmp3$tp),
                   rep.x=sapply(split(tmp3$norm_numCells, tmp3[,c(colMatch)]), mean)[["1"]],
                   rep.y=sapply(split(tmp3$norm_numCells, tmp3[,c(colMatch)]), mean)[["2"]],
                   rep.x.num=table(tmp3[,c(colMatch)])[["1"]],
                   rep.y.num=table(tmp3[,c(colMatch)])[["2"]],
                   type=repro)
        
      }, simplify=F)
      
      tmpRep <- do.call("rbind", tmpRep)
      rownames(tmpRep) <- NULL
      tmpRep
      
      
    } else if (repro=="donorRep"){
      
      vec <- c("1","2","3")
      names(vec) <- c("D11","D30","D52")
      tmp2[,c(colMatch)] <- paste0(unname(vec[tmp2$tp]),"-",as.character(tmp2[,c(colMatch)]))
      
      tmp2$newCol <- gsub("-.+","",tmp2[,c(colMatch)])
      realdonRep <- tmp2[,c(colMatch)][grepl("-2", tmp2[,c(colMatch)])]
      realdonRep <- c(unique(gsub("-2","-1", realdonRep))[unique(gsub("-2","-1", realdonRep)) %in% tmp2[,c(colMatch)]],
                      gsub("-1","-2", unique(gsub("-2","-1", realdonRep))[unique(gsub("-2","-1", realdonRep)) %in% tmp2[,c(colMatch)]]))
      tmp2 <- tmp2[tmp2[,c(colMatch)] %in% realdonRep,]
      
      tmpRep <- sapply(unique(tmp2$newCol), function(y){
        tmp3 <- subset(tmp2, newCol==y)
        tmp3[,c(colMatch)] <- gsub(".+-","", tmp3[,c(colMatch)])
        
        stopifnot(length(unique(tmp3[,c(colMatch)]))==2)
        data.frame(cell_lines=x,
                   tp=unique(tmp3$tp),
                   rep.x=sapply(split(tmp3$norm_numCells, tmp3[,c(colMatch)]), mean)[["1"]],
                   rep.y=sapply(split(tmp3$norm_numCells, tmp3[,c(colMatch)]), mean)[["2"]],
                   rep.x.num=table(tmp3[,c(colMatch)])[["1"]],
                   rep.y.num=table(tmp3[,c(colMatch)])[["2"]],
                   type=repro)
        
      }, simplify=F)
      
      tmpRep <- do.call("rbind", tmpRep)
      rownames(tmpRep) <- NULL
      tmpRep
      
    }
    
  }, simplify=F)
  
  tmpLine <- do.call("rbind", tmpLine)
  rownames(tmpLine) <- NULL
  return(tmpLine)
  
  
}


computeBurdenGeneSet <- function(allvariants, varCategory, geneset){
  
  if (varCategory=="all"){
    selectVariants <- allvariants
  } else if (varCategory=="deleterious"){
    selectVariants <- subset(allvariants, varCategoryFiltered=="missPatho" | varCategoryFiltered=="ptv")
  } else {
    selectVariants <- subset(allvariants, varCategoryFiltered==varCategory)
  }
  
  samples <- unique(allvariants$ips)
  genes <- unique(allvariants$gene_name)
  mtrixBurden <- matrix(0, nrow = length(genes), ncol = length(samples))
  rownames(mtrixBurden) <- genes
  colnames(mtrixBurden) <- unique(allvariants$ips)
  
  if (geneset!="none"){
    selectVariants <- subset(selectVariants, get(geneset)==T)
  }
  
  tmp <- sapply(samples, function(x) table(subset(selectVariants, selectVariants$ips==x)$gene_name))
  for (i in 1:length(tmp)){
    if (length(tmp[[i]])>0){
      mtrixBurden[match(names(tmp[[i]]), rownames(mtrixBurden)),
                  match(names(tmp)[i], colnames(mtrixBurden))] <- tmp[[i]]
    }
  }
  burden_obj <- as.data.frame(colSums(mtrixBurden))
  burden_obj$cell_lines <- rownames(burden_obj)
  rownames(burden_obj) <- NULL
  burden_obj <- burden_obj[,rev(colnames(burden_obj))]
  colnames(burden_obj)[2] <- paste0(varCategory)
  
  return(burden_obj)
  
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
