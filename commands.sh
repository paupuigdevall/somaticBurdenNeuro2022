
githubrepo="Path_to_Github_cloned_repository"

## Import libraries / SessionInfo()
cd ${githubrepo}
Rscript libraryInstall.R

##################
### CellRanger ###
##################

## Launch cellranger for each 10x sample
cd ${githubrepo}/cellranger
bsub -n32 -R'select[mem>200000] rusage[mem=200000] span[ptile=8]' -M200000 -J "cellrangerhg19[1-119]" -q long -o log/cellranger.%J.%I.hg19.out -e log/cellranger.%J.%I.hg19.err ./cellranger_launcher.sh pathTo10xFastQSamples.txt

################
### Demuxlet ###
################

## Generate VCF input for demuxlet
cd ${githubrepo}/demuxlet/vcfGenerator
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -o log/gnomadDownload.out -e log/gnomadDownload.err ./0.gnomAdRegionsGenerator.sh
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -J "filterVCF[1-22]" -o log/output.%J.%I -e log/error.%J.%I ./1.filterVCF.sh demuxlet/vcfGenerator/chrlist.txt 
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -o log/concatenate.out -e log/concatenate.err ./2.concatenateVCF.sh
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -o log/finalVCF.out -e log/finalVCF.err ./3.finalVCF.sh

## Run demuxlet (cell-donor deconvolution)
cd ${githubrepo}/demuxlet/deconvolution
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -J "demuxlethg19[1-119]" -q long -o log/demuxlet.%J.%I.hg19.out -e log/demuxlet.%J.%I.hg19.err ./demuxlet_launcher_ALL_hg19.sh pathToPossortedBam_And_PoolID_per10xSample.txt
cut -f1 pathToPossortedBam_And_PoolID_per10xSample.txt > inputQC.txt
bsub -R'select[mem>2000] rusage[mem=2000]' -M2000 -J "flagQC[1-119]" -o log/flagstatQC.%J.%I.out -e log/flagstatQC.%J.%I.err ./flagstat.sh inputQC.txt 


##############################
### Single-cell processing ###
##############################

cd ${githubrepo}/singleCellProcessing
### Merge and QC individual 10x samples
bsub -R'select[mem>80000] rusage[mem=80000]' -M80000 -o log/mergeScanpy.out -e log/mergeScanpy.err ./1_merge_scanpy_strictTP.sh   
### Run PCA with merged object
bsub -R'select[mem>80000] rusage[mem=80000]' -M80000 -o log/runPCA.out -e log/runPCA.err ./2_pca_scanpy_strictTP.sh  
### Run Harmony (with R)
bsub -R'select[mem>25000] rusage[mem=25000]' -M25000 -o log/harmony.out -e log/harmony.err ./3_run_harmony_strictTP.sh
### UMAP+Clustering
bsub -R'select[mem>80000] rusage[mem=80000]' -M80000 -q long -o log/clust.out -e log/clust.err ./4_cluster_cells_strictTP.sh
### Saving h5ad + annotation
bsub -R'select[mem>120000] rusage[mem=120000]' -M120000 -q long -o log/annot.out -e log/annot.err ./5_clustAnnot_toSeurat_strictTP.sh

##########################################
### Analysis & Figures reproducibility ###
##########################################

##R version 4.0.4
## install packages and sessionInfo()


cd ${githubrepo}/analysis
Rscript 1-donor_ctype_perPool.R
Rscript 2-curationMutAcqInVivo.R
Rscript 3-diffEfficiency_tab.R
Rscript 4-geneUniverseAnnotated.R
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -J "cLineVarInfo[1-832]" -q long -o log/cLineVarInfo.%J.%I.out -e log/cLineVarInfo.%J.%I.err ./5-cellLine_varInfo.sh listOFfilesVEP.txt
Rscript 6-varCountPerLineAndCategoryBuilder.R
Rscript 7-mutTab_logRegBuilder.R
Rscript 8-enrichmentOntology.R
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -J "pseudoBulkBCOR[1-3]" -q long -o log/pseudoBulkBCOR.%J.%I.out -e log/pseudoBulkBCOR.%J.%I.err ./9a-pseudoBulkTransform.sh timepoint.txt
Rscript 9b-bcorMutationsInFailed.R
Rscript 9c-repliTimeScript.R
bsub -R'select[mem>5000] rusage[mem=5000]' -M5000 -J "pseudoBulkLine[1-3]" -q long -o log/pseudoBulkLine.%J.%I.out -e log/pseudoBulkLine.%J.%I.err ./9d-calculateAvgPseudoExpPerLine.sh timepoint.txt

Rscript 10-reproTab_cfracPerBatch.R
Rscript 11-growthRate_outlierDef.R
Rscript 12-QC_check.R

## Main figures 2,3
cd ${githubrepo}
Rscript figure2.R
Rscript figure3.R

## Supplementary figures 1,2,3
Rscript suppFigure1.R
Rscript suppFigure2.R
Rscript suppFigure3.R


cd ${githubrepo}/analysis
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "DEperTimePoint[1-3]" -o log/DEperTimePoint.%J.%I.out -e log/DEperTimePoint.%J.%I.err ./13-masterOutcomeDE.sh timepoint.txt 
bsub -R'select[mem>80000] rusage[mem=80000]' -M80000 -J "seuratList[1-3]" -o log/seuratList.%J.%I.out -e log/seuratList.%J.%I.err ./14-matrixGeneDonorGeneratorTP.sh timepoint.txt 
Rscript 15-geneUniverseBuilder_scRNAseq.R
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "cosmicEnrich[1-7]" -o log/cosmicEnrich.%J.%I.out -e log/cosmicEnrich.%J.%I.err ./16-cosmicEnrichment.sh cosmicEnrichedTP.txt 


## Main figure 4
cd ${githubrepo}
Rscript figure4.R

## Supplementary figure 4
Rscript suppFigure4.R

cd ${githubrepo}/analysis
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "zscoCorrelationD11[1-26]" -o log/zscoCorrelationD11.%J.%I.out -e log/zscoCorrelationD11.%J.%I.err ./17-corr_expression.sh D11_indices_run.txt 
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "zscoCorrelationD30[1-29]" -o log/zscoCorrelationD30.%J.%I.out -e log/zscoCorrelationD30.%J.%I.err ./17-corr_expression.sh D30_indices_run.txt 
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "zscoCorrelationD52[1-30]" -o log/zscoCorrelationD52.%J.%I.out -e log/zscoCorrelationD52.%J.%I.err ./17-corr_expression.sh D52_indices_run.txt 
bsub -R'select[mem>60000] rusage[mem=60000]' -M60000 -J "exampleCorr[1-4]" -o log/exampleCorr.%J.%I.out -e log/exampleCorr.%J.%I.err ./18-filterGeneTPAnnot.sh exampleCorrRuns.txt 

## Main figure 5
cd ${githubrepo}
Rscript figure5.R

## Supplementary figure 5
Rscript suppFigure5.R






