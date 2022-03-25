library(harmony)
library(dplyr)

fileFolder="output/"

## 50 first PCA per sample (before Harmony)
pca_file=paste0(fileFolder,"allpools.scanpy.dimreduction.PCA.scaled.tsv")

## 50 first PCA (after Harmony)
out_file=paste0(fileFolder,"allpools.scanpy.dimreduction.harmonyPCA.scaled.tsv")

## Metadata table (columns to take into account with batch correction)
metadata_file=paste0(fileFolder,"allpools.scanpy.dimreduction.scaled.obs_df.tsv")

## sample-id indicates the different 10x runs
## we pretend to correct the batch effect among different 10xruns
merge_column="sample_id"
theta=2

pca_df = read.table(pca_file, sep='\t', header=TRUE, row.names=1)
metadata_df = read.table(metadata_file, sep='\t', header=TRUE)
rownames(metadata_df) = metadata_df$index

#max.iter.harmony = 10, max.iter.cluster = 200
harmony_embeddings <- HarmonyMatrix(pca_df, metadata_df, merge_column, theta=theta, do_pca=FALSE, max.iter.harmony = 25, max.iter.cluster = 500, plot_convergence=TRUE)
write.table(harmony_embeddings, out_file, sep='\t', quote=FALSE)

sessionInfo()