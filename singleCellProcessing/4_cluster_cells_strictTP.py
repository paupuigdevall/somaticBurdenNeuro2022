
import os
import numpy as np
import pandas as pd
import scanpy as sc
import h5py
import matplotlib.pyplot as plt
import yaml
import fa2

from collections import Counter
from matplotlib import rcParams

sc.set_figure_params(dpi=200, color_map='viridis')
sc.settings.verbosity = 2
sc.logging.print_header()
sc.settings.figdir = "figuresQC/"

yamlfile ="allpools_mito5.yaml"
with open(yamlfile, "r") as f:
    config_dict = yaml.safe_load(f)


## variables to be used
pathToDir="output/"
in_file=str(pathToDir)+"allpools.scanpy.dimreduction.PCA.scaled.h5"
pca_file=str(pathToDir)+"allpools.scanpy.dimreduction.harmonyPCA.scaled.tsv"
print(pca_file)
pca_before_harmony=str(pathToDir)+"allpools.scanpy.dimreduction.PCA.scaled.tsv"
print(pca_before_harmony)
resolution = config_dict["resolution"]

pca_df = pd.read_csv(pca_file, sep='\t', index_col=0)

##h5 file with the 3,000 highly variable genes
adata = sc.read(in_file)

adata.obs['batch'] = adata.obs['batch'].astype('str').astype('category')
adata.obs["sample_id"]=adata.obs["sample_id"].astype("str").astype('category')
adata.obs["donor_id"]=adata.obs["donor_id"].astype("str").astype('category')
adata.obs["pool_id"]=adata.obs["pool_id"].astype("str").astype('category')
adata.obs["time_point"]=adata.obs["time_point"].astype("str").astype('category')
adata.obs["treatment"]=adata.obs["treatment"].astype("str").astype('category')
adata.obs["tp_treat"]=adata.obs["tp_treat"].astype("str").astype('category')
adata.var["gene_ids-0"]=adata.var["gene_ids-0"].astype("str")
adata.obs.index=adata.obs.index.astype("str")
adata.var.index=adata.var.index.astype("str")


sc.pl.pca_variance_ratio(adata, log=True, save="varianceCorrected.png")
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.3)

## Comparative batch plot
plotname=str(sc.settings.figdir)+"/"+str("evaluateBatch.png")
fig, axs = plt.subplots(2, 1, figsize=(8,8),constrained_layout=True)
sc.pl.umap(adata, color="leiden", title="UMAP Corr w/ Harmony", ax=axs[0,0], show=False)
sc.pl.umap(adata, color="pool_id", title="Batch corr", ax=axs[1,0], legend_loc="none", show=False)
fig.savefig(plotname)


## Other QC plots for corrected batch
sc.pl.umap(adata, color="tp_treat", save="umap_tptreat.png")
sc.pl.umap(adata, color=['leiden','time_point','pool_id'], save="umap_clust_tp_poolid.png", legend_loc="none")
sc.pl.umap(adata, color=['leiden','tp_treat','pool_id'], save="umap_clust_tptreat_poolid.png", legend_loc="none")

## PAGA trajectories and UMAP with PAGA positions
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, color="leiden", save="paga_clustDistribution.png", title="PAGA connections")
sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color="leiden", save="suppfig2H.png")
sc.pl.umap(adata, color="time_point", save="suppfig2I.png")
sc.pl.umap(adata, color=["leiden", "time_point"], save="umap_posPaga_clust_tp_tptreat.png", legend_loc="none" )

## Coarse-grained graph
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden','time_point'], legend_loc='on data', save="graph_posPaga_clust_tp.png")
# sc.pl.draw_graph(adata, color=['leiden','tp_treat'], legend_loc='on data', save="graph_posPaga_clust_tptreat.png")


#table with leiden clusterId
cluster_outfile=str(pathToDir)+"allpools.scanpy.dimreduction.harmonyPCA.cell_clust_df.tsv"
cell_clustering_df = adata.obs[['leiden']]
cell_clustering_df.to_csv(cluster_outfile, sep='\t')

## Add metadata with initial scanpy_h5
h5_input=str(pathToDir)+"allpools.scanpy.init.h5ad"
h5_output=str(pathToDir)+"allpools.scanpy.harmonyPCA.clust.wmeta.notNorm.h5ad"
umap_outfile=str(pathToDir)+"umap_outfile.tsv"
graph_fa_df_outfile=str(pathToDir)+"graph_fa.tsv"

umap_mat = adata.obsm["X_umap"]
umap_cell_index = adata.obs.index
umap_df = pd.DataFrame(data=umap_mat, index=umap_cell_index)
umap_df.to_csv(umap_outfile, sep="\t")

graph_fa_mat = adata.obsm["X_draw_graph_fa"]
graph_fa_cell_index = adata.obs.index
graph_fa_df = pd.DataFrame(data=graph_fa_mat, index=graph_fa_cell_index)
graph_fa_df.to_csv(graph_fa_df_outfile, sep="\t")

init_h5 = sc.read(h5_input)
init_h5.obs.index=init_h5.obs.index.astype("str")
init_h5.var.index=init_h5.var.index.astype("str")

init_h5.obs['leiden_id'] = adata.obs[["leiden"]]
init_h5.obsm['X_umap'] = umap_df.loc[init_h5.obs.index,:].values
init_h5.obsm['X_draw_graph_fa'] = graph_fa_df.loc[init_h5.obs.index,:].values

umap_df=None
graph_fa_df=None

init_h5.write_h5ad(h5_output)



