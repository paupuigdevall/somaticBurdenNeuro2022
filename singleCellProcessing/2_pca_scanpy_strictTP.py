import os
import numpy as np
import pandas as pd
import scanpy as sc
import h5py
import matplotlib.pyplot as plt
import seaborn as sb

from collections import Counter
from matplotlib import rcParams

sc.set_figure_params(dpi=200, color_map='viridis')
sc.logging.print_header()

sc.settings.figdir = 'figuresQC/'
pathDir="output/"
filename=str(pathDir)+"allpools.scanpy.init.h5ad"
adata=sc.read(filename)
adata.var_names_make_unique()
print('{} cells, {} genes'.format(adata.n_obs, adata.n_vars))

######## Move from bytes to str

adata.obs['batch'] = adata.obs['batch'].astype('str').astype('category')
adata.obs["sample_id"]=adata.obs["sample_id"].astype("str").astype('category')
adata.obs["donor_id"]=adata.obs["donor_id"].astype("str").astype('category')
adata.obs["pool_id"]=adata.obs["pool_id"].astype("str").astype('category')
adata.obs["time_point"]=adata.obs["time_point"].astype("str").astype('category')
adata.obs["treatment"]=adata.obs["treatment"].astype("str").astype('category')
adata.obs["tp_treat"]=adata.obs["tp_treat"].astype("str")
adata.var["gene_ids-0"]=adata.var["gene_ids-0"].astype("str")
adata.obs.index=adata.obs.index.astype("str")
adata.var.index=adata.var.index.astype("str")
adata.obs["tp_treat"]=adata.obs["tp_treat"].replace({"D52-no":"D52","D52-yes":"D52+R", "D11-no":"D11","D30-no":"D30","D119-no":"D119"})

################
### QC plots ###
################

### ribo counts
ribo_genes = adata.var_names.str.startswith(("RPS","RPL"))
print("Number of ribosomal genes")
print(sum(ribo_genes))
adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mito','percent_ribo'],
             jitter=0.4, groupby = 'tp_treat', size=0.1, save="suppfig2D.png")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mito', color="tp_treat", size=1, alpha=0.5, save="scatter_totCounts_mito.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="tp_treat", size=1, alpha=0.5, save="scatter_totCounts_nGenes.png")
sc.pl.scatter(adata, x='total_counts', y='percent_ribo', color="tp_treat", size=1, alpha=0.5, save="scatter_totCounts_ribo.png")
sc.pl.scatter(adata, x='pct_counts_mito', y='percent_ribo', color="tp_treat", size=1, alpha=0.5, save="scatter_mito_ribo.png")

## plot highly variable by timepoint
tpoints=adata.obs["tp_treat"].unique()
for tpoint in tpoints:
    print(tpoint)
    tmp_adata = adata[adata.obs["tp_treat"]==tpoint]
    figname=str(tpoint)+"_top50expGenes.png"
    a=sc.pl.highest_expr_genes(tmp_adata, n_top=50, show=False)
    a.set_xlim(0,10)
    tt=a.get_figure()
    plotname=str(sc.settings.figdir)+str(figname)
    tt.savefig(plotname)


pathToSave="figuresQC/"

## Thresholding decision: counts
plt.figure()
filename=pathToSave+"distr_tcounts.png"
p3 = sb.distplot(adata.obs['total_counts'], kde=False)
p3.set_yticklabels([str(i) for i in p3.get_yticks()], fontsize = 8)
p3.set(title="Total counts")
p3.figure.savefig(filename)

plt.figure()
filename=pathToSave+"distr_tcounts_st8000.png"
p4 = sb.distplot(adata.obs['total_counts'][adata.obs['total_counts']<8000], kde=False, bins=60)
p4.set_yticklabels([str(i) for i in p4.get_yticks()], fontsize = 8)
p4.set(title="Total counts (<8000)")
p4.figure.savefig(filename)

plt.figure()
filename=pathToSave+"distr_tcounts_gt10000.png"
p5 = sb.distplot(adata.obs['total_counts'][adata.obs['total_counts']>10000], kde=False, bins=60)
p5.set_yticklabels([str(i) for i in p5.get_yticks()], fontsize = 8)
p5.set(title="Total counts (>10000)")
p5.figure.savefig(filename)

## Thresholding decision: genes
plt.figure()
filename=pathToSave+"distr_ngenes.png"
p6 = sb.distplot(adata.obs['n_genes_by_counts'], kde=False, bins=60)
p6.figure.savefig(filename)

plt.figure()
filename=pathToSave+"distr_ngenes_st1000.png"
p7 = sb.distplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts']<1000], kde=False, bins=60)
p7.figure.savefig(filename)

### Filter genes by counts
# only consider genes expressed in more than 0.1% of cells
min_cells=0.001*len(adata.obs.index)
sc.pp.filter_genes(adata, min_cells=min_cells)
print('{} observations, {} genes'.format(adata.n_obs, adata.n_vars))

## Calculate cell-cycle scores
filename="regev_lab_cell_cycle_genes.txt"
cell_cycle_genes = [x.strip() for x in open(filename)]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))

sc.pp.normalize_per_cell(adata, key_n_counts="total_counts")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)

print(Counter(adata.var['highly_variable']))

## plot highly variable
highVarPlot="highlyVarGenes.png"
sc.pl.highly_variable_genes(adata, save=highVarPlot)
adata_highVar = adata[:, adata.var.highly_variable]
print('{} observations, {} genes'.format(adata_highVar.n_obs, adata_highVar.n_vars))

sc.pp.scale(adata_highVar)

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pl.violin(adata, ['S_score', 'G2M_score'],
             jitter=0.4, groupby = 'tp_treat', size=0.1, save="cell_cycle_scores.png")

c = adata.obs.groupby(['tp_treat', 'phase']).size().rename("count")
c = (100*c) / c.groupby(level=0).sum()

print(c)
output_file=str(pathDir)+"cellCycle_perTimePoint.csv"
c.to_csv(output_file, sep="\t")
adata.obs['cycleTP'] = adata.obs[['tp_treat', 'phase']].apply(lambda x: ''.join(x), axis=1)
adata.obs['cycleTP'] = adata.obs['cycleTP'].astype("str").astype("category")

to_extract = ["batch","donor_id", "sample_id","pool_id","time_point","treatment", "total_counts","tp_treat","percent_ribo","phase"]
metadata_opt = pd.DataFrame(adata.obs[to_extract],
                           index=adata.obs_names,
                           columns=to_extract)
metadata_optfile=str(pathDir)+"metadata_withPercentRibo_CellCyclePhase.tsv"
metadata_opt.to_csv(metadata_optfile, sep="\t")


## produce plots
adata.obs['cycleTP'] = adata.obs[['time_point', 'phase']].apply(lambda x: ''.join(x), axis=1)
sc.set_figure_params(dpi=100, color_map='viridis', figsize=(10,4))
sc.pl.violin(adata, ['percent_ribo'],
             jitter=0.4, groupby = 'cycleTP', size=0.1, save="perc_ribo_time_point_cell_cycle.png")


#########
## PCA ##
#########

print('PCA')
n_pcs=50
sc.tl.pca(adata_highVar, n_comps=n_pcs, svd_solver='arpack')
pca_df = pd.DataFrame(adata_highVar.obsm['X_pca'], index=adata_highVar.obs_names, columns = ['PC{}'.format(x) for x in range(1,n_pcs+1)])

pca_outfile=str(pathDir)+"allpools.scanpy.dimreduction.PCA.scaled.tsv"
dataset_outfile=str(pathDir)+"allpools.scanpy.dimreduction.PCA.scaled.h5"

pca_df.to_csv(pca_outfile, sep='\t')
adata_highVar.write(dataset_outfile)

## Metadata table for Harmony

to_extract = ["batch","donor_id", "sample_id","pool_id","time_point","treatment", "total_counts","tp_treat"]
metadata_df = pd.DataFrame(adata_highVar.obs[to_extract],
                           index=adata_highVar.obs_names,
                           columns=to_extract)

metadata_outfile=str(pathDir)+"allpools.scanpy.dimreduction.scaled.obs_df.tsv"
metadata_df.to_csv(metadata_outfile, sep="\t")

