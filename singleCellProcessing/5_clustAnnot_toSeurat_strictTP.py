
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
sc.settings.verbosity = 2
sc.logging.print_header()
sc.settings.figdir = "figuresQC/"

pathToDir="output/"
filename=str(pathToDir)+"allpools.scanpy.harmonyPCA.clust.wmeta.notNorm.h5ad"
adata=sc.read_h5ad(filename)
print('{} cells, {} genes'.format(adata.n_obs, adata.n_vars))
adata.obs['batch'] = adata.obs['batch'].astype('str').astype('category')
adata.obs["sample_id"]=adata.obs["sample_id"].astype("str").astype('category')
adata.obs["donor_id"]=adata.obs["donor_id"].astype("str").astype('category')
adata.obs["pool_id"]=adata.obs["pool_id"].astype("str").astype('category')
adata.obs["time_point"]=adata.obs["time_point"].astype("str").astype('category')
adata.obs["treatment"]=adata.obs["treatment"].astype("str").astype('category')
adata.obs["tp_treat"]=adata.obs["tp_treat"].astype("str").astype('category')
adata.var["gene_ids-0"]=adata.var["gene_ids-0"].astype("str")
adata.obs["leiden_id"]=adata.obs["leiden_id"].astype("str").astype('category')
adata.obs.index=adata.obs.index.astype("str")
adata.var.index=adata.var.index.astype("str")
adata.obs["tp_treat"]=adata.obs["tp_treat"].replace({"D52-no":"D52","D52-yes":"D52R", "D11-no":"D11","D30-no":"D30","D119-no":"D119"})
adata.obs["clust_TP"]=adata.obs['leiden_id'].astype(str)+"-"+adata.obs['tp_treat'].astype(str)


sc.pp.filter_cells(adata, min_genes=200)
print('{} observations, {} genes'.format(adata.n_obs, adata.n_vars))
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

results_file=str(pathToDir)+"allpools.scanpy.wMetaClustUmapGraph.exprLogNormNotScaled.h5ad"
adata.write_h5ad(results_file)

## split dataset per timepoint (logNorm counts)
tpoints=adata.obs["tp_treat"].unique()
for tpoint in tpoints:
    print(tpoint)
    tmp=adata[adata.obs["tp_treat"] == tpoint, :]
    print('{} observations, {} genes'.format(tmp.n_obs, tmp.n_vars))
    new_h5file=str(pathToDir)+"allpools.scanpy."+str(tpoint)+".wMetaClustUmapGraph.exprLogNormNotScaled.h5ad"
    tmp.write_h5ad(new_h5file)
    # sc.pp.scale(tmp, max_value=1)
    # newnew_h5file=str(pathToDir)+"allpools.scanpy."+str(tpoint)+".wMetaClustUmapGraph.exprLogNormScaled.h5ad"
    # tmp.write_h5ad(newnew_h5file)

sc.pp.scale(adata, max_value=1)

############################
### Cell type annotation ###
############################

astrocyte=list(pd.read_csv("annotgenes/astrocyte.txt", header=None).loc[:,0])
cortical=list(pd.read_csv("annotgenes/cortical.txt", header=None).loc[:,0])
dopaminergic=list(pd.read_csv("annotgenes/dopaminergic.txt", header=None).loc[:,0])
ependymal=list(pd.read_csv("annotgenes/ependymal.txt", header=None).loc[:,0])
floorPlateProg=list(pd.read_csv("annotgenes/floorPlateProg.txt", header=None).loc[:,0])
glia=list(pd.read_csv("annotgenes/glia.txt", header=None).loc[:,0])
neuroblasts=list(pd.read_csv("annotgenes/neuroblasts.txt", header=None).loc[:,0])
neuronalProg=list(pd.read_csv("annotgenes/neuronalProg.txt", header=None).loc[:,0])
neuronal=list(pd.read_csv("annotgenes/neuronal.txt", header=None).loc[:,0])
proliferative_floorPlateProg=list(pd.read_csv("annotgenes/proliferative_floorPlateProg.txt", header=None).loc[:,0])
serotonergic=list(pd.read_csv("annotgenes/serotonergic.txt", header=None).loc[:,0])
targeted_floorPlate=list(pd.read_csv("annotgenes/targeted_floorPlate.txt", header=None).loc[:,0])

#average expression profile of canonical marker genes across cell types (expression scaled between 0 and 1)
#heatmap1

dict_heatmap1 = {   'Prolif': proliferative_floorPlateProg,
                    'Neur_Prog': neuronalProg,
                    'FP_prog': floorPlateProg,
                    'Glia': glia,
                    'FP': targeted_floorPlate,
                    'Neuroblasts':neuroblasts,
                    'Astrocyte':astrocyte,
                    'Ependyma':ependymal,
                    'Serotonergic':serotonergic,
                    'Neur':neuronal }

# sc.pl.matrixplot(adata, dict_heatmap1, 'leiden_id', dendrogram=False, cmap='Blues', standard_scale='var',
#                  save="celltypes_neuroseq_manual_cellAnnotation.png")

sc.pl.matrixplot(adata, dict_heatmap1, 'clust_TP', dendrogram=False, cmap='Blues', standard_scale='var',
                 save="suppfig2E.png")

#average expression profile of canonical marker for dopaminergic neurons (literature-curated)
#heatmap2, Dopaminergic

# sc.pl.matrixplot(adata, dopaminergic, 'leiden_id', dendrogram=False, cmap='Blues', standard_scale='var', title="Dopaminergic markers",
#                  save="dopaminergic_neuroseq_manual_cellAnnotation.png")

sc.pl.matrixplot(adata, dopaminergic, 'clust_TP', dendrogram=False, cmap='Blues', standard_scale='var', title="Dopaminergic markers",
                 save="suppfig2F.png")

#average expression profile of canonical marker for cortical hem/Cajal retzius
#heatmap3, Cortical

# sc.pl.matrixplot(adata, cortical, 'leiden_id', dendrogram=False, cmap='Blues', standard_scale='var', title="Cortical hem/Cajal retzius",
#                  save="cortical_neuroseq_manual_cellAnnotation.png")

sc.pl.matrixplot(adata, cortical, 'clust_TP', dendrogram=False, cmap='Blues', standard_scale='var', title="Cortical hem/Cajal retzius",
                 save="suppfig2G.png")

