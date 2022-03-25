import pandas as pd
import scanpy as sc
import numpy as np
import scipy
import glob
import os
import argparse
import yaml
import re
import sys
from functools import reduce

# memory efficiency step - exclude genes on the fly if we already know which genes to include
genelist_file="gene_list.csv"
gene_list = pd.read_csv(genelist_file, sep="\t", header=None)[1]

# loading 10x selected samples
selected_samples_file="pools_to_merge_strictTP.csv"
selected_samples = pd.read_csv(selected_samples_file)

raw_directory = "/projects/NeuroSeq/scRNAseqData/fastq/"
selected_samples['path_to_h5_file'] = selected_samples['path_to_files'].apply(lambda x: os.path.join(raw_directory, '{}filtered_feature_bc_matrix.h5'.format(x)))
selected_samples['path_to_demuxlet_file'] = selected_samples['path_to_files'].apply(lambda x: os.path.join(raw_directory, '{}output.demuxlet.doublet0.5.best'.format(x)))

selected_samples['file_exists'] = selected_samples['path_to_h5_file'].apply(lambda x: os.path.exists(x))
if selected_samples['file_exists'].sum()>0:
    print('h5 files missing for: {}'.format(selected_samples.query('not file_exists')['path_to_files'].tolist()))
# only keep files that exist
selected_samples = selected_samples.query('file_exists')

yamlfile = "allpools_mito5.yaml"

with open(yamlfile, "r") as f:
    config_dict = yaml.safe_load(f)

print(config_dict)
print(type(config_dict["donor_filter"]))

datasets = []
count = 1
total_cells_counter = []

# pool 10-D11 have already been excluded (not even processed by cellranger)
# exclude poor 10x samples (pool12 D52 (x2), pool1 D30 (x1), pool8 D30 (x1), )
excluded_samples = ['run29749/pool12-BIO1-TEC2-TEN1-D52/5245STDY7982705/cellranger-hg19/outs',
                    'run29749/pool12-BIO1-TEC1-TEN1-D52/5245STDY7982704/cellranger-hg19/outs/',
                    'run25528/pool1-BIO1-TEC2-TEN1-D30/5245STDY7387188/cellranger-hg19/outs/',
                    'run27007/pool8-BIO1-TEC1-TEN2-D30/5245STDY7654677/cellranger-hg19/outs/']
                    
selected_samples = selected_samples.query('path_to_files not in @excluded_samples')


# iterate over samples, checking whether it should be included, then adding it to the dataset with donor ID info
for idx,row in selected_samples.iterrows():
    
    h5_file = row['path_to_h5_file']
    demuxlet_file = row['path_to_demuxlet_file']
    print(count)
    
    sample_id = row['sample_id']
    print('sample: ' + sample_id)

    # read the donor mapping information
    donor_df = pd.read_csv(demuxlet_file, sep='\t', index_col=0)
    donor_df['TYPE'] = donor_df['BEST'].apply(lambda x: x.split('-')[0])
    # select only Single Cells - exclude doublets and ambiguous cells ("DBL","AMB")
    print('{} out of {} cells are singlets.'.format((donor_df['TYPE']=='SNG').sum(), donor_df.shape[0]))
    donor_df = donor_df.query('TYPE=="SNG"')
        
    if 'donor_filter' in config_dict.keys():
        # filter out specific donors
        pool_id = row['pool_id']
        if pool_id in config_dict['donor_filter'].keys():
            # filter selected donors from this sample
            donors_to_filter_out = config_dict['donor_filter'][pool_id]
            donor_df = donor_df[~donor_df['SNG.1ST'].isin(donors_to_filter_out)]
    selected_cells = donor_df.index

    # read the count info, and merge with donor mapping
    data = sc.read_10x_h5(h5_file, genome='hg19')
    data = data[selected_cells, :]
    data.var_names_make_unique()
    #data.obs_names_make_unique()
    # label mitochondrial genes
    data.var['mito'] = [x.startswith('MT-') for x in data.var_names]
    # calculate qc metrics
    sc.pp.calculate_qc_metrics(data, qc_vars=['mito'], inplace=True)
    

    if 'cell_qc' in config_dict.keys():
        n_cells_start = data.n_obs
        for filter_query in config_dict['cell_qc']:
            data = data[data.obs.query(filter_query).index, :]
        print('cell QC applied: {} cells dropped'.format(n_cells_start-data.n_obs))
    

    # drop everything except gene ids and mito from data.var
    data.var = data.var[['gene_ids','mito']]
    if count>1:
        # drop everything else from data.var (duplicated with first dataset)
        data.var = data.var.iloc[:,:0]
    
    
    #apply gene list filtering if specified    
    if gene_list is not None:
        
        print(data)
        print(data.var)
        data = data[:, gene_list]
        
        
    data.obs['sample_id'] = sample_id
    data.obs['donor_id'] = donor_df.loc[data.obs.index, 'SNG.1ST']
    data.obs['pool_id'] = pool_id
    data.obs['time_point'] = row['time_point']
    data.obs['treatment'] = row['treatment']
    data.obs['tp_treat'] = data.obs['time_point'].astype(str)+"-"+data.obs['treatment'].astype(str)

    data.obs["sample_id"]=data.obs["sample_id"].astype("str")
    data.obs["donor_id"]=data.obs["donor_id"].astype("str")
    data.obs["pool_id"]=data.obs["pool_id"].astype("str")
    data.obs["time_point"]=data.obs["time_point"].astype("str")
    data.obs["treatment"]=data.obs["treatment"].astype("str")
    data.obs["tp_treat"]=data.obs["tp_treat"].astype("str")
    
    data.var_names=data.var_names.astype("str")
    data.obs_names=data.obs_names.astype("str")

    d= {'Sampleid' : sample_id,
        'pool_id' : pool_id,
        'time_point' : row['time_point'],
        'pathTofile' : h5_file,
        'initCells' : n_cells_start,
        'QCedCells' : data.n_obs
       }
    total_cells_counter.append(d)
    
    print('{} obs, {} vars'.format(data.n_obs, data.n_vars))
    
    if data.n_obs>0:
        print(data)
        print(data.obs)
        data.var_names_make_unique()
        datasets.append(data)
        count+=1
        
cols = data.obs.columns.tolist()
cols = cols[-5:] + cols[:-5]
data.obs = data.obs[cols]
data.obs

total_cells_counter = pd.DataFrame(total_cells_counter)

merged_data = datasets[0].concatenate(*datasets[1:])
print('{} obs, {} vars'.format(merged_data.n_obs, merged_data.n_vars))

cols = merged_data.obs.columns.tolist()
print(cols)
newcols = cols[:2]+cols[-6:-3]+[cols[-1]]+cols[2:-6]+cols[-3:-1]
print(newcols)
merged_data.obs = merged_data.obs[newcols]

print(merged_data)
print(merged_data.obs)
print(merged_data.var)
print('Write merged initial h5 file')

path_to_h5ad_file = 'output/allpools.scanpy.init.h5ad'
merged_data.write_h5ad(path_to_h5ad_file)

print('Write file with cell counts per 10xinlet')
output_file2="output/10xinlet_cellCounter.mitoSt5.tsv"
total_cells_counter.to_csv(output_file2,sep="\t", index=False)



