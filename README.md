[![DOI](https://zenodo.org/badge/474143031.svg)](https://zenodo.org/badge/latestdoi/474143031)


# Somatic mutations alter the differentiation outcomes of iPSC-derived neurons

This repository contains the scripts to reproduce the analysis and the figures of this paper (DOI): https://doi.org/10.1101/2022.03.04.482992

It is divided in the following sections:
- **cellranger**: Processing of 119 10x samples (libraries) from Jerber et al.2021
- **demuxlet**: It deconvolutes the donor origin of pooled lines by generating a VCF with the common biallelic variants (gnomad, >5%) from HipSci lines (genotype data).
- **singleCellProcessing**: Using the Scanpy platform (Python 3.7.3), we merged the data, performed dimensionality reduction, batch correction, clustering and annotate them. h5 objects are deposited in https://zenodo.org/record/6079122#.Yj4p7S2l1QI
- **analysis**: It includes the scripts that analyse the impact of somatic variants in the pooled dopaminergic differentiation profiled at days 11, 30 and 52. It generates all the intermediate tables necessary to reproduce the main and supplementary figures. We included the input files required to run the array jobs.

In the main directory, there are also those scripts:
- **commands.sh**: indicates the order to run the entire pipeline and be able to reproduce the main figures (figure*.R) and the supplementary figures (suppFigure*.R) of the paper.
- **libraryInstall.R**: contains the list of required R packages to run the data analysis and to reproduce the figures of the paper.
- **functionsToImport.R**: list of functions to be imported to reproduce figures.
- **suppTabs/**: Folder with the supplementary tables (tab-separated values).
