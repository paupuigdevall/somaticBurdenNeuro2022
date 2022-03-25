#!/bin/bash

folder=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $folder

##need to download hg19 reference that comes with cellranger 3.1.0
TPATH=/projects/repo/refgenomes/cellranger/refdata-cellranger-hg19-3.0.0
sampleid=cellranger-hg19
cd $folder

/software/cellranger-3.1.0/cellranger count --id=${sampleid} --fastqs=${folder} --transcriptome=${TPATH} --jobmode=local --localcores=32 --localmem=200 --maxjobs=64

