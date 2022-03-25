#!/bin/bash

echo $LSB_JOBINDEX
echo $1
input=$(head -n $LSB_JOBINDEX $1 | tail -n1)
vcfFile=$(echo $input | awk '{print $1}')

echo $input
echo $vcfFile

Rscript 5-cellLine_varInfo.R ${vcfFile}

