#!/bin/bash

data=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $data
index0=$(echo $data | awk '{print $1}')
index1=$(echo $data | awk '{print $2}')
index2=$(echo $data | awk '{print $3}')

/lustre/scratch117/cellgen/kilpinen/shared/software/R-4.0.2/bin/Rscript corr_expression.R $index0 $index1 $index2

