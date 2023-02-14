#!/bin/bash

data=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $data
index0=$(echo $data | awk '{print $1}')
index1=$(echo $data | awk '{print $2}')
index2=$(echo $data | awk '{print $3}')

Rscript 17-corr_expression.R $index0 $index1 $index2

