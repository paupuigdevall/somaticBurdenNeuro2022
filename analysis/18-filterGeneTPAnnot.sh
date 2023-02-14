#!/bin/bash

data=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $data
index0=$(echo $data | awk '{print $1}') #timepoint
index1=$(echo $data | awk '{print $2}') #gene
index2=$(echo $data | awk '{print $3}') #ctype

Rscript 18-filterGeneTPAnnot.R $index0 $index1 $index2

