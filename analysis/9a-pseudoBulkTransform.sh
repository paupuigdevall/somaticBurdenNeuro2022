#!/bin/bash

echo $LSB_JOBINDEX
echo $1
input=$(head -n $LSB_JOBINDEX $1 | tail -n1)
tp=$(echo $input | awk '{print $1}')

echo $tp

Rscript 9a-pseudoBulkTransform.R ${tp}