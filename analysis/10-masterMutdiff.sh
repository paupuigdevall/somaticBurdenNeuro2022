#!/bin/bash

echo $LSB_JOBINDEX
echo $1
input=$(head -n $LSB_JOBINDEX $1 | tail -n1)
failed=$(echo $input | awk '{print $1}')
succ=$(echo $input | awk '{print $2}')

echo $input
echo $failed
echo $succ

Rscript 10-masterMutdiff.R ${failed} ${succ}
