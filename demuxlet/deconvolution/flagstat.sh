#!/bin/bash

echo $LSB_JOBINDEX
echo $1
data=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $data
file=$(echo $data | awk '{print $1}')
echo $file

summary=bam_flagstatQC.txt
file_name=$(echo $file | sed 's/\.bam//g')
tmp_file=tmp${file_name}.test
tmp_file2=$(echo $tmp_file | sed 's,\/,_,g')

/lustre/scratch117/cellgen/kilpinen/shared/software/samtools-1.9/bin/samtools flagstat $file > log/$tmp_file2 
echo $file_name >> $summary
grep "%" log/$tmp_file2 | cut -d'(' -f2 | cut -d':' -f1 >> $summary  

