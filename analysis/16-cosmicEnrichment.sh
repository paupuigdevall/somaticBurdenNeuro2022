echo $LSB_JOBINDEX
echo $1
input=$(head -n $LSB_JOBINDEX $1 | tail -n1)
tp=$(echo $input | awk '{print $1}')
cluster=$(echo $input | awk '{print $2}')

echo $tp
echo $cluster

Rscript 16-cosmicEnrichment.R ${tp} ${cluster}
