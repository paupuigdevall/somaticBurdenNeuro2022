#!/bin/bash

PATH=/software/bcftools-1.9
PATH2=/software/htslib-1.9

OUTDIR=demuxlet/vcfGenerator
OUTFILE=hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.gnomadMAF0.05.biallelic.ALLchr.vcf.gz

## We concatenated all "chr" files, ordering them as required for demuxlet input
for chr in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9; do echo ${OUTDIR}/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.gnomadMAF0.05.biallelic.chr${chr}.vcf.gz >> sorted_list_of_chr_vcfs.txt ; done

$PATH/bcftools concat -o ${OUTDIR}/$OUTFILE -f sorted_list_of_chr_vcfs.txt -O z
$PATH2/tabix ${OUTDIR}/$OUTFILE

