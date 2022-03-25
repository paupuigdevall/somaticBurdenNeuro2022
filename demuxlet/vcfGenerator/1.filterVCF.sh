#!/bin/bash

chr=$(head -n $LSB_JOBINDEX $1 | tail -n1)

## Path to HipSci VCF files (each VCF contains info for a given chromosome)
## VCF files are imputed and phased, and filtered by IMPUTE2 info score > 0.4
INDIR=pathToHipsci_gtarray_MergedFiles_REF-2018-01
OUTDIR=demuxlet/vcfGenerator

#Path_to_bcftools_installation
PATH=/software/bcftools-1.9

## We only keep the HipSci variants that overlap with the common biallelic variants from GATK Exomes (the regions file is the output from script "0.gnomAdRegionsGenerator.sh")

REGIONS=${OUTDIR}/gnomad.exomes.r2.1.1.sites.hg19.onlyINFOAF.aboveMAF05.biallelic.vcf.gz
INFILE=${INDIR}/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.${chr}.vcf.gz
OUTFILE=${OUTDIR}/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.gnomadMAF0.05.biallelic.${chr}.vcf.gz

$PATH/bcftools view -m2 -M2 -v snps -R $REGIONS $INFILE -Ou | $PATH/bcftools sort -o $OUTFILE -O z


