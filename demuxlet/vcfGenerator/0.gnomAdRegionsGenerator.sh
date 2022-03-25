
#!/bin/bash

## Download gnomad.exomes.r2.1.1.sites.vcf.bgz from https://gnomad.broadinstitute.org/downloads
## Exomes > All chromosomes sites VCF 
## 58.81 GiB, MD5: f034173bf6e57fbb5e8ce680e95134f2

## We retain biallelic variants with a MAF frequency above 5% (subset of common biallelic variants)

INPUT=demuxlet/vcfGenerator

/software/bcftools view -q 0.05:minor -m2 -M2 -v snps ${INPUT}/gnomad.exomes.r2.1.1.sites.vcf.bgz -Ou | /software/bcftools-1.9/bcftools annotate -x ^INFO/AF -o ${INPUT}/gnomad.exomes.r2.1.1.sites.hg19.onlyINFOAF.aboveMAF05.biallelic.vcf.gz -O z

/software/htslib-1.9/tabix ${INPUT}/gnomad.exomes.r2.1.1.sites.hg19.onlyINFOAF.aboveMAF05.biallelic.vcf.gz
