#!/bin/bash

OUTDIR=demuxlet/vcfGenerator
OUTFILE=${OUTDIR}/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.gnomadMAF0.05.biallelic.ALLchr.vcf.gz

## Exclude chromosomes X, Y and mithocondrial genome
OUTFILE2=${OUTDIR}/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20180102.genotypes.gnomadMAF0.05.biallelic.ALLchr.hg19.inputdemuxlet.vcf
zcat ${OUTDIR}/${OUTFILE} | grep -v -e '##contig=<ID=X' -e '##contig=<ID=Y' -e '##contig=<ID=MT' > ${OUTDIR}/${OUTFILE2}

