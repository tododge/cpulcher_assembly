#!/bin/bash

#########
# PLINK #
#########

ml biology plink

#plink version 1.90b6.5

plink --vcf crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
    --allow-extra-chr \
    --homozyg \
    --out crypul_v2024.1_rmY.fa.hifireads.default

plink --vcf crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
    --allow-extra-chr \
    --homozyg \
    --out crypul_v2024.1_rmY.fa.omniCreads.default
