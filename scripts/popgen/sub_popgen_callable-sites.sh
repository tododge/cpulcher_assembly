#!/bin/bash

#############################
# BEDTOOLS CALLABLE REGIONS #
#############################

module load biology java gatk bedtools bcftools

#java version 21.0.4
#gatk version 4.6.0.0
#bedtools version 2.30.0
#bcftools version 1.16

bcftools query -f '%CHROM\t%POS0\t%POS\n' gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz | bedtools merge -i - > regions/crypul_v2024.1_rmY.fa.hifireads.called_regions.test.bed

awk '($3 - $2) > 1' regions/crypul_v2024.1_rmY.fa.hifireads.called_regions.test.bed > regions/crypul_v2024.1_rmY.fa.hifireads.called_regions.nosingle.test.bed

bcftools query -f '%CHROM\t%POS0\t%POS\n' gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz | bedtools merge -i - > regions/crypul_v2024.1_rmY.fa.omniCreads.called_regions.test.bed

awk '($3 - $2) > 1' regions/crypul_v2024.1_rmY.fa.omniCreads.called_regions.test.bed > regions/crypul_v2024.1_rmY.fa.omniCreads.called_regions.nosingle.test.bed

bedtools intersect -a regions/crypul_v2024.1_rmY.fa.hifireads.called_regions.nosingle.test.bed -b regions/crypul_v2024.1_rmY.fa.omniCreads.called_regions.nosingle.test.bed > regions/crypul_v2024.1_rmY.fa.hifireads.omniCreads.common.nosingle.called_regions.test.bed

for vcf in gatk_output/crypul_v2024.1_rmY.fa.*.sorted.q20.bam.auto*.SNP_INVARIANT.FILT.PASS.vcf.gz; do \
prefix=$(basename ${vcf} .SNP_INVARIANT.FILT.PASS.vcf.gz); \

bedtools intersect \
-a gatk_output/${prefix}.SNP_INVARIANT.FILT.PASS.vcf.gz \
-b regions/crypul_v2024.1_rmY.fa.hifireads.omniCreads.common.nosingle.called_regions.test.bed \
-header | bgzip > \
gatk_output/${prefix}.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.vcf.gz; \

gatk IndexFeatureFile \
    -I gatk_output/${prefix}.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.vcf.gz; done
