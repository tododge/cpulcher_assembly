#!/bin/bash

#####################
# GATK COMBINE VCFs #
#####################

module load biology java gatk

#java version 21.0.4
#gatk version 4.6.0.0

#hifi

gatk GatherVcfs \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr1.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr2.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr3.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr4.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr5.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr6.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr7.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr8.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr9.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr10.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr11.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr12.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr13.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr14.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr15.COMBINED.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.COMBINED.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.COMBINED.FILT.PASS.vcf.gz

gatk GatherVcfs \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr1.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr2.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr3.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr4.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr5.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr6.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr7.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr8.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr9.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr10.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr11.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr12.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr13.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr14.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.chr15.SNP_INVARIANT.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk GatherVcfs \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr1.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr2.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr3.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr4.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr5.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr6.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr8.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr9.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr10.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr11.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr12.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr13.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr14.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.chr15.SNP_INVARIANT.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz

#omniC
gatk GatherVcfs \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr1.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr2.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr3.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr4.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr5.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr6.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr7.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr8.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr9.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr10.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr11.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr12.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr13.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr14.COMBINED.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr15.COMBINED.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.COMBINED.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.COMBINED.FILT.PASS.vcf.gz

gatk GatherVcfs \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr1.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr2.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr3.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr4.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr5.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr6.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr7.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr8.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr9.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr10.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr11.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr12.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr13.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr14.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr15.SNP_INVARIANT.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk GatherVcfs \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr1.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr2.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr3.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr4.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr5.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr6.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr8.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr9.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr10.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr11.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr12.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr13.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr14.SNP_INVARIANT.FILT.PASS.vcf.gz \
-I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.chr15.SNP_INVARIANT.FILT.PASS.vcf.gz \
-O gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz

gatk IndexFeatureFile \
    -I gatk_output/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz
