#!/bin/bash

########
# PIXY #
########

#pixy version 1.2.10

zcat crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.vcf.gz | head -n 100 | grep CHROM | sed 's/\t/\n/g' | grep pacbio_libraryprep > pops.hifi

zcat crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.vcf.gz | head -n 100 | grep CHROM | sed 's/\t/\n/g' | grep crypul-M > pops.omni

#then edit pops files to add another column
#pacbio_libraryprep	crypul
#crypul-M	crypul

#allsites
pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
--populations pops.hifi \
--output_prefix crypul_v2024.1_rmY.fa.hifireads.10000 \
--window_size 10000 \
--n_cores 32

pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
--populations pops.hifi \
--output_prefix crypul_v2024.1_rmY.fa.hifireads.1000000 \
--window_size 1000000 \
--n_cores 32

pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
--populations pops.omni \
--output_prefix crypul_v2024.1_rmY.fa.omniCreads.10000 \
--window_size 10000 \
--n_cores 32

pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.vcf.gz \
--populations pops.omni \
--output_prefix crypul_v2024.1_rmY.fa.omniCreads.1000000 \
--window_size 1000000 \
--n_cores 32

#callable sites
pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.vcf.gz \
--populations pops.hifi \
--output_prefix crypul_v2024.1_rmY.fa.hifireads.callable.1000000 \
--window_size 1000000 \
--n_cores 32

pixy --stats pi \
--vcf crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.allchrs.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.vcf.gz \
--populations pops.omni \
--output_prefix crypul_v2024.1_rmY.fa.omniCreads.callable.1000000 \
--window_size 1000000 \
--n_cores 32
