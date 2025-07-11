#!/bin/bash

##################
# GATK FILT HIFI #
##################

module load biology java gatk

#java version 21.0.4
#gatk version 4.6.0.0

REF=crypul_v2024.1_rmY.fa
chr=$1

#mean depth of C. pulcher autosomes is 26.4. (1.5x=40, 0.66x=17) for hifi
#mean depth of C. pulcher autosomes is 73.2. (1.5x=110, 0.5x=37) for omniC
#mean depth of C. egeriae autosomes is 42.4. (1.5x=64, 0.5x=21) for hifi

#####
#adjust depth filters for different datasets
#####

#mean depth of C. pulcher autosomes is 26.4. (1.5x=40, 0.66x=17) for hifi
#mean depth of C. pulcher autosomes is 73.2. (1.5x=110, 0.5x=37) for omniC
#mean depth of C. egeriae autosomes is 42.4. (1.5x=64, 0.5x=21) for hifi


for vcf in gatk_output/${REF}.hifireads.sorted.q20.bam.${chr}.vcf.gz; do \
prefix=$(basename ${vcf} .vcf.gz); \
gatk SelectVariants \
-R reference/${REF} \
-V gatk_output/${prefix}.vcf.gz \
--select-type-to-include SNP \
-O gatk_output/${prefix}.SNP.vcf.gz; \

gatk SelectVariants \
-R reference/${REF} \
-V gatk_output/${prefix}.vcf.gz \
--select-type-to-include INDEL \
-O gatk_output/${prefix}.INDEL.vcf.gz; \

gatk SelectVariants \
-R reference/${REF} \
-V gatk_output/${prefix}.vcf.gz \
--select-type-to-include NO_VARIATION \
-O gatk_output/${prefix}.INVARIANT.vcf.gz; \

gatk VariantFiltration \
    -V gatk_output/${prefix}.SNP.vcf.gz \
    -filter "DP > 40" --filter-name "DP40" \
    -filter "DP < 10" --filter-name "DP10" \
    -filter "QD < 10.0" --filter-name "QD10" \
    -filter "QUAL < 60" --filter-name "QUAL60" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O gatk_output/${prefix}.SNP.FILT.vcf.gz; \

gatk VariantFiltration \
    -V gatk_output/${prefix}.INDEL.vcf.gz \
    -filter "DP > 40" --filter-name "DP40" \
    -filter "DP < 10" --filter-name "DP10" \
    -filter "QD < 10.0" --filter-name "QD10" \
    -filter "QUAL < 60.0" --filter-name "QUAL60" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O gatk_output/${prefix}.INDEL.FILT.vcf.gz; \

gatk VariantFiltration \
    -V gatk_output/${prefix}.INVARIANT.vcf.gz \
    -filter "DP > 40" --filter-name "DP40" \
    -filter "DP < 10" --filter-name "DP10" \
    --genotype-filter-expression "RGQ < 20" --genotype-filter-name "RGQ20" \
    -O gatk_output/${prefix}.INVARIANT.FILT.vcf.gz; \

zgrep -v RGQ20 gatk_output/${prefix}.INVARIANT.FILT.vcf.gz | \
awk '!($10 ~ /\.\/\.$|\.\/\.:[0-9]{1,3}:0/)' | bgzip > \
gatk_output/${prefix}.INVARIANT.FILT2.vcf.gz; \

gatk IndexFeatureFile \
    -I gatk_output/${prefix}.INVARIANT.FILT2.vcf.gz; \

gatk --java-options '-Xmx7g' MergeVcfs \
    -I gatk_output/${prefix}.SNP.FILT.vcf.gz \
    -I gatk_output/${prefix}.INDEL.FILT.vcf.gz \
    -I gatk_output/${prefix}.INVARIANT.FILT2.vcf.gz \
    -O gatk_output/${prefix}.COMBINED.FILT.vcf.gz; \

gatk --java-options '-Xmx7g' MergeVcfs \
    -I gatk_output/${prefix}.SNP.FILT.vcf.gz \
    -I gatk_output/${prefix}.INVARIANT.FILT2.vcf.gz \
    -O gatk_output/${prefix}.SNP_INVARIANT.FILT.vcf.gz; \

gatk SelectVariants \
    -R reference/${REF} \
    -V gatk_output/${prefix}.SNP_INVARIANT.FILT.vcf.gz \
    --exclude-filtered TRUE \
    -O gatk_output/${prefix}.SNP_INVARIANT.FILT.PASS.vcf.gz; \

gatk SelectVariants \
    -R reference/${REF} \
    -V gatk_output/${prefix}.COMBINED.FILT.vcf.gz \
    --exclude-filtered TRUE \
    -O gatk_output/${prefix}.COMBINED.FILT.PASS.vcf.gz; \

gatk SelectVariants \
    -R reference/${REF} \
    -V gatk_output/${prefix}.COMBINED.FILT.PASS.vcf.gz \
    --select-type-to-include SNP \
    -O gatk_output/${prefix}.SNP.FILT.PASS.vcf.gz; \

done
