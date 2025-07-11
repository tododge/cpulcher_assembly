#!/bin/bash

#############################
# BEDTOOLS CALLABLE REGIONS #
#############################

ml biology bcftools viz gnuplot

#bcftools version 1.16
#gnuplot version 5.2.0

bcftools view crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz | awk -F"\t" '{ if ($0 ~ /^#/) {print $0} else if ($5 == ".") {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 ";FQ=-40""\t" $9 "\t" $10} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 ";FQ=40""\t" $9 "\t" $10}}' | vcfutils.pl vcf2fq | sed 's/[a,t,c,g]/\U&/g' | bgzip > input/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.fastq.gz

./psmc/utils/fq2psmcfa -q 20 input/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.fastq.gz > input/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa

#psmc all sites
./psmc/psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o ./results/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmc ./input/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa

./psmc/utils/splitfa ./input/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa > ./input/crypul_v2024_rmY.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS-split.psmcfa

#hifi & omni-C callable sites
./psmc/psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o ./results/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.psmc ./input/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.psmcfa

bcftools view crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.vcf.gz | awk -F"\t" '{ if ($0 ~ /^#/) {print $0} else if ($5 == ".") {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 ";FQ=-40""\t" $9 "\t" $10} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 ";FQ=40""\t" $9 "\t" $10}}' | vcfutils.pl vcf2fq | sed 's/[a,t,c,g]/\U&/g' | bgzip > input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.fastq.gz

./psmc/utils/fq2psmcfa -q 20 input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.fastq.gz > input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa

#psmc all sites
./psmc/psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o ./results/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmc ./input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa

./psmc/utils/splitfa ./input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.psmcfa > ./input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS-split.psmcfa

#hifi & omni-C callable sites
./psmc/psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o ./results/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.psmc ./input/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.psmcfa

for file in $(ls results/*.psmc); do \
PREFIX=$(basename ${file} .psmc); \
./psmc/utils/psmc2history.pl ./results/${PREFIX}.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh; \
./psmc/utils/psmc_plot.pl -u 6.25e-09 -g 3.5 -R ./results/${PREFIX} ./results/${PREFIX}.psmc; done

for file in $(ls results/bootstrap/*psmc); do \
PREFIX=$(basename ${file} .psmc); \
./psmc/utils/psmc2history.pl ./results/bootstrap/${PREFIX}.psmc | ./psmc/utils/history2ms.pl > ms-cmd.sh; \
./psmc/utils/psmc_plot.pl -u 6.25e-09 -g 3.5 -R ./results/bootstrap/${PREFIX} ./results/bootstrap/${PREFIX}.psmc; done


awk '{print $0, "\t", "C. pulcher","\t", "hifi autosomes", "\t", "primary"}' results/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.0.txt > Cpulcher_hifi_nonbootstrap.txt
for i in $(ls results/bootstrap/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.round*txt); do awk '{print $0, "\t", "C. pulcher bootstrap", "\t", "hifi autosomes", "\t", FILENAME}' $i; done > Cpulcher_hifi_bootstrap.txt

awk '{print $0, "\t", "C. pulcher","\t", "omniC autosomes", "\t", "primary"}' results/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.0.txt > Cpulcher_omniC_nonbootstrap.txt
for i in $(ls results/bootstrap/crypul_v2024.1_rmY.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.round*txt); do awk '{print $0, "\t", "C. pulcher bootstrap", "\t", "omniC autosomes", "\t", FILENAME}' $i; done > Cpulcher_omniC_bootstrap.txt

cat *bootstrap* | sed -i 's/ //' > pulcher_comb.txt
