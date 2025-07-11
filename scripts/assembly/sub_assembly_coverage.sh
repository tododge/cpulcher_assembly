#!/bin/bash

############
# BEDTOOLS #
############

ml biology bedtools

REF=crypul_v2024.1_rmY.fa

bedtools makewindows -w 100000 -g ${REF}.fai > ${REF}.windows_100000.bed

bedtools genomecov -ibam ${REF}.hifireads.sorted.q20.bam -bga > ${REF}.hifireads.sorted.q20.bam_whole-genome-coverage.bed

bedtools map -a ${REF}.windows_100000.bed -b ${REF}.hifireads.sorted.q20.bam_whole-genome-coverage.bed -c 4 -o mean > ${REF}.hifireads.sorted.q20.bam_whole-genome-coverage_windows_100000.bed
