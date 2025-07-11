#!/bin/bash

#######################
# MINIMAP2 ALIGNMENTS #
#######################

#minimap2 version

#for plotting purposes, taking reverse complement of C. egeriae

seqkit grep -r -p "chr" cryege_v2024.1.fa | seqkit seq -r -p -t DNA -o cryege_v2024.1_rc.fa
minimap2 -x asm20 -t 16 cryege_v2024.1_rc.fa crypul_v2024.1.fa > cryege_v2024.1_rc.fa_2_crypul_v2024.1.fa.paf

######################
# MUMMER4 ALIGNMENTS #
######################

#mummer version 4.0.0rc1
ml biology viz gnuplot

seqkit grep -p chr7 cryege_v2024.1.fa -o cryege_v2024.1_chr7.fa
seqkit grep -p chr7 crypul_v2024.1.fa -o crypul_v2024.1_chr7.fa

seqkit grep -p ctg_008_Y_32451619_45003156 cryege_V1.0.SM_rmY.fa -o cryege_v2024.1_ctg_008_Y_32451619_45003156.fa
seqkit grep -p chrY crypul_v2024.1.fa -o crypul_v2024.1_chrY.fa

for ASSEM in $(ls crypul_v2024.1_chr7.fa); \
do SEQ_X=cryege_v2024.1_chr7.fa; \
PREFIX=$(echo ${ASSEM}\_2_${SEQ_X});
echo ${ASSEM}; \
nucmer -t 16 ${SEQ_X} ${ASSEM} --prefix=${PREFIX}; \
delta-filter -m ${PREFIX}.delta > ${PREFIX}.delta.m; \
show-coords ${PREFIX}.delta.m -T > ${PREFIX}.delta.m.coords; \
done

for ASSEM in $(ls crypul_v2024.1_chrY.fa); \
do SEQ_X=cryege_v2024.1_ctg_008_Y_32451619_45003156.fa; \
PREFIX=$(echo ${ASSEM}\_2_${SEQ_X});
echo ${ASSEM}; \
nucmer -t 16 ${SEQ_X} ${ASSEM} --prefix=${PREFIX}; \
delta-filter -m ${PREFIX}.delta > ${PREFIX}.delta.m; \
show-coords ${PREFIX}.delta.m -T > ${PREFIX}.delta.m.coords; \
done
