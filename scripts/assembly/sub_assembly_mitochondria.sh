#!/bin/bash

#seqkit version 2.10.0
#mitohifi version 2.2

seqkit seq m64077_211121_052238.hifi_reads.q20.filt.fq.gz m64204e_211127_230633.hifi_reads.q20.filt.fq.gz | seqkit fq2fa -o crypul-M.hifi_reads.filt.fa

singularity exec ../mitohifi_2.2_cv1.sif mitohifi.py -r crypul-M.hifi_reads.filt.fa -f BTS_final.fasta -g BTS_final.gb -t 6 -o 2

#after manual curration of mito contig, combined it with reference genome

seqkit seq crypul_V1.0.fa crypul_mito_currated | seqkit sort -l -r | sed 's/>chr16/>chrY/g' > crypul_v2024.1.fa
