#!/bin/bash

###########
# HIFIASM #
###########

#hifiasm version 0.16.1
hifiasm -o crypul-M.q20.default -t 48 \
m64077_211121_052238.hifi_reads.q20.filt.fq.gz m64204e_211127_230633.hifi_reads.q20.filt.fq.gz

#convert gfa to fa
for FILE in crypul-M.q20.default*p_ctg.gfa; do \
PREFIX=$(basename ${FILE} .gfa); \
echo ${PREFIX}; \
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.gfa > ${PREFIX}.fa; done


##############
# PURGE DUPS #
##############

#purge_dups version 0.0.3
ml system zlib

ls m64077_211121_052238.hifi_reads.q20.filt.fq.gz m64204e_211127_230633.hifi_reads.q20.filt.fq.gz > pbfofn

purge_dups/scripts/pd_config.py -n config.crypul.json crypul-M.q20.default.bp.p_ctg.fa pbfofn
purge_dups/scripts/run_purge_dups.py -p bash config.crypul.json purge_dups/bin/ crypul-M.q20.default_purged
