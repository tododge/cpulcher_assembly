#!/bin/bash

############
# BAMTOOLS #
############

#bamtools version 2.5.1

ml biology bamtools

bamtools filter -in m64077_211121_052238.hifi_reads.bam -out m64077_211121_052238.hifi_reads.q20.bam -tag "rq":">=0.99"
bamtools filter -in m64204e_211127_230633.hifi_reads.bam -out m64204e_211127_230633.hifi_reads.q20.bam -tag "rq":">=0.99"


###################
# HIFIADAPTERFILT #
###################

#blastn version 2.16.0+
#hifiadapterfilt version 3.0.1

export PATH=$PATH:HiFiAdapterFilt/
export PATH=$PATH:HiFiAdapterFilt/DB

sh hifiadapterfilt.sh -t 16 -l 30
