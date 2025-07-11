#!/bin/bash

#ragtag version 2.1.0

#ragtag scaffold without putative Y chromosome (chr16)
ragtag.py scaffold crypul_v2024.1_rmY.fa cryege_V1.0.SM_rmY.fa --mm2-params '-x asm20 -t 8' -o ./crypul_v2024.1_rmY.fa_ragtag_output

#unfortunately, ragtag made 2 misjoins, for ctg022 and ctg045. added these short scaffolds to two t2t contigs
#manually changed the agp

ragtag.py agp2fa ragtag.scaffold_corrected.agp ../cryege_V1.0.SM_rmY.fa | sed 's/_RagTag//g' | seqkit sort -l -r -o cryege_v2024.1_rmY.fa
