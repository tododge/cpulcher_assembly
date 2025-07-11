#!/bin/bash

################
# DFAM MASKING #
################

singularity exec dfam-tetools-latest.sif BuildDatabase -name crypul_v2024.1.fa_repeat_db -engine ncbi crypul_v2024.1.fa

singularity exec dfam-tetools-latest.sif RepeatModeler -pa 62 -engine ncbi -database crypul_v2024.1.fa_repeat_db

REPEATLIB=./RM*/consensi.fa.classified #path to consensus file

singularity exec dfam-tetools-latest.sif RepeatMasker -pa 62 -gff -lib ${REPEATLIB} -xsmall crypul_v2024.1.fa -dir crypul_v2024.1.fa_RM_SM

cut -f 1,4,5 crypul_v2024.1.fa.out.gff | grep -v "#" > crypul_v2024.1.fa.out.bed
