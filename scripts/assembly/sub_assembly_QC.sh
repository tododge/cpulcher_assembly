#!/bin/bash

######################
# BUSCO COMPLETENESS #
######################

#busco version 5.3.2; vertebrata_odb10

singularity exec busco_v5.3.2_cv1.sif busco -i ../final_genome/crypul_v2024.1.fa -l vertebrata_odb10 -o crypul_v2024.1.fa_busco_vertebrata_odb10_26X24 -m geno -c 32 --long

singularity exec busco_v5.3.2_cv1.sif busco -i ../final_genome/cryege_v2024.1.fa -l vertebrata_odb10 -o cryege_v2024.1.fa_busco_vertebrata_odb10_26X24 -m geno -c 32 --long


###########################
# SEQTK SEQKIT CONTIGUITY #
###########################

for $genome in $(ls cryege_v2024.1.fa crypul_v2024.1.fa); do \
seqtk telo $genome > $genome_telo.bed; \
seqtk gap $genome > $genome_gap.bed; \
seqkit stat -a -T $genome > $genome.stats.txt; done
